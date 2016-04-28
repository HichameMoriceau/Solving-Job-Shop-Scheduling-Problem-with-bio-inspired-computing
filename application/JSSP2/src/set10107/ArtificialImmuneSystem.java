package set10107;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import modelP.JSSP;
import modelP.Problem;

public class ArtificialImmuneSystem {

	private Problem problem;
	
	public ArtificialImmuneSystem(){
		// default problem
		problem = JSSP.getProblem(90);
	}

	public ArrayList<Schedule> generate_population(int pop_size){
		ArrayList<Schedule> pop = new ArrayList<Schedule>();
		for(int i=0 ; i<pop_size ;++i){
			pop.add(new Schedule(this.problem));
		}
		return pop;
	}

	public void evaluate_affinity(ArrayList<Schedule> population){
		for(int i=0; i<population.size() ; ++i)
			population.get(i).update_fitness(this.problem);
		ScheduleComparator c = new ScheduleComparator(this.problem);
		Collections.sort(population, c);
	}
	
	private ArrayList<Schedule> generate_clones(Schedule s, double nb_clones) {
		ArrayList<Schedule> clone_pop = new ArrayList<Schedule>();
		for(int i=0; i<nb_clones; i++){
			clone_pop.add(new Schedule(s, problem));
		}
		return clone_pop;
	}

	// selects the <n> highest antibodies from affinity-ordered <population>
	private ArrayList<Schedule> select_antibodies(ArrayList<Schedule> population, int n) {
		// evaluate affinity and rank population
		evaluate_affinity(population);
		ArrayList<Schedule> antibodies = new ArrayList<Schedule>();
		for(int i=0; i<n; i++)
			antibodies.add(population.get(i));
		return antibodies;
	}
	
	// selects the <n> highest antibodies from affinity-ordered <population>
		private ArrayList<Schedule> select_antibodies(ArrayList<Schedule> population, ArrayList<Schedule> clones, int n) {
			ArrayList<Schedule> full_pop = new ArrayList<Schedule>();
			full_pop.addAll(population);
			full_pop.addAll(clones);
			evaluate_affinity(full_pop);
			ArrayList<Schedule> antibodies = new ArrayList<Schedule>();
			for(int i=0; i<n; i++)
				antibodies.add(population.get(i));
			return antibodies;
		}
	
	public Schedule mutate(Schedule s, double mutation_rate){
		// return variable
		Schedule mutated_schedule = new Schedule(s,this.problem);
		Random r = new Random();
		int genotype_height = mutated_schedule.get_genotype().length;
		int genotype_width = mutated_schedule.get_genotype()[0].length;
		if(r.nextDouble()<mutation_rate){
			int[][] mutated_genome = mutated_schedule.get_genotype();
			
			// for each machine
			for(int i=0; i<genotype_height; ++i){
				int random_genome_index_A = -1;
				int random_genome_index_B = -1;
				// select two different jobs
				do {
					random_genome_index_A = r.nextInt(mutated_schedule.get_genotype().length);
					random_genome_index_B = r.nextInt(mutated_schedule.get_genotype().length);
				}while(random_genome_index_A == random_genome_index_B);
				// swap selected jobs
				int tmp = mutated_genome[i][random_genome_index_A];
				mutated_genome[i][random_genome_index_A] = mutated_schedule.get_genotype()[i][random_genome_index_B];
				mutated_genome[i][random_genome_index_B] = tmp;
			}
			mutated_schedule.set_genotype(mutated_genome);	
			return mutated_schedule;
		}
		return mutated_schedule;
	}
	

	private double compute_nb_clones(double beta, int pop_size, int index) {
		return Math.round((beta*pop_size)/index);
	}
	
	public ArrayList<Schedule> replace(ArrayList<Schedule> pop, Schedule offspring) {
		evaluate_affinity(pop);
		double current_worst_score = pop.get(pop.size()-1).get_fitness(problem);
		double offspring_score = offspring.get_fitness(problem);
		if(offspring_score < current_worst_score) {
			Schedule offspring_copy = new Schedule(offspring, this.problem);
			// replace
			pop.set(pop.size()-1, offspring_copy);
		}
		return pop;
	}
	
	public Schedule solve(Problem p, int nb_epochs, int population_size, Double[][] avrg_results, int replicate_index) {
		problem = p;
		Schedule optimized_schedule = new Schedule(problem);		
		// AIS settings
		int selection_size 				= (int) population_size*100/100;
		double scaling_factor			= 5;
		double mutation_decay_factor 	= 0.7;
		int nb_random_cells 			= (int) population_size*10/100;

		// create population of antibodies and antigens
		ArrayList<Schedule> population = generate_population(population_size);
		System.out.println("NEW AIS TRAINING");
		for(int i=0; i<nb_epochs; i++){
			evaluate_affinity(population);
			// aggregate fitness values in vector
			double[] population_scores = new double[population.size()];
			for(int s=0; s<population.size(); s++)
				population_scores[s] = population.get(s).get_fitness(problem);
			
			int res_index = i;// = i+(nb_epochs*(p.getId()-1));
			// record performance metrics for current replicate
			avrg_results[res_index][0] += (double) problem.getId();
			avrg_results[res_index][1] += (double) problem.getNumberOfJobs();
			avrg_results[res_index][2] += (double) problem.getNumberOfMachines();
			avrg_results[res_index][3] += (double) 1; // OPTIMIZATION ALGORITHM: AIS
			avrg_results[res_index][4] += (double) i;
			avrg_results[res_index][5] += population.get(0).get_fitness(this.problem);
			avrg_results[res_index][6] += (double) problem.getLowerBound();
			avrg_results[res_index][7] += ((Double)JsspExperiment.compute_mean(population_scores, population.size())).intValue();
			avrg_results[res_index][8] += ((Double)JsspExperiment.compute_variance(population_scores, population_scores.length)).intValue();
			
			Schedule local_search_res = neighborhoods_search(population.get(0), nb_epochs);
			if(local_search_res.get_fitness(problem) < population.get(population.size()-1).get_fitness(problem)){
				population.set(population.size()-1, local_search_res);
			}
			
			// instantiate sub-populations
			ArrayList<Schedule> antibodies_pop = select_antibodies(population, selection_size);
			ArrayList<ArrayList<Schedule>> clones_pop = new ArrayList<ArrayList<Schedule>>();
			
			// generate clone population
			for(int s=0; s<antibodies_pop.size(); s++){
				double nb_clones = compute_nb_clones(scaling_factor, antibodies_pop.size(), s+1);
				//System.out.println("Cloning " + nb_clones + " antibodies for individual of rank " + s);
				clones_pop.add( generate_clones(antibodies_pop.get(s), nb_clones));
			}
			
			// hypermutate clones
			for(ArrayList<Schedule> clones_subpopulation : clones_pop){
				for(Schedule s : clones_subpopulation){
					// best antibodies mutate at a lower degree
					double mutation_rate = compute_mutation_rate(mutation_decay_factor, s.get_affinity(problem));
					s = mutate(s, mutation_rate);
					// update fitness/affinity
					s.get_affinity(problem);
				}
			}

			ArrayList<Schedule> full_clones_pop = new ArrayList<Schedule>();
			for(ArrayList<Schedule> clones : clones_pop)
				full_clones_pop.addAll(clones);
			// re-select best n cells
			population = select_antibodies(population, full_clones_pop, selection_size);
			
			// introduce random cells
			ArrayList<Schedule> random_pop = generate_population(nb_random_cells);
			for(int r_idx=0; r_idx<nb_random_cells; r_idx++)
				population.set(population.size()-1-r_idx, random_pop.get(r_idx));
			
			population = select_antibodies(population, selection_size);
			
			String line = "";
			// record performance metrics for current replicate
			line = line + "PB.ID=" + problem.getId();
			line = line + "\t, " + " PB.NB.JOBS=" + problem.getNumberOfJobs();
			line = line + "\t, " + " PB.NB.MACHINES=" + problem.getNumberOfMachines();
			line = line + "\t, " + " 1=" + 1; // OPTIMIZATION ALGORITHM: AIS
			line = line + "\t, " + " GEN=" + i;
			line = line + "\t, " + " FITNESS=" + population.get(0).get_fitness(this.problem);
			line = line + "\t, " + " LB=" + problem.getLowerBound();
			line = line + "\t, " + " MEAN.FIT=" + JsspExperiment.compute_mean(population_scores, population.size());
			line = line + "\t, " + " VAR.FIT=" + JsspExperiment.compute_variance(population_scores, population.size());
			System.out.println(line);
		}
		optimized_schedule = population.get(0);
		return optimized_schedule;
	}

	private double compute_mutation_rate(double decay_factor, double affinity) {
		return Math.exp((-decay_factor)*affinity);
	}
	
	private Schedule neighborhoods_search(Schedule s, int nb_epochs){
		// return variable
		Schedule best_model  		= new Schedule(s, problem);
		Schedule mutated_schedule 	= new Schedule(s, problem);
		//System.out.println("Old solution:" + mutated_schedule.get_fitness(problem));
		
		Random r = new Random();
		int genotype_height = mutated_schedule.get_genotype().length;
		int[][] mutated_genome = mutated_schedule.get_genotype();
		
		for(int i=0; i<nb_epochs; i++){
			// for a random machine
			int random_machine_index  = r.nextInt(genotype_height);
			int random_genome_index_A = -1;
			int random_genome_index_B = -1;
			// select a job at random
			random_genome_index_A = r.nextInt(mutated_schedule.get_genotype().length-1);
			// select index of right neighbour
			random_genome_index_B = (random_genome_index_A+1);
			
			//System.out.println("job indexes to be swapped: " + random_genome_index_A + " and " + random_genome_index_B + " for machine " + random_machine_index);
			// swap jobs
			int tmp = mutated_genome[random_machine_index][random_genome_index_A];
			mutated_genome[random_machine_index][random_genome_index_A] = mutated_schedule.get_genotype()[random_machine_index][random_genome_index_B];
			mutated_genome[random_machine_index][random_genome_index_B] = tmp;

			mutated_schedule.set_genotype(mutated_genome);
			if(mutated_schedule.get_fitness(problem) < best_model.get_fitness(problem)){
				best_model = new Schedule(mutated_schedule, problem);
			}
		}
		return mutated_schedule;
	}
}
