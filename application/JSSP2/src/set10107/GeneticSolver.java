package set10107;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import modelP.JSSP;
import modelP.Problem;

public class GeneticSolver {
	
	private Problem problem;

	public GeneticSolver(){
		this.problem = JSSP.getProblem(90);
	}
	
	public GeneticSolver(Problem p){
		this.problem = p;
	}
	
	public Schedule solve_elitist(Problem p, int nb_generations, int population_size, double mutation_rate, int crossover_type){
		// return variable
		Schedule best_schedule = new Schedule(p);
		this.problem = p;
		ArrayList<Schedule> population = this.generate_population(population_size);
		
		String header = "Optimization method, Problem id, Generation, GA_elitist - best indiv fitness, Lower Bound";
		System.out.println(header);
		
		for(int i=0; i<nb_generations; ++i) {
			evaluate(population);

			String line = "GA_elitist, " + problem.getId() + ", " + i + ", " + population.get(0).get_fitness(this.problem) + ", " +problem.getLowerBound();
			System.out.println(line);

			Schedule offspring=null;
			switch(crossover_type){
			case 0:
				offspring = one_point_crossover(new Schedule(population.get(0), this.problem), new Schedule (population.get(1), this.problem));
				break;
			case 1:
				offspring = two_point_crossover(new Schedule(population.get(0), this.problem), new Schedule (population.get(1), this.problem));
				break;
			case 2:
				offspring = three_parents_crossover(population.get(0), population.get(1), population.get(2));
				break;
			default:
				offspring = three_parents_crossover(population.get(0), population.get(1), population.get(2));
			}

			Schedule mutated_offspring = mutate(new Schedule(offspring, this.problem), mutation_rate);
			population = replace(population, mutated_offspring);
		}

		evaluate(population);
		best_schedule = population.get(0);
		return best_schedule;
	}
	
	
	public Schedule solve_hybrid_elitist_neighborhood(Problem p, int nb_generations, int population_size, double mutation_rate, int crossover_type, Double[][] avrg_results, int replicate_index){
		// return variable
		Schedule best_schedule = new Schedule(p);
		this.problem = p;
		ArrayList<Schedule> population = this.generate_population(population_size);
		
		System.out.println("NEW HYBRID TRAINING");
		for(int i=0; i<nb_generations; ++i) {
			evaluate(population);
			
			//
			// PERFORM LOCAL SEARCH
			//
			Schedule local_search_res = neighborhoods_search(population.get(0), nb_generations);
			if(local_search_res.get_fitness(problem) < population.get(population_size-1).get_fitness(problem)){
				population.set(population_size-1, local_search_res);
			}
			evaluate(population);
			
			if(population.get(0).get_fitness(problem) <= best_schedule.get_fitness(problem))
				best_schedule = population.get(0);

			// aggregate fitness values in vector
			double[] population_scores = new double[population.size()];
			for(int s=0; s<population.size(); s++)
				population_scores[s] = population.get(s).get_fitness(problem);
			
			int res_index = i ;//= i+(nb_generations*(p.getId()-1)) - nb_generations;
			// record performance metrics for this replicate
			avrg_results[res_index][0] += (double) problem.getId();
			avrg_results[res_index][1] += (double) problem.getNumberOfJobs();
			avrg_results[res_index][2] += (double) problem.getNumberOfMachines();
			avrg_results[res_index][3] += (double) 0; // OPTIMIZATION ALGORITHM: GA
			avrg_results[res_index][4] += (double) i;
			avrg_results[res_index][5] += best_schedule.get_fitness(this.problem);
			avrg_results[res_index][6] += (double) problem.getLowerBound();
			avrg_results[res_index][7] += JsspExperiment.compute_mean(population_scores, population.size());
			avrg_results[res_index][8] += JsspExperiment.compute_variance(population_scores, population_size);
			
			Schedule offspring=null;
			switch(crossover_type){
			case 0:
				offspring = one_point_crossover(new Schedule(population.get(0), this.problem), new Schedule (population.get(1), this.problem));
				break;
			case 1:
				offspring = two_point_crossover(new Schedule(population.get(0), this.problem), new Schedule (population.get(1), this.problem));
				break;
			case 2:
				offspring = three_parents_crossover(population.get(0), population.get(1), population.get(2));
				break;
			default:
				offspring = three_parents_crossover(population.get(0), population.get(1), population.get(2));
			}
			
			Schedule mutated_offspring = mutate(new Schedule(offspring, this.problem), mutation_rate);
			population = replace(population, mutated_offspring);
			
			String line = "";
			// record performance metrics for current replicate
			line = line + "PB.ID=" + problem.getId();
			line = line + "\t, " + " PB.NB.JOBS=" + problem.getNumberOfJobs();
			line = line + "\t, " + " PB.NB.MACHINES=" + problem.getNumberOfMachines();
			line = line + "\t, " + " 1=" + 0; // OPTIMIZATION ALGORITHM: AIS
			line = line + "\t, " + " GEN=" + i;
			line = line + "\t, " + " FITNESS=" + best_schedule.get_fitness(this.problem);
			line = line + "\t, " + " LB=" + problem.getLowerBound();
			line = line + "\t, " + " MEAN.FIT=" + JsspExperiment.compute_mean(population_scores, population.size());
			line = line + "\t, " + " VAR.FIT=" + JsspExperiment.compute_variance(population_scores, population.size());
			System.out.println(line);
		}
		evaluate(population);
		best_schedule = population.get(0);
		return best_schedule;
	}

	public Schedule solve_tournament(Problem p, int nb_generations, int population_size, int tournament_size, int selection_pressure, double mutation_rate, int crossover_type, Double[][] avrg_results, int replicate_index){
		// return variable
		Schedule best_schedule = new Schedule(p);
		this.problem = p;
		ArrayList<Schedule> population 		= this.generate_population(population_size);
		ArrayList<Schedule> tournament_pop 	= new ArrayList<Schedule>();
		ArrayList<Schedule> parents 		= new ArrayList<Schedule>();
		
		String header = "Optimization method, Problem id, Generation, GA_tnmt - best indiv fitness, Lower Bound";
		System.out.println(header);
		
		for(int i=0; i<nb_generations; ++i) {
			//
			// PERFORM LOCAL SEARCH
			//
			Schedule local_search_res = neighborhoods_search(population.get(0), nb_generations);
			if(local_search_res.get_fitness(problem) < population.get(population_size-1).get_fitness(problem)){
				population.set(population_size-1, local_search_res);
			}
			evaluate(population);

			if(population.get(0).get_fitness(problem) <= best_schedule.get_fitness(problem))
				best_schedule = population.get(0);
			// aggregate fitness values in vector
			double[] population_scores = new double[population.size()];
			for(int s=0; s<population.size(); s++)
				population_scores[s] = population.get(s).get_fitness(problem);
			
			int res_index = i ;//= i+(nb_generations*(p.getId()-1)) - nb_generations;
			// record performance metrics for this replicate
			avrg_results[res_index][0] += (double) problem.getId();
			avrg_results[res_index][1] += (double) problem.getNumberOfJobs();
			avrg_results[res_index][2] += (double) problem.getNumberOfMachines();
			avrg_results[res_index][3] += (double) 0; // OPTIMIZATION ALGORITHM: GA
			avrg_results[res_index][4] += (double) i;
			avrg_results[res_index][5] += best_schedule.get_fitness(this.problem);
			avrg_results[res_index][6] += (double) problem.getLowerBound();
			avrg_results[res_index][7] += JsspExperiment.compute_mean(population_scores, population.size());
			avrg_results[res_index][8] += JsspExperiment.compute_variance(population_scores, population_size);
			
			tournament_pop 	= generate_tournament_population(tournament_size, population);
			parents			= select_tournament_parents(tournament_pop, selection_pressure);
			
			for(int k=0; k< parents.size(); ++k) {
				Schedule offspring=null;
				switch(crossover_type){
				case 0:
					offspring = one_point_crossover(new Schedule(population.get(0), this.problem), new Schedule (population.get(1), this.problem));
					break;
				case 1:
					offspring = two_point_crossover(new Schedule(population.get(0), this.problem), new Schedule (population.get(1), this.problem));
					break;
				case 2:
					offspring = three_parents_crossover(population.get(0), population.get(1), population.get(2));
					break;
				default:
					offspring = three_parents_crossover(population.get(0), population.get(1), population.get(2));
				}
				
				Schedule mutated_offspring = mutate(new Schedule(offspring, this.problem), mutation_rate);
				population = replace(population, mutated_offspring);
			}
			String line = "";
			// record performance metrics for current replicate
			line = line + "PB.ID=" + problem.getId();
			line = line + "\t, " + " PB.NB.JOBS=" + problem.getNumberOfJobs();
			line = line + "\t, " + " PB.NB.MACHINES=" + problem.getNumberOfMachines();
			line = line + "\t, " + " 1=" + 0; // OPTIMIZATION ALGORITHM: AIS
			line = line + "\t, " + " GEN=" + i;
			line = line + "\t, " + " FITNESS=" + best_schedule.get_fitness(this.problem);
			line = line + "\t, " + " LB=" + problem.getLowerBound();
			line = line + "\t, " + " MEAN.FIT=" + JsspExperiment.compute_mean(population_scores, population.size());
			line = line + "\t, " + " VAR.FIT=" + JsspExperiment.compute_variance(population_scores, population.size());
			System.out.println(line);
		}
		evaluate(population);
		best_schedule = population.get(0);
		return best_schedule;
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
	
	public ArrayList<Schedule> generate_population(int pop_size){
		ArrayList<Schedule> pop = new ArrayList<Schedule>();
		for(int i=0 ; i<pop_size ;++i){
			pop.add(new Schedule(this.problem));
		}
		return pop;
	}
	
	public void evaluate(ArrayList<Schedule> population){
		for(int i=0; i<population.size() ; ++i)
			population.get(i).update_fitness(this.problem);
		ScheduleComparator c = new ScheduleComparator(this.problem);
		Collections.sort(population, c);
	}
	
	public Schedule one_point_crossover(Schedule s1, Schedule s2){
		// return variable
		Schedule offspring = new Schedule(this.problem);
		Schedule s1_copy = new Schedule(s1, problem);
		Schedule s2_copy = new Schedule(s2, problem);
		
		int genotype_height = s1_copy.get_genotype().length;
		int genotype_width  = s1_copy.get_genotype()[0].length;
		int[][] offspring_genotype = new int[genotype_height][genotype_width];
		
		for(int i=0; i<genotype_height; ++i){
			if(i < genotype_height/2){
				offspring_genotype[i] = s1_copy.get_genotype()[i];
			}else{
				offspring_genotype[i] = s2_copy.get_genotype()[i];
			}
		}
		offspring.set_genotype(offspring_genotype);
		return offspring;
	}
	
	public Schedule two_point_crossover(Schedule s1, Schedule s2){
		// return variable
		Schedule offspring = new Schedule(this.problem);
		int genotype_height = s1.get_genotype().length;
		int genotype_width  = s1.get_genotype()[0].length;
		int[][] offspring_genotype = new int[genotype_height][genotype_width];
		
		for(int i=0; i<genotype_height; ++i) {
			if(i < genotype_height/3){
				offspring_genotype[i] = s1.get_genotype()[i];
			}else if(i < (genotype_height*(2/3))){
				offspring_genotype[i] = s2.get_genotype()[i];
			}else{
				offspring_genotype[i] = s1.get_genotype()[i];
			}
		}
		offspring.set_genotype(offspring_genotype);
		return offspring;
	}
	
	public Schedule three_parents_crossover(Schedule s1, Schedule s2, Schedule s3){
		// return variable
		Schedule offspring = new Schedule(this.problem);
		int genotype_height = s1.get_genotype().length;
		int genotype_width  = s1.get_genotype()[0].length;
		int[][] offspring_genotype = new int[genotype_height][genotype_width];
		
		for(int i=0; i<genotype_height; ++i){
			if(i < genotype_height/3){
				offspring_genotype[i] = s1.get_genotype()[i];
			}else if(i < (genotype_height*(2/3))){
				offspring_genotype[i] = s2.get_genotype()[i];
			}else{
				offspring_genotype[i] = s3.get_genotype()[i];
			}
		}
		offspring.set_genotype(offspring_genotype);
		return offspring;
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
	
	public ArrayList<Schedule> replace(ArrayList<Schedule> pop, Schedule offspring) {
		evaluate(pop);
		double current_worst_score = pop.get(pop.size()-1).get_fitness(problem);
		double offspring_score = offspring.get_fitness(problem);
		if(offspring_score < current_worst_score) {
			Schedule offspring_copy = new Schedule(offspring, this.problem);
			// replace
			pop.set(pop.size()-1, offspring_copy);
		}
		return pop;
	}
	
	public ArrayList<Schedule> generate_tournament_population(int tournament_size, ArrayList<Schedule> pop){
		// return variable
		ArrayList<Schedule> tournament_pop = new ArrayList<Schedule>();
		
		// select k individuals at random from pop
		for(int k=0; k<tournament_size; ++k) {
			Random rand = new Random();
			// copy obj
			Schedule s = new Schedule(this.problem);
			s.set_genotype(pop.get(rand.nextInt(pop.size())).get_genotype());
			tournament_pop.add(s);
		}
		return tournament_pop;
	}

	private ArrayList<Schedule> select_tournament_parents(ArrayList<Schedule> tournament_pop, int selection_pressure) {
		ArrayList<Schedule> parents = new ArrayList<Schedule>();
		evaluate(tournament_pop);
		for(int i=0; i<selection_pressure; i++){
			parents.add(tournament_pop.get(i));
		}
		return parents;
	}
}
