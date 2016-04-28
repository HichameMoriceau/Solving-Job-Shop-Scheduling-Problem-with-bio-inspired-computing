package set10107;

import java.io.File;
import java.io.PrintWriter;

import modelP.JSSP;
import modelP.Problem;

public class JsspExperiment {
	
	private static void save_performances_on_all_problems(){
		
		// average results over <n> replicates
		GeneticSolver GA = new GeneticSolver();
		ArtificialImmuneSystem AIS = new ArtificialImmuneSystem();
		
		// --------------------------------------------------------
		//					EXPERIMENT SETTINGS
		//--------------------------------------------------------
		
		int nb_replicates 				= 		10;
		int nb_generations 				= 		1;
		int population_size 			= 		10;
		double mutation_rate 			= 		0.65f;
		int tournament_size 			= 		population_size;
		int selection_pressure 			= 		tournament_size;
		int nb_recorded_metrics 		= 		9;
		int nb_problems 				= 		10;
		int start_pb 					= 		1;
		// --------------------------------------------------------
		
		PrintWriter writer = null;
		try {
			writer = new PrintWriter("results.csv", "UTF-8");
			writer.println("problem id, NB jobs, NB machines, OptAlgorithm(0=GA-1=AIS), generation, fitness, Lower Bound, mean fitness pop, variance fitness pop");
			
			int nb_line_res_file = nb_generations*nb_problems;
			// result matrices
			Double[][] avrg_results_GA = new Double[nb_line_res_file][nb_recorded_metrics];
			Double[][] avrg_results_AIS = new Double[nb_line_res_file][nb_recorded_metrics];
			for(int i=0; i<nb_line_res_file; i++){
				for(int j=0; j<nb_recorded_metrics; j++){
					avrg_results_GA[i][j]  = (double) 0;
					avrg_results_AIS[i][j] = (double) 0;
				}
			}

			//
			// RUNNING ALGORITHMS ON MANY REPLICATES
			//

			System.out.println("GA experiment:");
			for(int p=start_pb; p<(start_pb+nb_problems); p++) {
				System.out.println("getting pb "+ p);
				// select problem
				Problem problem = JSSP.getProblem(p);
				for(int r=0; r<nb_replicates; r++) {
					//Schedule optimized_hybrid_neighborhood_solution = GA.solve_hybrid_elitist_neighborhood(problem, nb_generations, population_size, mutation_rate, 2, avrg_results_GA, r+1);
					Schedule optimized_hybrid_neighborhood_solution = GA.solve_tournament(problem, nb_generations, population_size, tournament_size, selection_pressure, mutation_rate, 2, avrg_results_GA, r+1);
				}
			}

			System.out.println("AIS experiment:");
			for(int p=start_pb; p<(start_pb+nb_problems); p++) {
				// select problem
				Problem problem = JSSP.getProblem(p);
				for(int r=0; r<nb_replicates; r++){
					Schedule optimized_ais_solution = AIS.solve(problem, nb_generations, population_size, avrg_results_AIS, r+1);
				}
			}
			System.out.println("done");

			//
			// AVERAGE && PRINT OUT RESULTS
			//

			for(int i=0; i<nb_line_res_file; i++){
				for(int j=0; j<nb_recorded_metrics; j++){
					avrg_results_GA[i][j]  	/= 	nb_replicates;
					avrg_results_AIS[i][j] 	/= 	nb_replicates;
				}
			}

			System.out.println("NB REPLICATES= "+ nb_replicates);
			System.out.println("Averaged GA results");
			for(int i=0; i<avrg_results_GA.length; i++) {
				for(int j=0; j<avrg_results_GA[0].length; j++) {
					writer.print(avrg_results_GA[i][j] + ", ");
					System.out.print(avrg_results_GA[i][j] + ", ");
				}
				writer.println();
				System.out.println();
			}

			System.out.println("NB REPLICATES= "+ nb_replicates);
			System.out.println("Averaged AIS results");
			for(int i=0; i<avrg_results_AIS.length; i++) {
				for(int j=0; j<avrg_results_AIS[0].length; j++) {
					writer.print(avrg_results_AIS[i][j] + ", ");
					System.out.print(avrg_results_AIS[i][j] + ", ");
				}
				writer.println();
				System.out.println();
			}
		} catch (Exception e) {
			System.out.println("Error loading " + "res.csv");
			e.printStackTrace();
		} finally {
			try {writer.flush();writer.close();} catch (Exception ex) {/*ignore*/}
		}
	}

	public static void main(String[] args) {
		save_performances_on_all_problems();
	}

	public static double mean(double[] p) {
	    double sum = 0;  // sum of all the elements
	    for (int i=0; i<p.length; i++) {
	        sum += p[i];
	    }
	    return sum / p.length;
	}
	
	public static double compute_variance(double[] fitness_vector, int length) {
		double mean = compute_mean(fitness_vector, length);
		double tmp = 0;
		for(int i=0; i<length; i++)
			tmp += (mean-fitness_vector[i]) * (mean-fitness_vector[i]);
		return tmp/length;
	}

	public static double compute_mean(double[] fitness_vector, int length) {
		double total = 0;
		for(int i=0; i<length; i++)
			total += fitness_vector[i];
		return total/length;
	}
}
