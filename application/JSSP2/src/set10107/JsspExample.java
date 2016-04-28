package set10107;

import modelP.JSSP;
import modelP.Problem;

public class JsspExample {

	public static void main(String[] args) {

		/**
		 * Get a Problem instance using the static JSSP class in the range [1-142]
		 */
		Problem problem = JSSP.getProblem(86);		
	
		//print problem to std.out
		JSSP.printProblem(problem);
		
		/**
		 * Get a randomly initialised solution for the problem
		 */
//		int[][] solution = JSSP.getRandomSolution(problem);
		int[][] solution = {{0,1,4,3,2,5},{3,2,0,1,5,4},{0,2,4,3,5,1},{0,3,1,5,4,2},{2,0,3,4,1,5},{4,0,5,1,3,2}};
		
		
		
		//print to std.out (also prints fitness)
		JSSP.printSolution(solution, problem);
		
		/**
		 * Check the solution and return its fitness (invalid solutions return Integer.MaxValue)
		 */
		int fitness = JSSP.getFitness(solution, problem);
		System.out.println("Fitness = " + fitness);
		
		/**
		 * Save the solution (defaults to the project directory saving to a .txt file prefixed with the current computer time in milliseconds) 
		 * The filename is returned
		 */
		String filename = JSSP.saveSolution(solution, problem);
		
		/**
		 * load solution from default directory
		 */
		int[][] solution2 = JSSP.loadSolution(filename); 
	
		/**
		 * get the problem Id from the saved solution 
		 */
		int id = JSSP.getProblemIdFromSolution(filename);
		/**
		 * load the associated problem
		 */
		Problem problem2 = JSSP.getProblem(id);
		
		/**
		 * Check the solution and print to std.out 
		 */
		JSSP.printSolution(solution2, problem2);
		
		
		/**
		 * Display a saved solution graphically 
		 */
		JSSP.displaySolution(filename);	
		
		/**
		 * Create population 
		 */
		int populationSize = 100;
		int[][][] population = new int[populationSize][problem.getNumberOfMachines()][problem.getNumberOfJobs()];
		for(int i = 0; i < 100; i++){
			population[i] = JSSP.getRandomSolution(problem);
		}
		
		
		
	}
}
