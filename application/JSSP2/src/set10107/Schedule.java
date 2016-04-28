package set10107;

import modelP.JSSP;
import modelP.Problem;

public class Schedule {

	private int [][] genotype;
	private double fitness;
	
	// ctor
	public Schedule(Problem p) {
		this.fitness = 0;
		this.genotype = JSSP.getRandomSolution(p);
	}
	
	// copy ctor
	public Schedule(Schedule to_be_replicated, Problem p){
		this.fitness = to_be_replicated.get_fitness(p);
		int genotype_height = to_be_replicated.get_genotype().length;
		int genotype_width = to_be_replicated.get_genotype()[0].length;
		int[][] genotype_replica = new int[genotype_height][genotype_width];
		for(int i=0; i<genotype_height; ++i){
			for(int j=0; j<genotype_width; ++j){
				genotype_replica[i][j] = to_be_replicated.get_genotype()[i][j];
			}
		}
		this.genotype = genotype_replica;
	}
	public void update_fitness(Problem p){
		this.fitness = JSSP.getFitness(this.genotype, p);
	}
	
	public double get_fitness(Problem p){
		update_fitness(p);
		return fitness;
	}
	
	public double get_affinity(Problem p){
		update_fitness(p);
		// use sigmoid function for output belongs to [0, 1]
		// double normalized_affinity = 1/(1+Math.exp(-fitness));
		return fitness;//normalized_affinity;
	}

	public void set_fitness(double fitness) {
		this.fitness = fitness;
	}

	public int[][] get_genotype() {
		return genotype;
	}

	public void set_genotype(int[][] genotype) {
		this.genotype = genotype;
	}
}
