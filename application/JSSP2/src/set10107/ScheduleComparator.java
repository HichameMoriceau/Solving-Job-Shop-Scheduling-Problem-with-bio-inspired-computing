package set10107;

import java.util.Comparator;
import modelP.JSSP;
import modelP.Problem;

public class ScheduleComparator implements Comparator {

	Problem p;
	
	public ScheduleComparator(Problem problem){
		this.p = problem;
	}
	
	public int compare(Object obj1, Object obj2) {
		Schedule s1 = (Schedule) obj1;
		Schedule s2 = (Schedule) obj2;
		// compare up to 3 decimal point
		return (int) ((s1.get_fitness(p) - s2.get_fitness(p)) * 1000 );
	}

	public void set_problem(Problem prob){
		this.p = prob;
	}
}
