package motifs;

public interface StructProbFromDegreeEdgeCalculator {
	/**
	 * given the degrees of nodes in a subgraph, number of edges, and motif type, compute the probability of 
	 * @param deg
	 * @param numEdge
	 * @param t
	 * @return
	 */
	public double computeProbForStruct(double[][] deg, int numEdge, int t);
}
