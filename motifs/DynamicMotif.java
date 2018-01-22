package motifs;

public interface DynamicMotif {
	public long[][] getMotifFreq(int motifSize);
	public int[][] getMotfiTransition(int motifSize);
	public double[] getAvgMotifFreq(int motifSize);
}
