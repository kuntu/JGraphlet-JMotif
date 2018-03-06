package motifs;

public interface RandomGraphMotif {
	public double[] getMotifFreq(int motifSize);
	public double[][] getMotifFreqFromSampledGraphs(int motifSize, int numOfGraphs);
	
}
