package randomgraph;

public interface TemporalRandomGraph {
	public double[][] getProperties();
	public String[] getPropoertyNames();
	public long[][] generateRanomEdges();
	public void selfShuffle();
	public void saveToFile(String fileName);
	public void saveRandomGraphs(String filePrefix, int num);
}
