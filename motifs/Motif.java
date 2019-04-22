package motifs;

public interface Motif {
	public long[] getMotifFreq(int motifSize);
	public long[][] getNodeMotifFreq(int motifSize);
}