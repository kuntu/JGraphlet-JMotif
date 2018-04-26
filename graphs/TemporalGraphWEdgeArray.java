package graphs;
import motifs.*;

public class TemporalGraphWEdgeArray extends TemporalNetwork implements DynamicMotif {
	public int[][][] edges;
	int numNullSnapshot;
	public TemporalGraphWEdgeArray(int s, int t, int[][][] e){
		size =s;
		time =t;
		edges = e;
		numNullSnapshot  =0;
	}
	
	public GraphOfEdgeArray getSnapshot(int idx){
		if(idx>edges.length) {
			System.out.printf("\nThe %d-th snapshot does not exits, return the first snapshot instead!\n", idx);
			idx = 0;
		}
		GraphOfEdgeArray g = new GraphOfEdgeArray(edges[idx], true, size);
		return g;
	}
	
	
	@Override
	public long[][] getMotifFreq(int motifSize) {
		long[][] res = new long[time][];
		if(motifSize == 3){
			for(int t=0; t< time; t++){
				if(edges[t].length ==0){
					res[t] = new long[16];
					res[t][0] = size; 
					res[t][0] = res[t][0] * (size-1) /2 * (size -2) / 3;
					continue;
				}
				GraphOfEdgeArray g = new GraphOfEdgeArray(edges[t], true, size);
				res[t] = g.getMotifFreq(motifSize);
			}
		}else if(motifSize == 4){
			res = new long[0][199];
		}
		
		return res;
	}
	@Override
	public int[][] getMotfiTransition(int numNode) {
		// TODO Auto-generated method stub
		return null;
	}



	@Override
	public double[] getAvgMotifFreq(int motifSize) {
		// TODO Auto-generated method stub
		return null;
	}
	
	
}
