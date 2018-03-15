package randomgraph;

import graphs.GraphFactory;
import graphs.GraphOfEdgeArray;
import mathFunctions.MathFun;
import motifs.RandomGraphMotif;

public class RandomGraphNumNodeEdge implements RandomGraphModel,
		RandomGraphMotif {
	private final static int[] coordinator = new int[]{1, 6, 3, 3, 3, 6, 6, 6, 6, 2, 3, 3, 3, 6, 6, 1};
	private final static boolean[] multiplyP = new boolean[]{false, true, true, false, false, false, true, false, false,false, true, false, false, false, true, true};
	public int numNode;
	public int numEdge;
	public RandomGraphNumNodeEdge(int nnode, int nedge){
		numNode = nnode;
		numEdge = nedge;
	}
	@Override
	public double[] getMotifFreq(int motifSize) {
		double[] res = null;
		if(motifSize == 3){
			res = new double[16];
			if(numNode <3) return res;
			double N = numNode;
			N = N * (N-1);
			double p = (numEdge)/N;
			double q = 1-p;
			res[0] = Math.pow(q, 6);
			for(int i = 1; i< 16; i++){
				if(multiplyP[i]){
					res[i] = res[i-1] / q * p;
				}else res[i] = res[i-1];
			}
			for(int i = 0; i< 16; i++){
				res[i] *= coordinator[i] * N / 6 * (numNode -2);
			}
		}else{
			res = new double[0];
		}
		return res;
	}

	@Override
	public double[][] getMotifFreqFromSampledGraphs(int motifSize,
			int numOfGraphs) {
		double[][] res = null;
		long[] freq = null;
		if(motifSize == 3){
			res = new double[numOfGraphs][];
			for(int t = 0; t<numOfGraphs; t++){
				GraphOfEdgeArray gea = generateRandomGraph();
				freq = gea.getMotifFreq(motifSize);
				res[t] = new double[16];
				for(int i = 0; i<freq.length; i++){
					res[t][i] = freq[i];
				}
			}
		}else if(motifSize == 4){
			return res;
		}
		return res;
	}

	@Override
	public GraphOfEdgeArray generateRandomGraph() {
		long[] mapping = MathFun.resorviorSampling(((long) numNode) * (numNode-1), numEdge);
		int[][] edges = new int[numEdge][2];
		int r = 0, c = 0, n = numNode -1;
		for(int i = 0; i< numEdge; ++i){
			r = (int) (mapping[i] / n);
			c = (int) (mapping[i] % n);
			if(c>=r) c++;
			edges[i][0] = r;
			edges[i][1] = c;
		}
		return new GraphOfEdgeArray(edges, true, numNode);
	}
	@Override
	public double[][] getGraphInfo() {
		double[][] res = new double[][]{{numNode}, {numEdge}};
		return res;
	}

}
