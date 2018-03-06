package randomgraph;

import java.util.Random;

import graphs.GraphOfEdgeArray;
import motifs.RandomGraphMotif;

public class RandomGraphNumMANPair implements RandomGraphModel, RandomGraphMotif {
	public int numNode;
	public int numMutual;
	public int numAsym;
	public RandomGraphNumMANPair(int nNode, int numReciprocal, int numAsymmetric) {
		// TODO Auto-generated constructor stub
		numNode = nNode;
		numMutual = numReciprocal;
		numAsym = numAsymmetric;
	}
	
	
	@Override
	public double[] getMotifFreq(int motifSize) {
		double[] res = null;
		if(motifSize == 3){
			double numEdge = numNode;
			numEdge *= numNode -1;
			double m = 2 / numEdge * numMutual;	// prob of mutual pair
			double a =  numAsym / numEdge;	//dyad prob with specific direction
			double n = 1 - m - 2 * a;	// prob of null pair
			res = new double[16];
			res[0] =  n * n * n ;
			res[1] = 6* a * n * n;
			res[2] = 3* m * n * n;
			res[3] = 3* a * a * n;
			res[4] = res[3];
			res[5] = 2 * res[3];
			res[6] = 6* m * a * n;
			res[7] = res[6];
			res[8] = 6 * a * a * a;
			res[9] = res[8] / 3;
			res[10] = 3 * m * m * n;
			res[11] = 3 * m * a * a;
			res[12] = res[11];
			res[13] = res[12] * 2 ;
			res[14] = 6 * m * m * a;
			res[15] = m * m * m;
			for(int i = 0; i< 16; i++){
				res[i] *= numEdge / 6 * (numNode -2);
			}
		}else if(motifSize == 4){
			return new double[0];
		}
		return res;
	}

	@Override
	public double[][] getMotifFreqFromSampledGraphs(int motifSize, int numOfGraphs) {
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
		Random rnd = new Random();
		int[][] edges = new int[numMutual * 2 + numAsym][];
		int[][] nodePairs = RandomGraphToolBox.sampleK_NodePairs(numNode, numMutual + numAsym);
		int idx = 0;
		for(int i = 0; i<numMutual; i++){
			edges[idx] = nodePairs[i];
			++idx;
			edges[idx] = new int[]{nodePairs[i][1], nodePairs[i][0]};
			++idx;
		}
		int tmp = 0;
		for(int i = 0; i< numAsym; i++){
			edges[idx] = nodePairs[numMutual + i];
			if(rnd.nextInt()%2 ==1){
				tmp = edges[idx][1];
				edges[idx][1] = edges[idx][0];
				edges[idx][0] = tmp;
			}
			++idx;
		}
		return new GraphOfEdgeArray(edges, true, numNode);
	}
	
	public double[][] getGraphInfo(){
		double[][] info = new double[4][1];
		info[0][0] = numNode;
		info[1][0] = numMutual;
		info[2][0] = numAsym;
		info[2][0] = info[0][0] * (info[0][0]-1)/2 - numMutual - numAsym;
		return info;
	}
}
