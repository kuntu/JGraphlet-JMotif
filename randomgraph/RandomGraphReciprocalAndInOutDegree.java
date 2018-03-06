package randomgraph;

import graphs.GraphFactory;
import graphs.GraphOfEdgeArray;
import motifs.RandomGraphMotif;

public class RandomGraphReciprocalAndInOutDegree implements RandomGraphModel,
		RandomGraphMotif {
	private final static int RESAMPLE_REPEAT = 1000;
	public int[][] trippletSeq;
	public int numNode;
	public int numReciprocal;
	public int numAsymmetricEdge;
	
	/**
	 * 
	 * @param ts tripplet sequence of reciprocal edge, incoming edge and outgoing edge
	 */
	public RandomGraphReciprocalAndInOutDegree(int[][] ts){
		numNode = ts[0].length;
		trippletSeq = new int[3][numNode];
		trippletSeq[0] = ts[0];
		trippletSeq[1] = ts[1];
		trippletSeq[2] = ts[2];
		for(int i = 0; i< trippletSeq[0].length; i++){
			numReciprocal += trippletSeq[0][i];
			numAsymmetricEdge += trippletSeq[1][i];
		}
		if((numReciprocal % 2 ) == 1){
			numNode = 0;
			numReciprocal = 0;
			numAsymmetricEdge = 0;
			trippletSeq = new int[4][0];
		}else numReciprocal /= 2;
	}
	@Override
	public double[] getMotifFreq(int motifSize) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[][] getMotifFreqFromSampledGraphs(int motifSize, int numOfGraphs) {
		double[][] res = null;
		long[] freq = null;
		int[][] edges = null;
		if(motifSize == 3){
			res = new double[numOfGraphs][];
			for(int t = 0; t<numOfGraphs; t++){
				edges = RandomGraphToolBox.generateDirectedEdgesWithReciprocalAndInOutDegreeTripplet(trippletSeq, motifSize, numReciprocal, numAsymmetricEdge, RESAMPLE_REPEAT);
				edges = GraphFactory.removeLoopAndMultiEdges(edges);
				GraphOfEdgeArray gea = new GraphOfEdgeArray(edges, true, numNode);
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
		//
		int[][] edges = RandomGraphToolBox.generateDirectedEdgesWithReciprocalAndInOutDegreeTripplet(trippletSeq, numNode, numReciprocal, numAsymmetricEdge, RESAMPLE_REPEAT);	
		return new GraphOfEdgeArray(edges, true, numNode);
	}
	@Override
	public double[][] getGraphInfo() {
		double[][] res = new double[6][];
		res[0] = new double[]{numNode};
		res[1] = new double[]{numReciprocal};
		res[2] = new double[]{numAsymmetricEdge};
		for(int i = 3; i < 6; i++){
			res[i] = new double[trippletSeq[i-3].length];
			for(int j =0; j< res[i].length; j++) res[i][j] = trippletSeq[i-3][j];
		}
		return res;
	}

}
