package randomgraph;

import graphs.GraphFactory;
import graphs.GraphOfEdgeArray;

public class RandomGraphJoinInOutDegreeRemoveLoopMultiEdge extends
		RandomGraphJointInOutDegree {

	public RandomGraphJoinInOutDegreeRemoveLoopMultiEdge(int[][] degreeSeq) {
		super(degreeSeq);
	}
	@Override
	public GraphOfEdgeArray generateRandomGraph(){
		int[][] edges = RandomGraphToolBox.generateEdgesFromInOutDegreeSeq(jointIODegreeSequence[0], jointIODegreeSequence[1], numEdge, true);
		edges = GraphFactory.removeLoopAndMultiEdges(edges);
		return new GraphOfEdgeArray(edges, true, numNode);
	}
	
	@Override
	public double[][] getMotifFreqFromSampledGraphs(int motifSize,
			int numOfGraphs) {
		double[][] res = null;
		long[] freq = null;
		//int[][] edges = null;
		if(motifSize == 3 || motifSize == -3){
			res = new double[numOfGraphs][];
			for(int t = 0; t<numOfGraphs; t++){
				//edges = RandomGraphToolBox.generateEdgesFromInOutDegreeSeq(jointIODegreeSequence[0], jointIODegreeSequence[1], numEdge, true);
				GraphOfEdgeArray gea = generateRandomGraph();	//new GraphOfEdgeArray(edges, true, numNode);
				freq = gea.getMotifFreq(motifSize);
				res[t] = new double[freq.length];
				for(int i = 0; i<freq.length; i++){
					res[t][i] = freq[i];
				}
			}
		}else if(motifSize == 4){
			return res;
		}
		return res;
	}
}
