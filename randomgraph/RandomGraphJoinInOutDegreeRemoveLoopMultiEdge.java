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
}
