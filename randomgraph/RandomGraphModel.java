package randomgraph;
import graphs.*;

public interface RandomGraphModel {
	public GraphOfEdgeArray generateRandomGraph();
	public double[][] getGraphInfo();
}
