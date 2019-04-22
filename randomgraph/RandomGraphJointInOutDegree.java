package randomgraph;

import graphs.*;
import motifs.Motif;
import motifs.MotifAlg;
import motifs.RandomGraphMotif;

public class RandomGraphJointInOutDegree implements RandomGraphModel, RandomGraphMotif{
	/**
	 * 
	 */
	public int[][] jointIODegreeSequence;
	int[][] degFreq;
	protected int numEdge;
	protected int numNode;
	public RandomGraphJointInOutDegree(int[][] degreeSeq){
		numEdge = 0;
		int tmp =0;
		for(int i=0; i< degreeSeq[0].length; i++){
			numEdge += degreeSeq[0][i];
			tmp += degreeSeq[1][i];
		}
		if(numEdge != tmp){
			System.out.println("\n In/Out Degree not matched during random graph construction. Null Graph Created Instead");
			numEdge = 0;
			jointIODegreeSequence = new int[0][2];
			numNode = 0;
		}
		jointIODegreeSequence = degreeSeq;
		numNode = degreeSeq[0].length;
		degFreq = graphs.GraphPropertiesToolBox.getJointInOutFreqFromInOutDegree(jointIODegreeSequence);
		System.out.printf("RandomGraph Created:\n\tnum nodes: %d\n\t num Edge: %d\n\t", numNode, numEdge);
		RandomGraphToolBox.getExpectedConnectProbWJointInOutDegSeq(jointIODegreeSequence[0], jointIODegreeSequence[1]);
	}
	
	public int getNumNode(){	// only allow getter function. number of nodes are decided by the degree sequence.
		return numNode;
	}
	
	public int getNumEdge(){
		return numEdge;
	}
	
	public double[][] getGraphInfo(){
		double[][] res = new double[1][]; // new double[8][];
		res[0] = RandomGraphToolBox.getExpectationOfProperties(degFreq[0], degFreq[1], degFreq[2]);
//		double[] tmp = null;
		// number of nodes
//		res[0] = new double[]{numNode};
		// number of edges
//		res[1] = new double[]{numEdge};
		// in degree sequence
//		tmp = new double[]{0.0};
//		tmp = new double[jointIODegreeSequence[0].length];
//		for(int i=0; i< tmp.length; i++) tmp[i] = jointIODegreeSequence[0][i];
//		res[2] = tmp;
		// out degree sequence
//		tmp = new double[jointIODegreeSequence[0].length];
//		for(int i=0; i< tmp.length; i++) tmp[i] = jointIODegreeSequence[1][i];
//		res[3] = tmp;
		// expected degree properties
		//tmp numNode, numEdge, expIn=expOut, expInSecMoment = 0,expOutSecMoment =0, expConnect = 0, expNunJointDeg
//		tmp = RandomGraphToolBox.getExpectationOfProperties(degFreq[0], degFreq[1], degFreq[2]);
//		res[4] = tmp;
//		res[5] = new double[degFreq[0].length];
//		res[6] = new double[degFreq[1].length];
//		res[7] = new double[res[6].length];
//		for(int i = 0; i<res[6].length; ++i){
//			res[5][i] = degFreq[0][i];
//			res[6][i] = degFreq[1][i];
//			res[7][i] = degFreq[2][i];
//		}
		return res;
	}
	
	@Override
	public GraphOfEdgeArray generateRandomGraph() {
		int[][] edges = RandomGraphToolBox.generateEdgesFromInOutDegreeSeq(jointIODegreeSequence[0], jointIODegreeSequence[1], numEdge, false);
		return new GraphOfEdgeArray(edges, true, numNode);
	}

	@Override
	public double[] getMotifFreq(int motifSize) {	// use default configuration model probability (code 0) to compute
		if(motifSize == 3){
			//return MotifAlg.getTriadFreqRandomGraphFromJoinInOutDegFreq(degFreq[0], degFreq[1], degFreq[2], numNode);
			return MotifAlg.getTriadFreqFromConfigurationModel(degFreq[0], degFreq[1], degFreq[2], numNode, 0);
		}else if(motifSize == 4){
			return new double[0];
		}else return new double[0];
	}
	public double[] getMotifFreq(int motifSize, int opt) {
		if(motifSize == 3){
			//return MotifAlg.getTriadFreqRandomGraphFromJoinInOutDegFreq(degFreq[0], degFreq[1], degFreq[2], numNode);
			return MotifAlg.getTriadFreqFromConfigurationModel(degFreq[0], degFreq[1], degFreq[2], numNode, opt);
		}else if(motifSize == 4){
			return new double[0];
		}else return new double[0];
	}

	@Override
	public double[][] getMotifFreqFromSampledGraphs(int motifSize,
			int numOfGraphs) {
		double[][] res = null;
		long[] freq = null;
		int[][] edges = null;
		if(motifSize == 3|| motifSize == -3){
			res = new double[numOfGraphs][];
			for(int t = 0; t<numOfGraphs; t++){
				edges = RandomGraphToolBox.generateEdgesFromInOutDegreeSeq(jointIODegreeSequence[0], jointIODegreeSequence[1], numEdge, false);
				edges = GraphFactory.removeLoopAndMultiEdges(edges);
				GraphOfEdgeArray gea = new GraphOfEdgeArray(edges, true, numNode);
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
