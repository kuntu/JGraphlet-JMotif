package randomgraph;

import graphs.GraphIO;
import mathFunctions.MathFun;

public class TemporalTimeStampRandomGraph implements TemporalRandomGraph {
	public long[][] tripletEdges;
	private long[][] backupEdges;
	public TemporalTimeStampRandomGraph(long[][] triEdges){
		tripletEdges = triEdges;
	}
	@Override
	public double[][] getProperties() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String[] getPropoertyNames() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public long[][] generateRanomEdges() {
		if(backupEdges == null){
			if(tripletEdges.length == 0 || tripletEdges[0].length == 0) return new long[0][0];
			backupEdges = new long[tripletEdges.length][tripletEdges[0].length];
		}
		TemporalRandomGraphToolBox.shuffleTripletEdgeTimeStamp(backupEdges);
		return backupEdges;
	}
	@Override
	public void selfShuffle() {
		TemporalRandomGraphToolBox.shuffleTripletEdgeTimeStamp(tripletEdges);		
	}
	@Override
	public void saveToFile(String fileName) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void saveRandomGraphs(String filePrefix, int num) {
		long[][] tmpEdges = backupEdges;
		if(tmpEdges == null) tmpEdges = tripletEdges;
		for(int i=0; i<num; ++i){
			TemporalRandomGraphToolBox.shuffleTripletEdgeTimeStamp(tmpEdges);
			GraphIO.outputMatrix(filePrefix + i +".txt", tmpEdges);
		}
	}

}
