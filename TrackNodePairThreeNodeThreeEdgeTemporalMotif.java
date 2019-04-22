import java.util.*;

import graphs.GraphFactory;
import graphs.GraphIO;
import graphs.TemporalGraphEdgeListTimeStamp;

public class TrackNodePairThreeNodeThreeEdgeTemporalMotif {

	public static void main(String[] args) {
		String inFile = null, outFile = "outFile.txt";
		long delta = 1;
		for(String str: args) {
			if(str.startsWith("-i:")){
				inFile = str.substring(3);
			}else if(str.startsWith("-delta:")){
				delta = Long.parseLong(str.substring(7));
			}else if(str.startsWith("-o:")){
				outFile = str.substring(3);
			}
		}
		if(inFile == null){
			System.out.println("[Error]: no infile: use '-i:input_file_path' to indicate inputfile");
		}else{
			System.out.println("\tinput File = " + inFile );
			System.out.println("\toutput File = " + outFile );
			System.out.println("\ttime interval = " + delta);
		}
		TemporalGraphEdgeListTimeStamp tg = GraphFactory.getTemporalGraphWithTimeStampFromFile(inFile);
		List<long[]> trackedNodePairs = tg.trackNodePairThreeEdgeMotif(delta);
		long[][] matrix = new long[trackedNodePairs.size()][];
		int idx =0;
		for(long[] nodePair: trackedNodePairs){
			matrix[idx] = nodePair;
			++idx;
		}
		GraphIO.outputMatrix(outFile, matrix);
	}

}
