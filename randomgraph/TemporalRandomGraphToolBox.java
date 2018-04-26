package randomgraph;

import java.util.*;

import mathFunctions.*;
public class TemporalRandomGraphToolBox {
	
	/**
	 * return an array of triplet edges with time stamps shuffled so as to generate a new random graph.
	 * A triplet edge is an array with source node, target node, time stamp.
	 * @param tripetEdges	N*3 matrix representing triplet edges.
	 */
	public static void shuffleTripletEdgeTimeStamp(long[][] tripletEdges){
		MathFun.durstenfeldShuffleMatrixMultiColumn(tripletEdges, new int[]{0, 1}, tripletEdges.length);
	}
	
	public static void sortTripletEdges(long[][] edges){
		Arrays.sort(edges, new ArrayLongComparator(new int[]{2}, false));
	}
	
	public static long[][] removeSelfLoopMultiEdges(long[][] tripletEdges){
		int len  = tripletEdges.length;
		long curTime = Long.MIN_VALUE, edgeKey= 0;
		HashSet<Long> set = new HashSet<Long>();
		List<Integer> validIdx = new LinkedList<Integer>();
		for(int i = 0; i< len; ++i){
			if(tripletEdges[i][0] == tripletEdges[i][1]) continue;
			if(curTime<tripletEdges[i][2]){
				//update current time step
				curTime = tripletEdges[i][2];
				set.clear();
			}
			//check if there are duplicated edges
			edgeKey = (tripletEdges[i][0] << 32) + tripletEdges[i][1];
			if(!set.contains(edgeKey)){
				validIdx.add(i);
				set.add(edgeKey);
			}
		}
		if(validIdx.size() == tripletEdges.length) return tripletEdges;
		else{
			long[][] res = new long[validIdx.size()][];
			int idx = 0;
			for(int i: validIdx){
				res[idx] = tripletEdges[i];
				++idx;
			}
			return res;
		}
	}
}
