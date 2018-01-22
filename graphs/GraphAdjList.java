package graphs;
import java.util.*;
import motifs.*;

public class GraphAdjList extends BasicGraph implements Motif {
	public HashMap<Integer, HashSet<Integer>> nodes;
	public GraphAdjList(int s, boolean dir, int[][] edges){
		size = s;
		directed = dir;
		nodes = new HashMap<Integer, HashSet<Integer>>();
		HashSet<Integer> tmpSet = null;
		for(int[] edge: edges){
			if(nodes.containsKey(edge[0])){
				nodes.get(edge[0]).add(edge[1]);
			}else{
				tmpSet = new HashSet<Integer>();
				tmpSet.add(edge[1]);
				nodes.put(edge[0], tmpSet);
			}
		}
		if(!directed){
			for(int[] edge: edges){
				if(nodes.containsKey(edge[1])){
					nodes.get(edge[1]).add(edge[0]);
				}else{
					tmpSet = new HashSet<Integer>();
					tmpSet.add(edge[0]);
					nodes.put(edge[1], tmpSet);
				}
			}
		}
		if(size ==-1) size = nodes.size();
	}
	@Override
	public long[] getMotifFreq(int numNode) {
		if(numNode == 3) {
			if(directed) return triadFreq();
		}
		return null;
	}

	
	private long[] triadFreq(){	// from PingHui
		long[] triadFreq = new long[16];
		for(int node: nodes.keySet()){
			for(int u: nodes.get(node)){
				for(int v: nodes.get(node)){
					if(v<=u || ((nodes.containsKey(u) && nodes.get(u).contains(v)) || (nodes.containsKey(v) && nodes.get(v).contains(u) ) && node > v ) ) continue;
					else{	// count triads
						
					}
				}
			}
		}
		return triadFreq;
	}

}
