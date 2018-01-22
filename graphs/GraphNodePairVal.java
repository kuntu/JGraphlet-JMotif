package graphs;
import java.util.*;
import java.util.Map.Entry;
import motifs.*;

public class GraphNodePairVal extends BasicGraph implements Motif {
	//nodes.key is the first node a, node.hm.key is the second node b of an edge, 
	//	nodes.hm.val is the dyad type of the edge, 1 represent a->b, 2 preprsent a<-b
	//	and 3 represent a<-->b.
	public HashMap<Integer, HashMap<Integer, Integer>> nodes;
	private int[] tmpArray;
	private static final int[][] triadIODegIdx = {{0, 1}, {0, 2}, {1, 2}};
	private static int[] triadIODeg = new int[3];
	private static int[] triadStructure = new int[3];
	public GraphNodePairVal(int s, boolean dir, int[][] edges){
		size = s;
		directed = dir;
		nodes = new HashMap<Integer, HashMap<Integer, Integer>>();
		HashMap<Integer, Integer> neig = null;
		for(int[] edge: edges){
			if(nodes.containsKey(edge[0])) neig = nodes.get(edge[0]);
			else{
				neig = new HashMap<Integer, Integer>();
				nodes.put(edge[0], neig);
			}
			//assume there is no duplicated edges, if edge[0] already has edge[1] as neighbor, then there must be an edge edge[1]->edge[2] 
			if(neig.containsKey(edge[1])) neig.put(edge[1], 3);
			else neig.put(edge[1], 1);	// 1 represent edge[0] -> edge[1]
			// 
			if(nodes.containsKey(edge[1])) neig = nodes.get(edge[1]);
			else{	//
				neig = new HashMap<Integer, Integer>();
				nodes.put(edge[1], neig);
			}
			if(neig.containsKey(edge[0])) neig.put(edge[0], 3);
			else neig.put(edge[0], 2);	// 2 represent edge[1] has a neighbor edge[1] <-- edge[0]
		}
		//if(size < nodes.size()) size = nodes.size(); 
	}
	@Override
	public long[] getMotifFreq(int motifSize) {
		if(tmpArray==null || tmpArray.length!= motifSize) tmpArray = new int[motifSize];
		if(directed){
			if(motifSize ==3) {
//				MotifGraph.initializeTriadCode();
				// tmpArray use as in/out degree array for tiple nodes
				return getTriadFreq();
			}else if(motifSize ==4){
				return getFourNodeMotifFreq(tmpArray);
			}
		}
		return null;
	}


	
	public int getPairVal(int u, int v){
		if(!nodes.containsKey(u) || !nodes.get(u).containsKey(v)) return 0;
		else return nodes.get(u).get(v);
	}
	
	/**
	 * count triad frequency
	 * @return
	 */
	private long[] getTriadFreq(){
		// from Pinghui
		long[] res = new long[16];
		HashMap<Integer, Integer> neig = null;
		int edgeVal = -2, motifType = -1;
		int node = -1, u=-1, v = -1, nu = -1, reCnt=0;
		boolean isTriangle = false;
		for(Entry<Integer, HashMap<Integer, Integer>> en: nodes.entrySet()){
			node = en.getKey();
			neig = en.getValue();
			for(Entry<Integer, Integer> eu: neig.entrySet()){
				u = eu.getKey();
				nu = eu.getValue();
				reCnt = 0;
				for(Entry<Integer, Integer> ev: neig.entrySet()){
					v = ev.getKey();
					if(u==v) continue;
					isTriangle = nodes.get(u).containsKey(v);
					if(isTriangle) reCnt++;
					if(u<v || (isTriangle && node > v ) ) continue;	// only count LEADER node u and skip node 
					else{
						//decide motifType according to triad structure
						Arrays.fill(triadIODeg, 0);
						if((nu & 1) >0){
							triadIODeg[0] +=1;
							triadIODeg[1] +=4;
						}
						if((nu & 2)>0) {
							triadIODeg[0] +=4;
							triadIODeg[1] +=1;
						}
						edgeVal = ev.getValue();
						if((edgeVal& 1) >0){
							triadIODeg[0] +=1;
							triadIODeg[2] +=4;
						}	
						if((edgeVal & 2)>0) {
							triadIODeg[0] +=4;
							triadIODeg[2] +=1;
						}
						if(isTriangle){
							edgeVal = nodes.get(v).get(u);
							if((edgeVal& 1) >0){
								triadIODeg[2] +=1;
								triadIODeg[1] +=4;
							} 
							if((edgeVal & 2)>0) {
								triadIODeg[2] +=4;
								triadIODeg[1] +=1;
							}
						}
						motifType = MotifGraph.triadCodeIdx.get(MotifGraph.getTriadHashKey(triadIODeg));
						res[motifType]++;
					}
				}
				if(u>node){
					if(neig.get(u)==3) motifType = 2;
					else motifType =1;
					res[motifType] += size - neig.size() - nodes.get(u).size() + reCnt;
				}
			}
		}
		long total = (long) size;
		total = total *(total-1)/2 * (total -2) / 3;
		for(long i: res) total -= i;
		res[0] = total;
		return res;
	}
	
	public int[] updateMotifFreq(GraphNodePairVal nextGraph, int[] curFreq){
		int[] nextFreq = new int[curFreq.length];
		HashMap<Integer, Integer> neig = null;
		HashMap<Integer, Integer> preNeig = null;
		for(int node: nodes.keySet()){
			neig = nodes.get(node);
			preNeig = nextGraph.nodes.get(node);
			for(int u: neig.keySet()){
				//TODO if(preNeig.containsKey(u) && preNeig.get(u))
			}
		}
		return nextFreq;
	}
	
	/**
	 * count 4-node motifs
	 * @param motifTypeArray
	 * @return
	 */
	private long[] getFourNodeMotifFreq(int[] motifTypeArray){
		long[] res = new long[199];
		
		return res;
	}
	
	/**
	 * set the triad structure for node, u and v. where u and v are in the neigborhood of node;
	 * @param node
	 * @param u
	 * @param v
	 */
	public void setTriadStructure(int node, int u, int v){
		Arrays.fill(triadStructure, 0);
		triadStructure[0] = nodes.get(node).get(u);
		triadStructure[1] = nodes.get(nodes).get(v);
		if(nodes.get(u).containsKey(v)) triadStructure[2] = nodes.get(u).get(v);
	}
	
	public static void setTriadIODegree(){
		Arrays.fill(triadIODeg, 0);
		for(int i=0; i<triadStructure.length; i++){
			if(triadStructure[i] == 1){
				triadIODeg[triadIODegIdx[i][0]] += 1;
				triadIODeg[triadIODegIdx[i][1]] += 4;
			}else if(triadStructure[i] == 2){
				triadIODeg[triadIODegIdx[i][0]] += 4;
				triadIODeg[triadIODegIdx[i][1]] += 1;
			}else if(triadStructure[i] == 3){
				triadIODeg[triadIODegIdx[i][0]] += 5;
				triadIODeg[triadIODegIdx[i][1]] += 5;
			}
		}
	}
	
	public static int getTriadTypeFromStructureVal(int uv, int uw, int vw){
		triadStructure[0] = uv;
		triadStructure[1] = uw;
		triadStructure[2] = vw;
		setTriadIODegree();
		return MotifGraph.triadCodeIdx.get(MotifGraph.getTriadHashKey(triadIODeg));
	}
}
