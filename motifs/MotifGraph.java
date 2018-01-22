package motifs;

import java.util.*;

import javax.sound.midi.Transmitter;

public class MotifGraph{
	public static final long EDGEBASE = 1l<<32;
	public static final String[] TRIADNAME = {"003", "012","102","021D", "021U", "021C", "111D", "111U",
		"030T", "030C", "201", "120D", "120U", "120C", "210", "300"};
	public static HashMap<Integer, Integer> triadCodeIdx = initializeTriadCode();
	public static int[][] TRIAD_EDIT_DISTANCE ={
		{0, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 6},	//1-003
		{1, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 5},	//2-012
		{2, 1, 0, 2, 2, 2, 1, 1, 3, 3, 2, 2, 2, 2, 3, 4},	//3-102
		{2, 1, 2, 0, 2, 1, 3, 1, 1, 3, 2, 2, 2, 2, 3, 4},	//4-021D
		{2, 1, 2, 2, 0, 1, 1, 3, 1, 2, 2, 2, 2, 2, 3, 4},	//5-021u
		{2, 1, 2, 1, 1, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 4},	//6-021c
		{3, 2, 1, 3, 1, 1, 0, 2, 2, 2, 1, 1, 3, 1, 2, 3},	//7-111d
		{3, 2, 1, 1, 3, 1, 2, 0, 2, 2, 1, 3, 1, 1, 2, 3}, 	//8-111u
		{3, 2, 3, 1, 1, 1, 2, 2, 0, 2, 3, 1, 1, 1, 2, 3},	//9-030t
		{3, 2, 3, 3, 2, 1, 2, 2, 2, 0, 3, 3, 3, 1, 2, 3},	//10-030c
		{4, 3, 2, 2, 2, 2, 1, 1, 3, 3, 0, 2, 2, 2, 1, 2},	//11-201
		{4, 3, 2, 2, 2, 2, 1, 3, 1, 3, 2, 0, 2, 1, 1, 2},	//12-120d
		{4, 3, 2, 2, 2, 2, 3, 1, 1, 3, 2, 2, 0, 1, 1, 2},	//13-120u
		{4, 3, 2, 2, 2, 2, 1, 1, 1, 1, 2, 1, 1, 0, 1, 2},	//14-120C
		{5, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 1}, 	//15-210
		{6, 5, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 0},	//16-300
	};
	public static final int[][] TRIAD_EDIT_DISTANCE_REPLACE ={
		{0, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 6},	//1-003
		{1, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 5},	//2-012
		{2, 1, 0, 2, 2, 2, 1, 1, 3, 3, 2, 2, 2, 2, 3, 4},	//3-102
		{2, 1, 2, 0, 2, 1, 2, 1, 1, 2, 2, 2, 2, 2, 3, 4},	//4-021D
		{2, 1, 2, 2, 0, 1, 1, 2, 1, 2, 2, 2, 2, 2, 3, 4},	//5-021u
		{2, 1, 2, 1, 1, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 4},	//6-021c
		{3, 2, 1, 2, 1, 1, 0, 1, 2, 2, 1, 1, 2, 1, 2, 3},	//7-111d
		{3, 2, 1, 1, 2, 1, 1, 0, 2, 2, 1, 2, 1, 1, 2, 3}, 	//8-111u
		{3, 2, 3, 1, 1, 1, 2, 2, 0, 1, 3, 1, 1, 1, 2, 3},	//9-030t
		{3, 2, 3, 2, 2, 1, 2, 2, 1, 0, 3, 2, 2, 1, 2, 3},	//10-030c
		{4, 3, 2, 2, 2, 2, 1, 1, 3, 3, 0, 2, 2, 2, 1, 2},	//11-201
		{4, 3, 2, 2, 2, 2, 1, 2, 1, 2, 2, 0, 2, 1, 1, 2},	//12-120d
		{4, 3, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 0, 1, 1, 2},	//13-120u
		{4, 3, 2, 2, 2, 2, 1, 1, 1, 1, 2, 1, 1, 0, 1, 2},	//14-120C
		{5, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 1}, 	//15-210
		{6, 5, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 0},	//16-300
	};
	
	public int numNode;
	public HashSet<Long> edges;
	public boolean directed;
	
	
	public MotifGraph(int num, int[][] es, boolean dir){
		numNode = num;
		edges = new HashSet<Long>();
		long tmp = 0;
		for(int i=0; i< es.length; i++){
			tmp = es[i][0];
			edges.add((long) ((tmp<<32) + es[i][1]));
		}
		directed = dir;
	}
	
	public void resetEdges(int[][] es){
		edges.clear();
		long tmp = 0;
		for(int i=0; i< es.length; i++){
			tmp = es[i][0];
			edges.add((long) ((tmp<<32) + es[i][1]));
		}
	}
	
	public static long edgeHashKey(int s, int t){
		long key = s;
		key = (key<<32) + t;
		return key;
	}
	
	public static int[] hashKeyToEdge(long k, int[] res){
		if(res==null || res.length<2)
			res = new int[2];
		res[0] = (int) (k>>32);
		res[1] = (int) (k % EDGEBASE);
		return res;
	}

	public int[][] getEdges(){
		int[][] res = new int[edges.size()][2];
		int i = 0;
		for(long e: edges){
			res[i][1] =(int) (e % EDGEBASE);
			res[i++][0] = (int) (e>>32);
		}
		return res;
	}
	
	public double[] getTriadCountFromRandomGraphWSameAsymDyads(){
		double[] freq = new double[16];
		if(numNode ==0 || edges.size()==0) return freq;
		double q = Math.log(edges.size());
		long base = numNode*(numNode-1)/2;
		double l2= Math.log(2), l3 = Math.log(3), l6 = l2 + l3, lbase = Math.log(base);
		
		q -= lbase;
		double p = Math.log(1 - Math.exp(q));
		freq[0] =  6 * p;
		freq[1] = l6 + p + q*5;
		freq[2] = l3 + 2* p + q*4;
		freq[3] = freq[2];
		freq[4] = freq[2];
		freq[5] = l6 + 2* p + q*4;;
		freq[6] = l6 + p*3 + q*3;
		freq[7] = freq[6];
		freq[8] = freq[6];
		freq[9] = l2 + p* 3 + q * 3;
		freq[10] = l3 + p*4 + q * 2;
		freq[11] = freq[10];
		freq[12] = freq[10];
		freq[13] = l6 + p * 4 + q * 2;
		freq[14] = l6 + p*5 + q;
		freq[15] = q * 6;
		lbase += Math.log(numNode-2) - l3;
		for(int i= 0; i<freq.length; i++){
			freq[i] = Math.exp(freq[i] + lbase); 
		}
		return freq;
	}
	
	public double[] getLogTriadDistrFromRandomGraphWSameDyads(){
		double[] dist = new double[16];
		double mutualNum = (double) countMutualDyads();
		double asymNum = (double) edges.size() - 2 * mutualNum;
		double nullNum = (double) numNode*(numNode-1)/2 - (mutualNum + asymNum);
		// use log
		double l2 = Math.log(2), l3 = Math.log(3);
		mutualNum = Math.log(mutualNum) - Math.log(numNode) - Math.log(numNode-1) + l2; 
		asymNum = Math.log(asymNum) - Math.log(numNode) - Math.log(numNode-1) + l2; 
		nullNum = Math.log(nullNum) - Math.log(numNode) - Math.log(numNode-1) + l2;
		dist[0] = nullNum * 3;
		dist[1] = l3 + asymNum + nullNum + nullNum;
		dist[2] = l3 + mutualNum + nullNum + nullNum;
		dist[3] = l3-4 + asymNum + asymNum + nullNum;
		dist[4] = dist[3];
		dist[5] = dist[3] + l2;
		dist[6] = l3 + mutualNum + asymNum + nullNum;
		dist[7] = dist[6];
		dist[8] = asymNum * 3 + l3 - l2 * 2;
		dist[9] = dist[8] - l3;
		dist[10] = l3 + mutualNum + mutualNum + nullNum;
		dist[11] = l3 + mutualNum+ asymNum + asymNum - l2 * 2;
		dist[12] = dist[11];
		dist[13] = dist[11] + l2;
		dist[14] = l3 + mutualNum + mutualNum + asymNum;
		dist[15] = mutualNum *3;
		return dist;
	}
	
	public double[] getTriadCountFromRandomGraphWSameDyads(){
		double[] freq = getLogTriadDistrFromRandomGraphWSameDyads();
		double lbase = Math.log(numNode) + Math.log(numNode-1) + Math.log(numNode-2)  - Math.log(6);
		for(int i=0; i<freq.length; i++) freq[i] = Math.exp(freq[i] + lbase);
		return freq;
	}
	
	
	public double[] getTriadDistrFromRandomGraphWSameDyads(){
		double[] dist = getLogTriadDistrFromRandomGraphWSameDyads();
		/*
		mutualNum /= numNode*(numNode-1)/2;
		asymNum /= numNode*(numNode-1)/2;
		nullNum /= numNode*(numNode-1)/2; // 1 -mutualNum - asymNum;
		//System.out.println(mutualNum+", "+ asymNum +", " + nullNum);
		dist[0] = Math.pow(nullNum, 3);
		dist[1] = 3 * asymNum * nullNum * nullNum;
		dist[2] = 3 * mutualNum * nullNum * nullNum;
		dist[3] = 3.0/4 * asymNum * asymNum * nullNum;
		dist[4] = dist[3];
		dist[5] = dist[3] * 2;
		dist[6] = 3 * mutualNum * asymNum * nullNum;
		dist[7] = dist[6];
		dist[8] = Math.pow(asymNum, 3) * 3 / 4;
		dist[9] = dist[8] / 3;
		dist[10] = 3 * mutualNum * mutualNum * nullNum;
		dist[11] = 3 * mutualNum* asymNum * asymNum /4;
		dist[12] = dist[11];
		dist[13] = dist[11] * 2;
		dist[14] = mutualNum * mutualNum * asymNum;
		dist[15] = Math.pow(mutualNum, 3);
		*/
		// change to probability
		for(int i = 0; i<dist.length; i++) dist[i] = Math.exp(dist[i]);
		return dist;
	}
	public int countMutualDyads(){
		int res = 0;
		for(long e: edges){
			if(edges.contains( (e<<32) + (e>>>32) )) res++;
		}
		return res/2;
	}
	
	/**
	 * obtain an edge set that contains edges appearing or disappearing at graph2 at current time step
	 * @param g1 hashed edge network in previous time step
	 * @param g2 hased edge network in current time step
	 * @param edgeSet HashSet to store the changing edges
	 * @return edgeSet, set of  changing edges
	 */
	public static HashSet<Long> getChangingEdges(MotifGraph g1, MotifGraph g2, HashSet<Long> edgeSet){
		if(edgeSet == null) edgeSet = new HashSet<Long>();
		else edgeSet.clear();
		HashSet<Long> s1=null, s2 = null;
		HashSet<Long> removed = new HashSet<Long>();
		if(g1.edges.size()>g2.edges.size()){
			s1 = g2.edges;
			s2 = g1.edges;
		}else{
			s1 = g1.edges;
			s2 = g2.edges;
		}
		
		for(long e: s2){
			if(s1.contains(e)){
				s1.remove(e);
				removed.add(e);
			}else{
				edgeSet.add(e);
			}
		}
		edgeSet.addAll(s1);
		if(s1==g2.edges){// if future graph contains few edges, need to recover the edges
			s1.addAll(removed);
		}
		return edgeSet;
	}
	
	/**
	 * triad is consider a kind of isomorphism class of 3-node motif here
	 * @param ioDegrees: int[3], each element i, denoted as ioDegrees[i] in ioDegree contains 
	 * in degree and outdegree of node i, use four bit to store the degrees, the left two bit 
	 * the indegree, the right two bit for out degree, eg, ioDegrees[i] = 10 01 (in binary), indicates
	 * that i has in-degree 2(=10) and out degree 1(=01). 
	 * ioDegree = {01 00, 01 00, 00 10} = {4, 4, 2} represent a subgraph with 2 nodes of indegree 1, outdegree 0 (with 00 01)
	 * and the other node of indegree 0, outdegree 2 (with 00 10). it is a triad 021D
	 * @return
	 */
	public static int getTriadHashKey(int[] ioDegrees){
		Arrays.sort(ioDegrees);
		int res = ioDegrees[0];
		res = (res<<4) + ioDegrees[1];
		res = (res<<4) + ioDegrees[2];
		return res;
	}
	
	
	/**
	 * An edge (s,t) is hashed to a long type e, where the left 32 bits store s and the right 32 bits store t.
	 * This function decode the hashcode e into the array sgAr at position s and t.
	 * @param e : hashcode for an edge
	 * @param sgAr : an array holds the node IDs for a subgraph
	 * @param s	: node ID for the source node of edge e
	 * @param t : node ID for the target node of edge e
	 * @return sgAr: the array holds the node IDS 
	 */
	public static int[] edgeKeyToSubgraphArray(long e, int[] sgAr, int s, int t){
		sgAr[s] = (int) (e>>32);
		sgAr[t] = (int) (e % EDGEBASE);
		return sgAr;
	}
	/**
	 * A triad has its own code as hash key, this function obtain the structure of a subgraph and detect its triad type,
	 * Then return the triad hash key.
	 * @param subGraph
	 * @param graph
	 * @param iodegrees
	 * @return
	 */
	public static int subGraphToTriadHashKey(int[] subGraph, MotifGraph graph, int[] iodegrees){
		long e = 0;
		if(iodegrees== null || iodegrees.length<3) iodegrees = new int[3];
		else Arrays.fill(iodegrees, 0);
		//obtain all edges that contained in subgraph and calculate the triad code
		for(int s=0; s<subGraph.length-1;s++){
			for(int t=s+1; t<subGraph.length; t++){
				e = edgeHashKey(subGraph[s], subGraph[t]);
				if(graph.edges.contains(e)){//s out, t in
					iodegrees[s] += 1;	// out degree ++
					iodegrees[t] += 4;	// in degree ++
				}
				e = edgeHashKey(subGraph[t], subGraph[s]);
				if(graph.edges.contains(e)){
					iodegrees[s] += 4;	// in degree ++
					iodegrees[t] += 1;	// out degree ++
				}
			}
		}
		return getTriadHashKey(iodegrees);
	}
	
	/**
	 * 
	 * @param nodes
	 * @return
	 */
	public static Set<Integer> getSubgraphkey(int[] nodes){
		//Arrays.sort(nodes);
		HashSet<Integer> hs = new HashSet<Integer>();
		for(int n: nodes) hs.add(n);
		return Collections.unmodifiableSet(hs);
	}
	
	/**
	 * 
	 * @return
	 */
	public static HashMap<Integer, Integer> initializeTriadCode(){
		HashMap<Integer, Integer> triadhm = new HashMap<Integer, Integer>();
		int[][] subgraph= {
				{0, 0, 0}, //00 00, 00 00, 00 00 -0
				{0, 1, 4}, //00 00, 00 01, 01 00 -1
				{0, 5, 5}, //00 00, 01 01, 01 01 -2
				{2, 4, 4}, //00 10, 01 00, 01 00 -3
				{1, 1, 8}, //00 01, 00 01, 10 00 -4
				{1, 4, 5}, //00 01, 01 00, 01 01 -5
				{1, 5, 9}, //00 01, 01 01, 10 01 -6
				{4, 5, 6}, //01 00, 01 01, 01 10 -7
				{2, 5, 8}, //00 10, 01 01, 10 00 -8
				{5, 5, 5}, //01 01, 01 01, 01 01 -9
				{5, 5, 10}, //01 01, 01 01, 10 10 -10
				{2, 9, 9}, //00 10, 10 01, 10 01 -11
				{6, 6, 8}, //01 10, 01 10, 10 00 -12
				{5, 6, 9}, //01 01, 01 10, 10 01 -13
				{6, 9, 10}, //01 10, 10 01, 10 10 -14
				{10, 10, 10} //10 10, 10 10, 10 10 -15
		};
		int tmp = 0;
		for(int i = 0; i<subgraph.length; i++){
			tmp = 0;
			Arrays.sort(subgraph[i]);
			for(int n: subgraph[i]) tmp = (tmp<<4) + n;
			triadhm.put(tmp, i);
			//System.out.print(tmp+"<-"+i+")\t");
		}
		if(triadCodeIdx == null) triadCodeIdx = triadhm;
		return triadhm;
	}
}