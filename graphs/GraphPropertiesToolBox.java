package graphs;
import java.util.*;

import randomgraph.RandomGraphToolBox;

public class GraphPropertiesToolBox {
	/**
	 * obtain in/out degree sequences from edge array, the result is represented as a 2 by N matrix m. m[0][i] is the in degree of node i, where i\in [1, N], i=0 is consider as tmpvarialbe and has no meaning
	 * @param size num of nodes in the graph
	 * @param edges, array ofdirected edges
	 * @return a array int[2][N], inOut[0][i] is  in degree of node i+1,  inOut[1][i] is the out degree of node i+1. 
	 */
	public static int[][] getInOutDegreeFromEdges(int size, int[][] edges){
		int[][] inOut = new int[2][size];
		for(int[] edge:edges){
			inOut[0][edge[1]]++;
			inOut[1][edge[0]]++;
		}
		return inOut;
	}
	
	/**
	 * return the frequencies of joint in/out degree
	 * @param inOutDeg int[2][N+1], inOutDeg[0][i] represent the in degree of node i and inOutDeg[1][i] is the out degree of node i
	 * 
	 * @return int[3][M] array, M is the number of unique in/out pair, inOutJointFreq[0] is inDegree, [1] is out and [2] is
	 * 		the frequency of (inOutJointFreq[0], inOutJointFreq[1]) pair
	 */
	public static int[][] getJointInOutFreqFromInOutDegree(int[][] inOutDeg){
		HashMap<Long, Integer> hm = new HashMap<Long, Integer>();
		long e = 0;
		for(int i=0; i<inOutDeg[0].length; i++){
			e =(long) inOutDeg[0][i];
			e = (e<<32) + inOutDeg[1][i];
			if(!hm.containsKey(e)) hm.put(e, 1);
			else hm.put(e, hm.get(e)+1);
		}
		int[][] inOutJointFreq = new int[3][hm.size()];
		int idx= 0;
		for(Map.Entry<Long, Integer> en: hm.entrySet()){
			e = en.getKey();
			inOutJointFreq[1][idx] = (int) e;
			inOutJointFreq[0][idx] = (int) (e >> 32);
			inOutJointFreq[2][idx++] = en.getValue();
		}
		return inOutJointFreq;
	}
	
	public static double getExpectedProductInOut(int[][] jointInOutDegFreq){
		double res = 0;
		double size = 0;
		for(int i = 0; i< jointInOutDegFreq.length; i++){
			size += jointInOutDegFreq[i][0]*jointInOutDegFreq[i][2];
		}
		return res;
	}
	
	/**
	 * compute the sequence of  # reciprocal edge, # asymmetric in degree and # asymmetric out degree
	 * @param edges
	 * @return an array of 4*n size for nodes (# reciprocal edge, # asymmetric incoming edge, # asymmetric outgoing edge, node ID )
	 */
	public static int[][] getReciprocalAndInOutDegreeSequence(int[][] edges){
		int[][] trippleSeq = null;
		HashMap<Integer, HashSet<Integer>> graph = new HashMap<Integer, HashSet<Integer>>();
		HashMap<Integer, Integer> numReciprocal = new HashMap<Integer, Integer>(), inDeg = new HashMap<Integer, Integer>();
		HashSet<Integer> set = null;
		for(int[] edge: edges){
			set = graph.get(edge[0]);
			if(set == null) {
				set = new HashSet<Integer>();
				graph.put(edge[0], set);
				numReciprocal.put(edge[0], 0);
				inDeg.put(edge[0], 0);
			}
			if(!graph.containsKey(edge[1])){
				graph.put(edge[1], new HashSet<Integer>());
				numReciprocal.put(edge[1], 0);
				inDeg.put(edge[1], 0);
			}
			//reciprocal
			if(graph.get(edge[1]).contains(edge[0])){
				numReciprocal.put(edge[0], numReciprocal.get(edge[0]) + 1);
				numReciprocal.put(edge[1], numReciprocal.get(edge[1]) + 1);
				graph.get(edge[1]).remove(edge[0]);
			}else{
				//out degree
				set.add(edge[1]);
			}
			//in degree
			inDeg.put(edge[1], inDeg.get(edge[1]) + 1);
		}
		trippleSeq = new int[4][graph.size()];
		int idx = 0;
		for(int i: graph.keySet()){
			trippleSeq[3][idx] = i;		//node ID
			trippleSeq[0][idx] = numReciprocal.get(i);	// # reciprocal 
			trippleSeq[1][idx] = inDeg.get(i) - trippleSeq[0][idx];	// in degree excluding receiprocal
			trippleSeq[2][idx] = graph.get(i).size();	// out degree
			idx++;
		}
		return trippleSeq;
	}
	
	
	
	/**
	 * count the number of reciprocal, asymmetric and null node-pair
	 * @param edges edges from the graph
	 * @return array, the first element if number of nodes, second is number of reciprocal pairs and the third is the number of null pair
	 */
	public static int[] getNumReciprocalAsymmetricNullPair(int[][] edges){
		int[] res = new int[3];
		HashMap<Integer, HashSet<Integer>> simpleGraph = new HashMap<Integer, HashSet<Integer>>();
		HashSet<Integer> neighbor = null, neigh2 = null;
		for(int[] e: edges){
			if(e[0] == e[1]) continue; 
			neighbor = simpleGraph.get(e[0]);
			if(neighbor == null){
				neighbor = new HashSet<Integer>();
				simpleGraph.put(e[0], neighbor);
			}
			if(!neighbor.contains(e[1])){
				if(simpleGraph.containsKey(e[1]) && simpleGraph.get(e[1]).contains(e[0])) res[1]++;
				neighbor.add(e[1]);
			}
			neighbor = simpleGraph.get(e[1]);
			if(neighbor == null){
				neighbor = new HashSet<Integer>();
				simpleGraph.put(e[1], neighbor);
			}
		}
		res[0] = simpleGraph.size();
		for(HashSet<Integer> s: simpleGraph.values()) res[2] +=s.size();
		res[2] -= res[1] * 2;
		return res;
	}
	
	/**
	 * 
	 * @param edges
	 * @param sizeThreshold
	 * @return
	 */
	public static int[][][] obtainConnectComponentEdges(int[][] edges, int sizeThreshold){
		int[][][] res = null;
		HashMap<Integer, Integer> node2ID = new HashMap<Integer, Integer>();
		int s= 0, t = 0;
		for(int[] edge:edges){
			if(!node2ID.containsKey(edge[0])){
				s = node2ID.size() + 1;
				node2ID.put(edge[0], s);
			}else s = node2ID.get(edge[0]);
			if(!node2ID.containsKey(edge[1])){
				t = node2ID.size() + 1;
				node2ID.put(edge[1], t);
			}else t = node2ID.get(edge[1]);
			edge[0] = s;
			edge[1] = t;
		}
		//unin find to get components	// excluding 0!
		int[] par = new int[node2ID.size() + 1];
		for(int i = 1; i<par.length; i++) par[i] = i;
		int[] cnt = new int[par.length];
		for(int[] edge: edges){
			union(edge[0], edge[1], par, cnt);
		}
		HashMap<Integer, List<int[]>> components = new HashMap<Integer, List<int[]>>();
		List<int[]> ls = null;
		for(int[] edge: edges){
			s = find(edge[0], par);
			ls = components.get(s);
			if(ls == null){
				ls = new LinkedList<int[]>();
				components.put(s, ls);
			}
			ls.add(edge);
		}
		//construct edges array for temproal network
		if(sizeThreshold>0){
			List<Integer> toRemove = new LinkedList<Integer>();
			for(Map.Entry<Integer, List<int[]>> en: components.entrySet()){
				if(en.getValue().size()<sizeThreshold) toRemove.add(en.getKey());
			}
			for(int i: toRemove) components.remove(i);
		}
		
		res = new int[components.size()][][];
		s = 0;
		for(List<int[]> comp: components.values()){
			t = comp.size();
			res[s] = new int[t][];
			t = 0;
			for(int[] edge: comp){
				res[s][t] = edge;
				++t;
			}
			++s;
		}
		return res;
	}
	private static int find(int a, int[] par){
		int p = a;
		while(p!=par[p]){
			par[a] = par[p];
			a = p;
			p = par[p];
		}
		return p;
	}
	private static void union(int a, int b, int[] par, int[] cnt){
		a = find(a, par);
		b = find(b, par);
		if(a!=b){
			if(cnt[a] > cnt[b]){
				par[b] = a;
				cnt[a] += cnt[b];
			}else{
				cnt[b] += cnt[a];
				par[a] = b;
			}
		}
	}
	
	public static int[][] directedEdgeToUndirect(int[][] edges){
		int[][] res = null;
		HashSet<Long> set = new HashSet<Long>();
		for(int[] e: edges){
			set.add(RandomGraphToolBox.getUndirectedEdgeKey(e[0], e[1]));
		}
		res = new int[set.size()][];
		int idx = 0;
		int[] e = null;
		for(long l: set){
			res[idx] = new int[2];
			RandomGraphToolBox.getEdgeFromKey(l, res[idx]);
			++idx;
		}
		return res;
	}
}
