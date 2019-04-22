package randomgraph;
import graphs.GraphFactory;
import graphs.GraphPropertiesToolBox;

import java.io.*;
import java.util.*;

import mathFunctions.*;
import motifs.*;

public class RandomGraphToolBox {
	public static Comparator<int[]> arraycp = new Comparator<int[]>(){
		public int compare(int[] a, int[] b){
			if(a[0] == b[0]) return a[1] - b[1];
			else return a[0] - b[0];
		}
	};
	
	/**
	 * calculate triad distribution in a random network given the in/out degree and the number of edges
	 * @param inDeg
	 * @param outDeg
	 * @param numEdge
	 * @return
	 */
	public static double[] getTriadDistributionFromInOutDegree(double[] inDeg, double[] outDeg, int numEdge){
		double expIn = 0, expOut = 0;
		for(int i=1; i<inDeg.length; i++){
			expIn += inDeg[i] * i ;
		}
		for(int i=1; i<outDeg.length; i++){
			expOut += + outDeg[i] * i;
		}
		double edgeProb = expIn / numEdge * expOut / numEdge;
		//use model 1 to calculate triads
		
		return getLogTriadDistributionFromEdgeProb(edgeProb);
	}
	
	public static double[] getLogTriadDistributionFromEdgeProb(double edgeProb){
		double[] freq = new double[16];
		double l2= Math.log(2), l3 = Math.log(3), l6 = l2 + l3;
		double q = Math.log(edgeProb);
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
		return freq;
	}
	
	/**
	 * Assume all inputs are valid. generte 
	 * @param inDeg
	 * @param outDeg
	 * @param numEdge
	 * @return
	 */
	public static int[][] generateEdgesFromInOutDegreeSeq(int[] inDeg, int[] outDeg, int numEdge, boolean allowLoopMultiEdge){
		//assume all inputs are valid,
		//initialize edges, allow loop and multi-edges
		int[][] res = initialDirectedEdgesFromInOutDegreeSequence(inDeg, outDeg, numEdge);
		//MathFun.durstenfeldShuffleMatrixColumn(res, 0, numEdge);
		if(!allowLoopMultiEdge) repairRandomEdgeGraph(res);
		return res;
	}
	
	/**
	 * Generate random graph from in/out degree sequences. A random graph is 
	 * 		in the form of int[][] array. edges[i][0] is the source node of the i-th
	 * 		edge and edges[i][1] is the target node of the i-th edge.
	 * @param inFreq	in degree array, inFreq[i] is the frequence of in degree i. 
	 * 		i\in [0, max_degree]
	 * @param outFreq	out degree array, outFreq[i] is hte frequence of degree i.
	 * 		i\in [0, max_out_degree]
	 * @param numNode	(int) number of nodes in the random graph.
	 * @return	int[][] res: res[i][0] is the source node of the i-th
	 * 		edge and res[i][1] is the target node of the i-th edge.
	 */
	public static int[][] generateEdgesFromInOutDegreeFrequencies(int[] inFreq, int[] outFreq, int numNode, int numEdge){
		int[][] res = new int[numEdge][2];
		//assign source node
		int id = 0, cnt = 0, idx = 0;
		for(int deg = 1; deg < outFreq.length; deg++){
			cnt = outFreq[deg];
			while(cnt > 0){
				id++;
				for(int d=0; d<deg; d++){
					res[idx++][0] = id; 
				}
				cnt--;
			}
		}
		//assign target nodes
		int[] targets = getShuffledEdgeEndpoint(inFreq, numNode, numEdge);
		for(int i=0; i<res.length; i++) res[i][1] = targets[i];
		repairRandomEdgeGraph(res);
		return res;
	}
	
	/**
	 * generate undirected edges given degree sequence using configuration model, allows self loop and multi-edges between nodes
	 * @param deg
	 * @param numEdge
	 * @return
	 */
	public static int[][] generateUndirectedEdgesWithConfigurationModel(int[] deg, int numEdge){
		return initialUndirectedEdgesFromDegreeSequence(deg, numEdge);
	}
	
	/**
	 * 
	 * @param tripplet
	 * @param numNode
	 * @param numReciprocalEdge
	 * @param numAsymmetricEdge
	 * @param reSample
	 * @return
	 */
	public static int[][] generateDirectedEdgesWithReciprocalAndInOutDegreeTripplet(int[][] tripplet, int numNode, int numReciprocalEdge, int numAsymmetricEdge, int reSample){
		int[][] undirectEdges = generateUndirectedEdgesWithConfigurationModel(tripplet[0], numReciprocalEdge);
		int[][] asymmetricEdges = initialDirectedEdgesFromInOutDegreeSequence(tripplet[1], tripplet[2], numAsymmetricEdge);
		boolean success = true;
		int cnt = 0;
		if(reSample > 0) {
			success &= repairRandomEdgeGraph(undirectEdges, false, 0, reSample);
			success &= repairRandomEdgeGraph(asymmetricEdges, false, 0, reSample);
			cnt = repeairEdgesRandomGraphWReciprocalInOutSeq(undirectEdges, asymmetricEdges, reSample);
		}
		success &= (cnt>=0) ;
		int[][] res = new int[numReciprocalEdge * 2 + numAsymmetricEdge][];
		int begIdx = undirectEdges.length * 2;
		for(int i=0;i<undirectEdges.length; i++){
			res[i*2] = undirectEdges[i];
			res[i*2+1] = new int[]{undirectEdges[i][1], undirectEdges[i][0]};
		}
		for(int i = begIdx; i < res.length; i++){
			res[i] = asymmetricEdges[i-begIdx];
		}
		if(reSample>0 &&!success) {
			//System.out.println("\t reparing...");
			return removeLoopAndMultiEdges(res);
		}
		return res;
	}
	private static int[][] removeLoopAndMultiEdges(int[][] edges){
		HashSet<Long> set = new HashSet<Long>();
		long edgeCode = 0L;
		LinkedList<int[]> ls = new LinkedList<int[]>();
		for(int[] e: edges){
			if(e[0] == e[1]) continue;
			edgeCode = e[0];
			edgeCode = (edgeCode<<32) + e[1];
			if(set.add( edgeCode)) ls.add(e);
		}
		if(set.size() == edges.length) return edges;

		int[][] res =  new int[set.size()][];
		{
			int idx = 0;
			for(int[] e: ls){
				res[idx] = e;
				++idx;
			}
		}
		return res;
	}
	public static int repeairEdgesRandomGraphWReciprocalInOutSeq(int[][] unDirEdge, int[][] dirEdge, int repeat){
		int res = -1;
		int sign = 1;
		HashMap<Long, Integer> unDirMap = new HashMap<Long, Integer>(), dirMap = new HashMap<Long, Integer>();
		HashSet<Long> duplicatedKeySet = new HashSet<Long>();
		long uKey = 0, dKey = 0;
		long[] newKeys = new long[2];
		//get undirected node pairs that overlap with directed edges
		LinkedList<int[]> toRepair = new LinkedList<int[]>();
		for(int i = 0; i< unDirEdge.length; ++i){
			uKey = getUndirectedEdgeKey(unDirEdge[i][0], unDirEdge[i][1]);
			unDirMap.put(uKey, i);
		}
		for(int i = 0; i< dirEdge.length; ++i){
			dKey = getUndirectedEdgeKey(dirEdge[i][0], dirEdge[i][0]);
			if(unDirMap.containsKey(dKey)){
				toRepair.add(new int[]{unDirMap.get(dKey), i});	// recoord the idx of reciprocal-pair and asymmetric pair
				duplicatedKeySet.add(dKey);
			}
			dirMap.put(dKey, i);
		}
		// repair with rewiring
		Random rnd = new Random();
		int r = 0, cand = 0;
		for(int[] idx: toRepair){
			r = 0;
			uKey = getUndirectedEdgeKey(unDirEdge[idx[0]][0], unDirEdge[idx[0]][1]);
			if(!duplicatedKeySet.contains(uKey)) continue;
			while(r < repeat){
				if(rnd.nextInt(2) == 1 && unDirEdge.length > 1){// rewire with reciprocal(undirected) edges
					cand = rnd.nextInt(unDirEdge.length - 1);
					if(cand>= idx[0]) ++cand;
					if(canRewire(unDirEdge[idx[0]], unDirEdge[cand], unDirMap, dirMap, newKeys, false)){
						dKey = getUndirectedEdgeKey(unDirEdge[cand][0], unDirEdge[cand][1]);
						rewireEdges(unDirEdge[idx[0]], unDirEdge[cand], unDirMap, newKeys, uKey, dKey, duplicatedKeySet);
						break;
					}
				}
				//rewire directed
				if(dirEdge.length >1){
					cand = rnd.nextInt(dirEdge.length - 1);
					if(cand>= idx[1]) ++cand;
//					if(cand>=dirEdge.length || idx[0] >=dirEdge.length || idx[1] >= dirEdge.length){
//						System.out.println("something wrong");
//					}
					if(canRewire(dirEdge[idx[1]], dirEdge[cand], unDirMap, dirMap, newKeys, false)){
						dKey = getUndirectedEdgeKey(dirEdge[cand][0], dirEdge[cand][1]);
						rewireEdges(dirEdge[idx[1]], dirEdge[cand], dirMap, newKeys, uKey, dKey, duplicatedKeySet);
						break;
					}
				}
				++r;
			}
			if(r>=repeat) sign = -1;
			res += r;
		}
		return res * sign;
	}
	
	/**
	 * rewire undirected/directed edges to remove loop and multi-edges, leave the first begIdx elements in [0, begIdx-1] unchanged.
	 * @param edges
	 * @param dir	
	 * @param begIdx
	 * @param reSampleRepeat
	 */
	public static boolean repairRandomEdgeGraph(int[][] edges, boolean dir, int begIdx, int reSampleRepeat){
		HashSet<Long> set = new HashSet<Long>();
		HashMap<Long, Integer> duplicateKeys = new HashMap<Long, Integer>();
		boolean success = true;
		long key = 0, tmp=0;
		LinkedList<Integer> toRepair = new LinkedList<Integer>();
		for(int i = 0; i< edges.length; i++){
			if(edges[i][0] == edges[i][1]){
				toRepair.add(i);
			}else{
				key = edgeHashCode(edges[i][0], edges[i][1], dir);
				if(set.contains(key)){
					toRepair.add(i);
					if(duplicateKeys.containsKey(key)) duplicateKeys.put(key, duplicateKeys.get(key)+1);
					else duplicateKeys.put(key, 1);
				}else set.add(key);
			}
		}
		if(toRepair.isEmpty()) return true;
		Random rnd = new Random();
		int size = edges.length - begIdx, candidateIdx = 0, cnt = 0;
		for(int i: toRepair){
			cnt = -1;
			do{
				key = edgeHashCode(edges[i][0], edges[i][1], dir);
				if(edges[i][0] != edges[i][1] && !duplicateKeys.containsKey(key)){
					cnt = reSampleRepeat + 1;
					break;
				}
				candidateIdx = begIdx + rnd.nextInt(size);
				if(edges[candidateIdx][0] != edges[i][0] && edges[candidateIdx][1] != edges[i][1] && edges[candidateIdx][1] != edges[i][0] && edges[candidateIdx][0] != edges[i][1]){
					key = edgeHashCode(edges[candidateIdx][0], edges[i][1], dir);
					if(set.contains(key)) {
						cnt++;
						continue;
					}else{
						tmp = key;	// store previous key
						key = edgeHashCode(edges[i][0], edges[candidateIdx][1], dir);
						if(!set.contains(key)) break;
					}
				}
				cnt++;
			}while( cnt <reSampleRepeat);
			if(cnt >= reSampleRepeat){
				if(cnt == reSampleRepeat) success = false;
					//System.out.printf("hard to repair loop and multi-edge with resample for %d times\n", cnt);
				continue;
			}
			set.add(key);
			set.add(tmp);
			//reduce multi-edge. selfloop by 1.
			key = edgeHashCode(edges[i][0], edges[i][1], dir);
			if(duplicateKeys.containsKey(key)){
				if(duplicateKeys.get(key) == 1) duplicateKeys.remove(key);
				else duplicateKeys.put(key, duplicateKeys.get(key) - 1);
			}else set.remove(key);
			key = edgeHashCode(edges[candidateIdx][0], edges[candidateIdx][1], dir);
			if(duplicateKeys.containsKey(key)){
				cnt = duplicateKeys.get(key);
				if(cnt == 1) duplicateKeys.remove(key);
				else duplicateKeys.put(key, cnt-1);
			}else set.remove(key);
			//rewire two edges
			cnt = edges[i][1];
			edges[i][1] = edges[candidateIdx][1];
			edges[candidateIdx][1] = cnt;
		}
		return success;
	}
	private static long edgeHashCode(int s, int t, boolean dir){
		long res = 0;
		if(!dir && s > t){
			int tmp = s;
			s = t;
			t = tmp;
		}
		res = s;
		return (res<<32) + t;
	}
	private static LinkedList<Integer> getIdxOfInValidEdges(int[][] edges, boolean dir, int begIdx){
		HashSet<Long> set = new HashSet<Long>();
		long key = 0, key2 = 0, tmp;
		LinkedList<Integer> toRepair = new LinkedList<Integer>();
		for(int i = begIdx; i< edges.length; i++){
			if(edges[i][0] == edges[i][1]){
				toRepair.add(i);
			}else{
				key = edges[i][0];
				key2 = edges[i][1];
				if(!dir && key > key2){
					tmp = key;
					key = key2;
					key = tmp;
				}
				key = (key<<32) + key2;
				if(set.contains(key)){
					toRepair.add(i);
				}else set.add(key);
			}
		}
		return toRepair;
	}
	
	/** 
	 * reSample if there are duplicated edges or loops
	 * @param edges
	 */
	public static void repairRandomEdgeGraph(int[][] edges){
		HashMap<Integer, HashSet<Integer>> graph = new HashMap<Integer, HashSet<Integer>>();
		LinkedList<Integer> edgeIdx = new LinkedList<Integer>();
		HashMap<Long, Integer> invalid_cnt = new HashMap<Long, Integer>();
		long invalid = 0;
		//find out the idx of edges that are duplicated or are loops
		HashSet<Integer> targetSet = null;
		for(int i=0; i<edges.length; i++){
			targetSet = graph.get(edges[i][0]);
			if(targetSet == null){
				targetSet = new HashSet<Integer>();
				graph.put(edges[i][0], targetSet);
			}
			if(targetSet.contains(edges[i][1]) || edges[i][0] ==edges[i][1]) {
				addInvalid(invalid_cnt, edges[i][0], edges[i][1]);
				edgeIdx.add(i);
			}else targetSet.add(edges[i][1]);
		}
		if(edgeIdx.isEmpty()) return;
		Random rnd = new Random();
		int sIdx = 0, cnt = 0;
		for(int idx: edgeIdx){
			sIdx = rnd.nextInt(edges.length);
			cnt = 0;
			while(edges[sIdx][0] == edges[idx][1] || edges[sIdx][0] == edges[idx][0] || graph.get(edges[sIdx][0]).contains(edges[idx][1]) || edges[idx][0] == edges[sIdx][1] || edges[idx][1] == edges[sIdx][1] || graph.get(edges[idx][0]).contains(edges[sIdx][1])){
				sIdx = rnd.nextInt(edges.length);
				if(++cnt ==  1000){
					System.out.println("difficult to repair edge: num of fail"+ cnt);
					return ;
				}
			}
			cnt = edges[sIdx][1];	//cnt store temporal target for chosen node pair
			edges[sIdx][1] = edges[idx][1];
			edges[idx][1] = cnt;
			targetSet = graph.get(edges[idx][0]);
			targetSet.add(edges[idx][1]);
			if(removeInvalid(invalid_cnt, edges[idx][0], edges[sIdx][1])){
				targetSet.remove(edges[sIdx][1]);
			}
			targetSet = graph.get(edges[sIdx][0]);
			targetSet.add(edges[sIdx][1]);
			if(removeInvalid(invalid_cnt, edges[sIdx][0], edges[idx][1])){
				targetSet.remove(edges[idx][1]);
			}
		}
	}
	private static void addInvalid(HashMap<Long, Integer> hm, int s, int t){
		long key = s;
		key = (key <<32) + t;
		if(hm.containsKey(key)) hm.put(key, hm.get(key) + 1);
		else hm.put(key, 1);
	}
	private static boolean removeInvalid(HashMap<Long, Integer> hm, int s, int t){
		long key = s;
		key = (key << 32) + t;
		if(hm.containsKey(key)){
			int cnt = hm.get(key);
			if(cnt == 1){
				hm.remove(key);
				return true;
			}else{
				hm.put(key, cnt-1);
				return false;
			}
		}
		return true;
	}
	
	/**
	 * 
	 * @param a
	 * @param repeat
	 * @return
	 */
	public static int repairEdgesAsOneArray(int[] a, int repeat){
		int cnt = 0;
		int len = a.length;
		if(len%2!=0){
			System.out.println("edge array not valid");
			return -1;
		}
		len /= 2 ;
		HashSet<Long> keySet = new HashSet<Long>();
		HashMap<Long, Integer> duplicatedCount = new HashMap<Long, Integer>();
		LinkedList<Integer> toResample = new LinkedList<Integer>();
		long key = 0, candKey;
		for(int i =0; i< len; ++i){
			if(a[i*2] == a[i*2+1]){	// loop edge
				toResample.add(i);
				continue;
			}
			//multi-edge
			key = getUndirectedEdgeKey(a[i*2], a[i*2+1]);
			if(keySet.contains(key)) {
				toResample.add(i);
				if(duplicatedCount.containsKey(key)) duplicatedCount.put(key, duplicatedCount.get(key) + 1);
				else duplicatedCount.put(key, 1);
			}else keySet.add(key);
		}
		// resample for duplicated edges
		Random rnd = new Random();
		int cand = -1, r = 0, sIdx= 0, tIdx = 0;
		--len;	//use for resample len-1 edge to rewire edges
		int[] curEdge = new int[2], candEdge = new int[2];
		long[] newKeys = new long[2];
		for(int i: toResample){
			r = 0;
			sIdx = 2*i;
			tIdx = 2 * i + 1;
			curEdge[0] = a[sIdx];
			curEdge[1] = a[tIdx];
			key = getUndirectedEdgeKey(curEdge[0], curEdge[1]);
			if(!duplicatedCount.containsKey(key)) continue;	// an invalid edge could rewired with others in previous step. no need to rewire  
			while(r < repeat){
				cand = rnd.nextInt(len);
				if(cand>= i) ++cand;
				candEdge[0] = a[sIdx];
				candEdge[1] = a[tIdx];
				if(!canRewire(curEdge, candEdge, keySet, newKeys, false)){
					++r;
				}else{
					candKey = getUndirectedEdgeKey(candEdge[0], candEdge[1]);
					rewireEdges(curEdge, candEdge, keySet, newKeys, key, candKey, duplicatedCount);
					a[tIdx] = curEdge[1];
					a[cand * 2 + 1] = candEdge[1];
					break;
				}
			}
			cnt += r;
			if(r>= repeat) System.out.println("[Warning]: Remove a loop or multidege-edge fail with resampling less than " + repeat + " times\n\t increase resampling time" );
		}
		return cnt;
	}
	public static long getEdgeKey(int s, int t){
		long res = s;
		return (res<<32) + t;
	}
	public static long getUndirectedEdgeKey(int s, int t){
		if(s> t) return getEdgeKey(t, s);
		else return getEdgeKey(s,t );
	}
	public static void getEdgeFromKey(long key, int[] edge){
		edge[1] = (int) (key%(1L<<32));
		edge[0] = (int) (key>>32);
	}
	/**
	 * check if two edges can be re-wired given the current (un)directed edges
	 * @param curEdge
	 * @param candEdge
	 * @param edgeKeys
	 * @param keys
	 * @param directedEdge
	 * @return
	 */
	private static boolean canRewire(int[] curEdge, int[] candEdge, HashSet<Long> edgeKeys, long[] keys, boolean directedEdge){
		if(curEdge[0] == candEdge[0] || curEdge[0] == candEdge[1] || curEdge[1] == candEdge[0] || curEdge[1] == candEdge[1]) return false;
		keys[0] = directedEdge? getEdgeKey(curEdge[0], candEdge[1]) : getUndirectedEdgeKey(curEdge[0], candEdge[1]);
		if(edgeKeys.contains(keys[0])) return false;
		keys[1] = directedEdge? getEdgeKey(candEdge[0], curEdge[1]) : getUndirectedEdgeKey(candEdge[0], curEdge[1]);
		if(edgeKeys.contains(keys[1])) return false;
		return true;
	}
	private static boolean canRewire(int[] curEdge, int[] candEdge, HashMap<Long, Integer> map1, HashMap<Long, Integer> map2, long[] keys, boolean directedEdge){
		if(curEdge[0] == candEdge[0] || curEdge[0] == candEdge[1] || curEdge[1] == candEdge[0] || curEdge[1] == candEdge[1]) return false;
		keys[0] = directedEdge? getEdgeKey(curEdge[0], candEdge[1]) : getUndirectedEdgeKey(curEdge[0], candEdge[1]);
		if(map1.containsKey(keys[0]) || map2.containsKey(keys[0])) return false;
		keys[1] = directedEdge? getEdgeKey(candEdge[0], curEdge[1]) : getUndirectedEdgeKey(candEdge[0], curEdge[1]);
		if(map2.containsKey(keys[1]) || map2.containsKey(keys[1])) return false;
		return true;
	}
	private static boolean rewireEdges(int[] curEdge, int[] candEdge, HashSet<Long> edgeKeys, long[] newKeys, long curKey, long candKey, HashMap<Long, Integer> duplicatedCount){
		//remove key of candidate edge
		int cnt = 0;
		if(duplicatedCount.containsKey(candKey)){
			cnt = duplicatedCount.get(candKey);
			if(cnt == 1) duplicatedCount.remove(candKey);
			else duplicatedCount.put(candKey, cnt - 1);
		}else edgeKeys.remove(candKey);
		//remove curKey
		if(duplicatedCount.containsKey(curKey)){
			cnt = duplicatedCount.get(curKey);
			if(cnt == 1) duplicatedCount.remove(curKey);
			else duplicatedCount.put(curKey, cnt - 1);
		}
		//add new keys
		for(long l: newKeys) edgeKeys.add(l);
		cnt = curEdge[1];
		curEdge[1] = candEdge[1];
		candEdge[1] = cnt;
		return true;
	}
	private static boolean rewireEdges(int[] curEdge, int[] candEdge, HashMap<Long, Integer> keyMap, long[] newKeys, long curKey, long candKey, HashSet<Long> duplicatedKeySet){
		// current Edge 
		int tmp  = keyMap.get(curKey);
		tmp = keyMap.remove(curKey);
		keyMap.put(newKeys[0], tmp);
		//candidate edge
		tmp = keyMap.get(candKey);
		keyMap.remove(candKey);
		keyMap.put(newKeys[1], tmp);
		duplicatedKeySet.remove(curKey);
		//rewire target in edges;
		tmp = curEdge[1];
		curEdge[1] = candEdge[1];
		candEdge[1] = tmp;
		return true;
	}
	
	/**
	 * get a degree sequence for a network given the degree frequences 
	 * @param degFreq degree frequences, degFreq[i] is the number of nodes with degree i
	 * @param isRandom	randomly assign the degree to nodes
	 * @param num number of nodes in the graph
	 * @return a degree sequence s: s[u] is the degree of node u. s[0] is null node
	 */
	public static int[] getDegSeqFromDegFreq(int[] degFreq, boolean isRandom, int num){
		//sum(degFreq) must equal to num
		int[] res = new int[num+1];
		int nID = 0;
		int[] idAr = new int[num];	// idx -> node mapping
		for(int i = 0; i<num; i++) idAr[i] = i + 1;
		if(isRandom){
			MathFun.durstenfeldShuffle(idAr, idAr.length);
			nID = degFreq[0];	// initialize the first node ID that has non-zero indegree
		}
		int numNodes = 0;
		for(int i=1; i<degFreq.length; i++){
			numNodes = degFreq[i];
			while(numNodes>0){
				res[idAr[nID++]] = i;
				numNodes--;
			}
		}
		return res;
	}
	
	/**
	 * given the size of a network and the edges, return the in/out degree frequences
	 * @param edges
	 * @param size
	 * @return
	 */
	public static int[][] getInOutDegreeFreqFromEdges(int[][] edges, int size){
		int[][] res = new int[2][];
		int[][] inoutDegs = new int[size+1][2];
		int maxIn=0, maxOut = 0;
		for(int[] edge: edges){
			inoutDegs[edge[0]][0]++;
			if(maxOut < inoutDegs[edge[0]][0]) maxOut = inoutDegs[edge[0]][0];
			inoutDegs[edge[1]][1]++;
			if(maxIn < inoutDegs[edge[1]][1]) maxIn = inoutDegs[edge[1]][1];
		}
		res[0] = new int[maxIn+1];
		res[1] = new int[maxOut+1];
		for(int i = 1; i< inoutDegs.length; i++){
			res[1][inoutDegs[i][0]] ++;
			res[0][inoutDegs[i][1]] ++;
		}
		return res;
	}
	/**
	 * return an array of randomly ordered edge endpoints of one side for a graph given their degree distribution
	 * @param degFreq 
	 * @param numNode
	 * @param numEdge
	 * @return
	 */
	private static int[] getShuffledEdgeEndpoint(int[] degFreq, int numNode, int numEdge){
		int[] res = new int[numEdge];
		int num2Shuffle = numNode - degFreq[0];
		int[] ids = new int[numNode];
		for(int i = 0; i<numNode; i++) ids[i] = i+1;
		MathFun.durstenfeldShuffle(ids, num2Shuffle);
		int idxEdge = 0, idxID = 0, cnt = 0;
		for(int d = 1; d<degFreq.length; d++){
			cnt = degFreq[d];	// num of nodes with degree d
			while(cnt>0){
				for(int i= d; i>0; i--){
					res[idxEdge++] = ids[idxID];
				}
				cnt--;
				idxID++;
			}
		}
		MathFun.durstenfeldShuffle(res, res.length);
		return res;
	}
	
	/** 
	 * generate k random graph given the in/out degree of a directed graph and save them as temporal networks into a file
	 * @param edges
	 * @param size
	 * @param k
	 * @param fileName
	 */
	public static void outputKRandomGraphWsameInOutDegreeFreq(int[][] edges, int size, int k, String fileName){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
			int[][] inoutDegFreq = getInOutDegreeFreqFromEdges(edges, size);
			int[][] graph = null;
			StringBuilder sb = new StringBuilder();
			bw.write(size + " " + k+"\n");
			for(int i=0; i<k; i++){
				graph = generateEdgesFromInOutDegreeFrequencies(inoutDegFreq[0], inoutDegFreq[1], size, edges.length);
				sb.setLength(0);
				for(int[] edge: graph){
					if(sb.length()>0) sb.append(' ');
					sb.append(edge[0]);
					sb.append(' ');
					sb.append(edge[1]);
				}
				sb.append("\n");
				bw.write(sb.toString());
			}
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static double getExpectedProdInOutDegWJointInOutDegFreq(int[][] jointInOutDegFreq, int size){
		double res  = 0, tmp = 0;
		double dv = 1.0 / size;
		for(int i = 0 ; i<jointInOutDegFreq.length; i++){
			tmp = 0;			
			for(int o = 0; o <jointInOutDegFreq.length; o++){
				if(o == i) continue;
				tmp += dv*jointInOutDegFreq[o][2] * jointInOutDegFreq[o][1];
			}
			res += dv * jointInOutDegFreq[i][0] * jointInOutDegFreq[i][2] * tmp;
		}
		return res;
	}
	
	public static double getExpectedConnectProbWJointInOutDegSeq(int[] in, int[] out){
		double res =0, tmp = 0;
		for(int i=1; i<in.length; i++){
			res += in[i];	//number of edges
			tmp += in[i]*out[i];
		}
		res = (res*res- tmp)/(in.length-1)/(in.length-2)/res;
		System.out.printf("Expected connect probability: %f\n", res);
		return res;
	}
	public static double[] getExpectationOfProperties(int[] in, int[] out, int[] cnt){
		HashSet<Integer> inHs = new HashSet<Integer>(), outHs = new HashSet<Integer>();
		double[] res = new double[12]; //see the printed output for meaning
		for(int i=0; i<in.length; i++) {
			res[0] += cnt[i];	//number of nodes
			res[1] += in[i] * cnt[i];	// number of edges
			res[7] = Math.max(res[7], in[i]);
			res[8] = Math.max(res[8], out[i]);
			inHs.add(in[i]);
			outHs.add(out[i]);
		}
		res[2] = res[1] / res[0];	//expected in/out degree
		for(int i=0; i<in.length; i++){
			res[3] += Math.min(1, in[i] * in[i]/res[0]) * cnt[i];	//second moment of in degree
			res[4] += Math.min(1, out[i] * out[i] / res[0]) * cnt[i];	// second moment of out degree
			res[11] += Math.min(1, out[i] * in[i] / res[0]) * cnt[i];	// to compute expected self-loop.
		}
		res[11] /= res[2];
		double tmp = res[0] * (res[0] -1);	//	max Num of edges;
		for(int i=0; i< in.length; i++){
			for(int j = 0; j<out.length; j++){
				if(in[i]!=0 && out[j]!=0){
					res[5] += in[i] * out[j]/res[1] *cnt[i] * cnt[j]/tmp;
				}
			}
		}
		res[6] = in.length;
		res[9] = inHs.size();
		res[10] = outHs.size();
		//print out info
		System.out.printf("\n Expectations:"
				+ "\n\t num Node: %f"
				+ "\n\t num Edge: %f"
				+ "\n\t expected In/Out Degree: %f" 
				+ "\n\t Second Moment In Degree: %f"
				+ "\n\t Second Moment out Degree: %f"
				+ "\n\t expected connection probability: %f"
				+ "\n\t num unique joint in/out degree: %f"
				+ "\n\t max In Degree: %f" 
				+ "\n\t max out degree: %f"
				+ "\n\t num of unique in degree: %f"
				+ "\n\t num of unique out degree: %f"
				+ "\n\t expected num of self loop :%f\n", res[0], res[1], res[2], res[3], res[4], res[5], res[6], res[7], res[8], res[9], res[10], res[11]);
		return res;
	}
	public static int[][] sampleK_NodePairs(int siz, int k){
		long N = siz;
		N = N*(N-1)/2;
		if(N > ((long) Integer.MAX_VALUE)){
			System.out.println("[error] num of Possible edges is larger than Max Integer. use method that supporst large scale graphs");
			return new int[0][0];
		}
		if(k > N) k = (int) N;
		int[] sampled = MathFun.sampleKIntfromN_withNoReplacement((int) N, k);
		int[][] res = new int[k][2];
		ArrayList<Integer> al = new ArrayList<Integer>();
		int tmp =0;
		al.add(tmp);
		while(siz>1){
			--siz;
			tmp += siz;
			al.add(tmp);
		}
		for(int i =0; i< k; i++){
			tmp = Collections.binarySearch(al, sampled[i]);
			if(tmp > -1){
				res[i][0] = tmp;
			}else {
				tmp = -tmp-1;	//garantee to be >=1, 
				res[i][0] = tmp-1;
				--tmp;
			}
			res[i][1] = sampled[i] - al.get(tmp) + res[i][0]+1;
		}
		return res;
	}
	
	//--------- supporting function
	/**
	 * 
	 * @param inSeq
	 * @param outSeq
	 * @return
	 */
	public static int[][] initialDirectedEdgesFromInOutDegreeSequence(int[] inSeq, int[] outSeq, int numEdge){
		int[][] res = new int[numEdge][2];
		int idx1 = 0, idx2= 0, deg = 0;
		//initialize edges, allow loop and multi-edges
		for(int i= 0; i<inSeq.length; i++){
			deg = outSeq[i];
			while(deg>0){
				res[idx1++][0] = i;
				deg--;
			}
			deg = inSeq[i];
			while(deg>0){
				res[idx2++][1] = i;
				deg--;
			}
		}
		MathFun.durstenfeldShuffleMatrixColumn(res, 1, numEdge);
		return res;
	}
	
	public static int[][] initialUndirectedEdgesFromDegreeSequence(int[] seq, int numEdge){
		int[] a = new int[numEdge * 2];
		int deg = 0, idx = 0;
		for(int i = 0; i< seq.length; ++i){
			deg = seq[i];
			while(deg > 0){
				a[idx] = i;
				++idx;
				--deg;
			}
		}
		MathFun.durstenfeldShuffle(a, a.length);
		int[][] res = new int[numEdge][2];
		for(int i = 0; i<numEdge; ++i){
			res[i][0] = a[i*2];
			res[i][1] = a[i*2 + 1];
			if(res[i][0] > res[i][1]){
				deg = res[i][0];
				res[i][0] = res[i][1];
				res[i][1] = deg;
			}
		}
		repairRandomEdgeGraph(res, true, 0, 1000);
		return res;
	}
	
	/**
	 * apply Havel-Hakimi algorithm to generate simple undirected graph:
	 * steps:
	 * 1. sort degree in non-decreasing order
	 * 2. remove the node with the largest degree d_i and connect it to the next d_i nodes (with idx=i-1,...,i-d_i)
	 * 3. each of the d_i node reduce its degree by one
	 * @param degSeq
	 * @return
	 */
	public static int[][] havelHakimiUndirect(int[] degSeq){
		if(degSeq.length == 1) return new int[0][0];
		int[][] edges = new int[degSeq.length][2];
		int numEdge = 0;
		for(int i = 0; i < degSeq.length; ++i){
			edges[i][0] = degSeq[i];
			edges[i][1] = i;
			numEdge += degSeq[i];
		}
		Arrays.sort(edges, arraycp);
		int[][] res = new int[numEdge][];
		//connect the node with largest degree to the next few nodes
		int swapIdx = 0, idx= 0, end = 0, beg = 0;
		int[] tmp = null;
		for(int i = 0; i<edges.length && edges[i][0]>0; ++i){
			//connect nodes;
			end = i + edges[i][0];
			beg = i+1;
			swapIdx = beg;
			while(swapIdx < edges.length && edges[beg] == edges[swapIdx]){
				++swapIdx;
			}
			--swapIdx;
			for(int j = beg; j<=end; ++j){
				--edges[j][0];
				res[idx] = new int[]{edges[i][1], edges[j][1]};
				++idx;
			}
			edges[i][0] = 0;
			//swap
			for(int j = beg; j<=end; ++j){
				if(j<swapIdx && edges[j][0] >= edges[swapIdx][0]) continue;
				swapIdx = Math.max(swapIdx, j);
				if(swapIdx == j){
					while(swapIdx<edges.length &&  edges[j][0] +1 == edges[swapIdx][0] + (swapIdx<=end?1: 0)){
						++swapIdx;
					}
					--swapIdx;
				}
				if(swapIdx <= j) continue;
				tmp = edges[swapIdx];
				edges[swapIdx] = edges[j];
				edges[j] = tmp;
				--swapIdx;
			}
		}
		return res;
	}
}
