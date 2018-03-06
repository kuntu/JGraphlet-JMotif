package randomgraph;
import graphs.GraphPropertiesToolBox;

import java.io.*;
import java.util.*;

import mathFunctions.*;
import motifs.*;

public class RandomGraphToolBox {
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
		//create edges here
		int[][] res = new int[numEdge][2];
		int idx1 = 0, idx2= 0, deg = 0;
		//initialize edges, allow loop and multi-edges
		for(int i= 1; i<inDeg.length; i++){
			deg = outDeg[i];
			while(deg>0){
				res[idx1++][0] = i;
				deg--;
			}
			deg = inDeg[i];
			while(deg>0){
				res[idx2++][1] = i;
				deg--;
			}
		}
		//MathFun.durstenfeldShuffleMatrixColumn(res, 0, numEdge);
		MathFun.durstenfeldShuffleMatrixColumn(res, 1, numEdge);
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
	 * @param numNode
	 * @param numEdge
	 * @return
	 */
	public static int[][] generateUndirectedEdgesWithConfigurationModel(int[] deg, int numNode, int numEdge){
		int[][] res = new int[numEdge][2];
		int idx = 0;
		for(int i = 0; i< deg.length; i++){
			for(int j =0 ; j< deg[i]; j++){
				res[idx/2][idx%2] = i;
				++idx;
			}
		}
		MathFun.durstenfeldShuffleMatrixColumn(res, 1, res.length);
		return res;
	}
	
	public static int[][] generateDirectedEdgesWithReciprocalAndInOutDegreeTripplet(int[][] tripplet, int numNode, int numReciprocalEdge, int numAsymmetricEdge, int reSample){
		int[][] undirectEdges = generateUndirectedEdgesWithConfigurationModel(tripplet[0], numNode, numReciprocalEdge);
		if(reSample > 0) repairRandomEdgeGraph(undirectEdges, false, 0, reSample);
		int[][] res = new int[numReciprocalEdge * 2 + numAsymmetricEdge][2];
		int begIdx = undirectEdges.length * 2;
		for(int i=0;i<undirectEdges.length; i++){
			res[i*2] = undirectEdges[i];
			res[i*2+1] = new int[]{undirectEdges[i][1], undirectEdges[i][0]};
		}
		int[][] asymmetricEdges = generateEdgesFromInOutDegreeSeq(tripplet[1], tripplet[2], numAsymmetricEdge, true);
		for(int i = begIdx; i < res.length; i++){
			res[i] = asymmetricEdges[i-begIdx];
		}
		if(reSample > 0) repairRandomEdgeGraph(res, true, begIdx, reSample);
		return res;
	}
	
	
	public static void repairRandomEdgeGraph(int[][] edges, boolean dir, int begIdx, int reSampleRepeat){
		HashSet<Long> set = new HashSet<Long>();
		HashMap<Long, Integer> duplicateKeys = new HashMap<Long, Integer>();
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
		if(toRepair.isEmpty()) return;
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
				if(cnt == reSampleRepeat)
					System.out.printf("hard to repair loop and multi-edge with resample for %d times\n", cnt);
				
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
			System.out.println("num of Possible edges is larger than Max Integer. user method that supporst large scale graphs");
			return new int[0][0];
		}
		if(k > siz* (siz-1)/2) k = siz * (siz-1)/2;
		int[] sampled = MathFun.sampleKIntfromN_withNoReplacement(siz * (siz-1)/2, k);
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
				res[i][0] = tmp + 1;
			}else {
				tmp = -tmp-1;	//garantee to be >=1, 
				res[i][0] = tmp;
				--tmp;
			}
			res[i][1] = sampled[i] - al.get(tmp) + res[i][0] + 1;
		}
		return res;
	}
}
