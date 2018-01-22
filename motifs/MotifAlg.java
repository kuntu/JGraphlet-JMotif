package motifs;
import java.util.Arrays;


public class MotifAlg {
	public static final int[][][] triadInOutDegs = new int[][][]{	// 16 triads: 	
		{{0,0}, {0,0},{0,0}},	//1 null 003
		{{0,1}, {1,0},{0,0}},	//2 012
		{{1,1}, {1,1},{0,0}},	//3 102
		{{0,2}, {1,0},{1,0}},	//4 021D
		{{2,0}, {0,1}, {0,1}},	//5 021U
		{{1,1}, {0,1}, {1,0}},	//6 021c
		{{0,1}, {2,1},{1,1}},	//7 111D
		{{1,0}, {1,2},{1,1}},	//8 111U
		{{2,0}, {0,2}, {1,1}},	//9 030T
		{{1,1}, {1,1},{1,1}},	//10 030C
		{{1,1}, {2, 2},{1,1}},	//11 201
		{{0,2}, {2,1},{2,1}},	//12 120D
		{{2,0}, {1,2}, {1,2}},	//13 120U
		{{1,1}, {1,2},{2,1}},	//14 120C
		{{2,1}, {1,2}, {2,2}},	//15 210
		{{2,2}, {2,2},{2,2}}	//16 300
	};
	private static final int[][] triadNodePermutation = new int[][]{
		{0},	// 1. only 1 possible combination {0, 1,2} for type 1{ (0,0) (0,0) (0,0)}, which is triadPermutation[0], (0,0) (0,0) (0,0)
		{0, 1, 2, 3, 4, 5}, // 2. 6 possible combination, see triadPermutation 
		{0, 3, 4},	// 3. 3 possible combination, see triadPermutation 
		{0, 3, 4},	// 4. 3 6 possible combination, see triadPermutation 
		{0, 3, 4},	// 5. 3 possible combination
		{0, 1, 2, 3, 4, 5},	// 6. 6 possible combination, see triadPermutation
		{0, 1, 2, 3, 4 ,5},	// 7. 6 possible combination, see triadPermutation
		{0, 1, 2, 3, 4, 5},	// 8. 6 possible combination, see triadPermutation
		{0, 1, 2, 3, 4, 5},	// 9. 6 possible combination, see triadPermutation
		{0, 1},	// 10. 2 possible combination, see triadPermutation
		{0, 3, 4},	// 11. only 1 possible combination {0, 1,2}, which is tirad permutation 0 (0,0) (0,0) (0,0)
		{0, 3, 4},	// 12. only 1 possible combination {0, 1,2}, which is tirad permutation 0 (0,0) (0,0) (0,0)
		{0, 3, 4},	// 13. only 1 possible combination {0, 1,2}, which is tirad permutation 0 (0,0) (0,0) (0,0)
		{0, 1, 2, 3, 4, 5},	// 14. 3 possible combination, see triadPermutation
		{0, 1, 2, 3, 4, 5},	// 15. 3 possible combination, see triadPermutation
		{0},	// 16. only 1 possible combination {0, 1,2}, which is tirad permutation 0 (0,0) (0,0) (0,0)
		
	};
	private static final int[][] triadPermutation = new int[][]{
		{0, 1, 2},
		{0, 2, 1},
		{1, 0, 2},
		{1, 2, 0},
		{2, 0, 1},
		{2, 1 ,0}
	};
	private static final int[] NumEdgeInTriads = new int[]{
			0, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5 ,6
	};
	
	/**
	 * computer Triad Freq for a random graph whose JOINT in/out degree frequencies is known.
	 * @param in	(int[] )joint in-degree
	 * @param out	(int[] )joint out-degree
	 * @param cnt	(int[] )frequency of join degree
	 * @param size	int		size of the graph.
	 * @return
	 */
	public static double[] getTriadFreqRandomGraphFromJoinInOutDegFreq(int[] in, int[] out, int[] cnt, int size){
		double[] dist = new double[16], tmpDist = new double[16];
		System.out.println("\ncalulating triadFreq\n");
		int numEdge =0;
		for(int i = 0; i<in.length; i++) numEdge += in[i] * cnt[i];
		//TODO: the logBase might be dependent on the motif type
		double logBase = Math.log(numEdge);
		double[][][] logProbFactor = new double[in.length][3][16];
		setLogProbFactor(logProbFactor, in, out, numEdge);
		/*	using null edge cause the probability to decrease by factor of O(N^2) when a triad has one more edge, when |E| << O(n^2), this cause inaccuracy. as a result abandon this method
		for(int d = 0; d < in.length; d++){	// compute the freq that a node with in degree in[d] and out degree oud[d], will have in degree triadInOutDegs[t][u][0] and out degree triadInOutDegs[t][u][1] in a triad. total frequency is (N-1)^2(N-2)^2
			for(int u=0; u<3; u++){	//The u-th node in triad,
				for(int t = 0; t<16; t++){	//the t-th type of triad
					logProbFactor[d][u][t] = getLogProbOfInOutTriadWNodeJoinInOutDeg(in[d], out[d], triadInOutDegs[t][u][0], triadInOutDegs[t][u][1],size);
				}
			}
		}
		*/
		int[] nIDs = new int[3];
		int sumTriadInOutDeg =0;
		double alpha = 0, log6 = Math.log(6), log2 = Math.log(2);
		for(int i=0;i <out.length; i++){	// first node 
			nIDs[0] = i;
			for(int j=i; j<out.length; j++){	// second node
				nIDs[1] = j;
				//isomorphic graph: check duplicate with j=i and make i node duplicate to 2
				if(j == i) {
					if(cnt[i] < 2) continue;
				}
				for(int k = j; k<out.length; k++){	//third node
					nIDs[2] = k;
					//isomorphic graph check duplicate with k=i and make i node duplicate to 3; else check i!=j && j=k, make j node duplicate to 2 and  calculate alpha
					if(k == i){
						if(cnt[i] <3) continue;
						alpha = 0;
						for(int l =0; l<3; l++) alpha += Math.log(cnt[i]-l);
						alpha -= log6;
//						duplicateContribution(logProbContributionOfNodeToTriadNode, 0, 2);
					}else if(i==j && j !=k){
						alpha  = Math.log(cnt[i]) + Math.log(cnt[i]-1) + Math.log(cnt[k]) - log2;
					}else if(i!=j && j == k){
						if(cnt[j] <2) continue;
						alpha = Math.log(cnt[i]) + Math.log(cnt[j]) + Math.log(cnt[j] -1) - log2;
//						duplicateContribution(logProbContributionOfNodeToTriadNode, 1, 2);
					}else{
						alpha = Math.log(cnt[i]) + Math.log(cnt[j]) + Math.log(cnt[k]);
					}
					sumTriadInOutDeg = in[k] + out[k] +in[j] + out[j] + in[i] + out[i];
					// 6 permutation 
					Arrays.fill(tmpDist, 0);
					double tmp1 =0, tmpSum = 0;	// store
					for(int t =0; t<16; t++){	// add expecte number of freq to triad
						//check if other nodes in the network have in degree large enough to receive out-going edge from i,j,k; also check if the sum of out-degree is large enough to than number of in-coming edges to i,j,k,
						//	example: a graph{a->b, a->c, d} extract in/out degree; your method also consider  {a->b, c, a->d} as valid random graph, however, d has 0/0 in/out degree, and that is wrong.  
						//		(a, b and c form a valid triad t --> Oa+Ob+Oc + Ia+Ib+Ic <= |E| + numberOfEdgeInTriad_t 
						if(numEdge + NumEdgeInTriads[t]< sumTriadInOutDeg){
							continue;
						}
						for(int p: triadNodePermutation[t]){
							tmp1 = 0;
							//per = triadPermutation[p];
							for(int l =0; l<3; l++){
								if(logProbFactor[nIDs[l]][triadPermutation[p][l]][t] == Double.NEGATIVE_INFINITY){
									tmp1 = Double.NEGATIVE_INFINITY;
									break;
								}
								tmp1 += logProbFactor[nIDs[l]][triadPermutation[p][l]][t];
							}
							if(tmp1 != Double.NEGATIVE_INFINITY){
								tmpDist[t] +=  Math.exp( tmp1 - NumEdgeInTriads[t]* logBase);
							}
						}
					}
					for(double v: tmpDist) tmpSum += v;
					System.out.printf("\n(%d,%d)*%d, (%d,%d)*%d, (%d,%d)*%d:\n", in[i],out[i],cnt[i], in[j],out[j],cnt[j], in[k], out[k],cnt[k]);
					for(int t=0; t<16; t++){
						System.out.printf("%.5f, ", tmpDist[t]);
						dist[t] += Math.exp(alpha) * tmpDist[t] / tmpSum;
					}
					System.out.println();
				}
			}
		}
		return dist;
	}
	/**
	 * 
	 * @param probMtr a 3*3*16 matrix, probMtr[nIdx][k][t] represent the contribution of log probability of node with nIdx 
	 * 					as the k-node in the t-th triad.
	 * @param triadType
	 * @param nIdx
	 * @param nOut
	 * @param nIn
	 */
	private static void setContributionOfNodeAsTriadNodeToLogProbability(double[][][] probMtr, int triadType, int nIdx, int nIn, int nOut, int size){
		for(int tnIdx = 0; tnIdx < 3; tnIdx++){
			//TODO: improvement check triadInOutDegs[triadType][tnIdx] is duplicated, then just copy previous tnIdx-1,
			probMtr[nIdx][tnIdx][triadType] = getLogProbOfInOutTriadWNodeJoinInOutDeg(nIn, nOut, triadInOutDegs[triadType][tnIdx][0], triadInOutDegs[triadType][tnIdx][1], size);
		}
	}
	private static void duplicateContribution(double[][][] probMtr, int idx1, int idx2){
		for(int i=0; i<3; i++){
			for(int t=0; t<16; t++) probMtr[idx2][i][t] = probMtr[idx1][i][t];
		}
	}
	/**
	 * given a node with join in/out degree and the node in a triad of type ?? has m' outgong edge and n' incomping edge, 
	 * compute the number of null edges and the in/out degree of null edges, then
	 * compute the node's contribution to the probability of this triad.
	 * @param nOut outdegree of a node
	 * @param nIn in degree of node
	 * @param triadOut outdegree of a node in triad ??
	 * @param triadIn in degree of a node in triad
	 * @param size size of the graph.
	 * @return a log of the contribution = \prod_{i=0}^{triadOut-1}(nOut -i)*\prod_{j=0}^{triadIn-1}(nIn-j) 
	 * 										*\prod_{k=0}{nullTriadOut-1}(nullNodeOut-k)*\prod...
	 * where nullTriadOut is the number of null edges going out from the triad node.
	 */
	private static double getLogProbOfInOutTriadWNodeJoinInOutDeg(int nIn, int nOut, int triadIn, int triadOut, int size){
		double res = 0;
		if(nOut<triadOut || nIn < triadIn) return Double.NEGATIVE_INFINITY;
		int nullNOut = size - 1 - nOut, nullTOut = 2 - triadOut;
		int nullNIn = size - 1 - nIn, nullTIn = 2 - triadIn;
		if(nullNOut < nullTOut || nullNIn < nullTIn) return Double.NEGATIVE_INFINITY;
		for(int i=0; i<triadOut; i++) res += Math.log(nOut - i);
		for(int i=0; i<triadIn; i++) res += Math.log(nIn - i);
		for(int i =0; i<nullTOut; i++) res += Math.log(nullNOut - i);
		for(int i=0; i<nullTIn; i++) res += Math.log(nullNIn - i);
		return res;
	}

	
	/**
	 * This ignore null edges, assuming there are too many null edges such that the probability of a null edge is close to 1
	 * @param nIn
	 * @param nOut
	 * @param triadIn
	 * @param triadOut
	 * @param size
	 * @return
	 */
	private static double getLogProbOfInOutTriadWNodeJoinInOutDeg2(int nIn, int nOut, int triadIn, int triadOut, int numEdge){
		double res = 0;
		if(nOut<triadOut || nIn < triadIn) return Double.NEGATIVE_INFINITY;
		for(int i=0; i<triadOut; i++) res += Math.log(nOut - i);
		for(int i=0; i<triadIn; i++) res += Math.log(nIn - i);
		return res;
	}
	public static void setLogProbFactor(double[][][] logProbFactor, int[] in, int[] out, int numEdge){		
		for(int d = 0; d < in.length; d++){	// compute the freq that a node with in degree in[d] and out degree oud[d], will have in degree triadInOutDegs[t][u][0] and out degree triadInOutDegs[t][u][1] in a triad. total frequency is (N-1)^2(N-2)^2
			for(int u=0; u<3; u++){	//The u-th node in triad,
				for(int t = 0; t<16; t++){	//the t-th type of triad
					logProbFactor[d][u][t] = getLogProbOfInOutTriadWNodeJoinInOutDeg2(in[d], out[d], triadInOutDegs[t][u][0], triadInOutDegs[t][u][1],numEdge);
				}
			}
		}
	}	

	
	public static double[] getTriadFreqFromConfigurationModel(int[] in, int out[], int[] cnt, int size, int opt){
		double[] res = new double[16], tmpDist = new double[16];
		int numEdge = 0;
		for(int i=0; i<in.length; i++){
			numEdge += in[i] * cnt[i];
		}
		int[] nIDs = new int[3];
		StructProbFromDegreeEdgeCalculator cal = null;
		if(opt == 4){	//configuration model allowing loop and multi-edges
			cal = new CalculatorTriadStructProbConfigModelIncreaseNullEdgeProb();
		}else if(opt == 1){	//allow loop and multi-edges, but probability of null edge are approximated by 1
			cal = new CalculatorTriadStructProbConfigModelNullEdgeProbIsOne();
		}else if(opt == 2){	//do not allow loop but allow multi-edge
			cal = new CalculatorTriadStructProbConfigModelNoLoop();
		}else if(opt == 3){	//	do not use first order approximation for null edge
			cal = new CalculatorTriadStructProbConfigModelNoFirstOrderApprox();
		}else {
			cal = new CalculatorTriadStructProbConfigModel();
		}
		int sumTriadInOutDeg =0;
		double alpha = 0, distSum = 0, log2 = Math.log(2), log3 = Math.log(3), log6 = Math.log(6), tmpVal = 0;
		for(int i=0; i<in.length; i++){	//first node
			nIDs[0] = i;
			for(int j = i; j<in.length; j++){	//seconde node
				if(j == i && cnt[j] < 2) continue;
				nIDs[1] = j;
				for(int k = j; k < in.length; k++){	//third node
					if(k == i && cnt[i] < 3) continue;
					else if(k == j && cnt[j] < 2) continue;
					nIDs[2] = k;
					//calculate combination / num of 3-node-graphs with identical in/out degrees
					alpha = Math.log(cnt[i]);
					if(j == i) alpha += Math.log(cnt[i] -1) - log2;
					else alpha += Math.log(cnt[j]);
					if(k == i) alpha += Math.log(cnt[i] - 2) - log3;
					else if(k == j) alpha += Math.log(cnt[j] - 1) - log2;
					else alpha += Math.log(cnt[k]);
					//cal distribution
					sumTriadInOutDeg = in[k] + out[k] +in[j] + out[j] + in[i] + out[i];
					Arrays.fill(tmpDist, 0);
//					System.out.printf("\n(%d,%d)*%d\t(%d,%d)*%d\t(%d,%d)*%d: %.4f\n", in[nIDs[0]], out[nIDs[0]], cnt[nIDs[0]], in[nIDs[1]], out[nIDs[1]], cnt[nIDs[1]], in[nIDs[2]], out[nIDs[2]], cnt[nIDs[2]], Math.exp(alpha));
					for(int t = 0; t< 16; t++){
						//check if other nodes in the network have in degree large enough to receive out-going edge from i,j,k; also check if the sum of out-degree is large enough to than number of in-coming edges to i,j,k,
						//	example: a graph{a->b, a->c, d} extract in/out degree; your method also consider  {a->b, c, a->d} as valid random graph, however, d has 0/0 in/out degree, and that is wrong.  
						//		(a, b and c form a valid triad t --> Oa+Ob+Oc + Ia+Ib+Ic <= |E| + numberOfEdgeInTriad_t 
						if(numEdge + NumEdgeInTriads[t]< sumTriadInOutDeg){
							continue;
						}
						tmpDist[t] = TriadStructProbability.getApproxProbConfigModel(nIDs, in, out, t, numEdge, cal);
					}
					distSum = 0;
					for(double d: tmpDist) distSum += d;
					
					for(int t = 0; t<16; t++){
//						System.out.printf("%.4f, ", tmpDist[t]/distSum);
						res[t] += Math.exp(alpha)* tmpDist[t]/distSum;
					}
				}				
			}
		}
		return res;
	}
	
	
	
	/*
	 * 
		//------configuration model for 3-node subgraph
	private static double[][] logInDegSeq;
	private static double[][] logOutDegSeq;
	private static double[] logNumEdge;
	
	private static void setLogInOutDegSeq(int[] in, int[] out, int[] cnt){
		logInDegSeq = new double[in.length][2];
		logOutDegSeq = new double[in.length][2];
		logNumEdge = new double[6];
		long numEdge = 0;
		for(int i =0; i< in.length; i++){
			logInDegSeq[i][0] = in[i]>0?Math.log(in[i]): Double.NEGATIVE_INFINITY;
			logInDegSeq[i][1] = in[i]>1?Math.log(in[i]-1): Double.NEGATIVE_INFINITY;
			logOutDegSeq[i][0] = out[i]>0?Math.log(out[i]): Double.NEGATIVE_INFINITY;
			logOutDegSeq[i][1] = out[i]>1?Math.log(out[i]-1): Double.NEGATIVE_INFINITY;
			numEdge += in[i]* cnt[i];
		}
		for(int i = 0; i< 6; i++){
			logNumEdge[i] = Math.log(numEdge - i);
		}
	}
	 
	 
	 
	public static double[] getTriadFreqFromConfigurationModel(int[] in, int out[], int[] cnt, int size){
		double[] res = new double[16], tmpDist = new double[16];
		setLogInOutDegSeq(in, out, cnt);	//in out degree to log, numEdge to log stored in static variables
		int numEdge = (int) Math.exp(logNumEdge[0]);
		int[] cacheIn = new int[3], cacheOut = new int[3];
		int[] tmpNodes = new int[3];
		int[] nIDs = new int[3];
		int sumTriadInOutDeg =0;
		double alpha = 0, distSum = 0, tmpVal=0;	//, log2 = Math.log(2), log3 = Math.log(3), log6 = Math.log(6), tmpVal = 0;
		for(int i=0; i<in.length; i++){	//first node
			nIDs[0] = i;
			for(int j = i; j<in.length; j++){	//seconde node
				if(j == i && cnt[j] < 2) continue;
				nIDs[1] = j;
				for(int k = j; k < in.length; k++){	//third node
					if(k == i && cnt[i] < 3) continue;
					if(k == j && cnt[j] < 2) continue;
					nIDs[2] = k;
					//calculate combination / num of 3-node-graphs with identical in/out degrees
					alpha = (cnt[i]);
					if(j == i) alpha = alpha * (cnt[i] -1)/ 2;
					else alpha *= (cnt[j]);
					if(k == i) alpha = alpha * (cnt[i] - 2) / 3;
					else if(k == j) alpha = alpha *(cnt[j] - 1) / 2;
					else alpha *= (cnt[k]);
					//cal distribution
					sumTriadInOutDeg = in[k] + out[k] +in[j] + out[j] + in[i] + out[i];
					Arrays.fill(tmpDist, 0);
					for(int t = 0; t< 16; t++){
						//check if other nodes in the network have in degree large enough to receive out-going edge from i,j,k; also check if the sum of out-degree is large enough to than number of in-coming edges to i,j,k,
						//	example: a graph{a->b, a->c, d} extract in/out degree; your method also consider  {a->b, c, a->d} as valid random graph, however, d has 0/0 in/out degree, and that is wrong.  
						//		(a, b and c form a valid triad t --> Oa+Ob+Oc + Ia+Ib+Ic <= |E| + numberOfEdgeInTriad_t 
						if(numEdge + NumEdgeInTriads[t]< sumTriadInOutDeg){
							continue;
						}
						//tmpDist[t] = TriadStructProbability.getApproxProbConfigModel(nIDs, in, out, t, numEdge);
						
						for(int idx: triadNodePermutation[t]){	// each motif type has one or more isomorphic graphs
							for(int m = 0; m<3; m++){
								tmpNodes[m] = nIDs[triadPermutation[idx][m]];
							}
								//v1. get probablity by assuming O_u*I_v/|E| <=1 (could has problem leading to less accuracy)
								//tmpVal = getProbForMotifTypeWNodes(t, tmpNodes, in, out, numEdge);
								//v2. get Probability and make sure probability all edge prob is less than 1
								tmpVal = getProb3NodeStructureWithCacheIODeg(t, tmpNodes, in, out, cacheIn, cacheOut, numEdge);
								if(tmpVal == Double.NEGATIVE_INFINITY) continue;
								tmpDist[t] += tmpVal;
							
						}
						
					}
					distSum = 0;
					for(double d: tmpDist) distSum += d;
					System.out.printf("\n(%d,%d)*%d\t(%d,%d)*%d\t(%d,%d)*%d: %.4f\n", in[nIDs[0]], out[nIDs[0]], cnt[nIDs[0]], in[nIDs[1]], out[nIDs[1]], cnt[nIDs[1]], in[nIDs[2]], out[nIDs[2]], cnt[nIDs[2]], alpha);
					for(int t = 0; t<16; t++){
						System.out.printf("%.4f, ", tmpDist[t]);
						res[t] += (alpha)* tmpDist[t]/distSum;
					}
				}				
			}
		}
		return res;
	}
	*/
	/*
	//---------------------------------------------
	// configuration model for directed graph [V2]: make sure O_uI_v/|E| <=1
	//---------------------------------------------
	final static int[][] DIREDGE = {	// {0, 1} represents node 0 -> node 1
		{0, 1}, {0, 2}, {1, 0}, {1, 2}, {2, 0}, {2,1} 
	};
	final static int[][][] TRIEDGE = {
		{ 	{},	//existing edge for type 0
				{0, 1, 2, 3, 4, 5}},	//null edge for type 1: TRIEDGE[0][1][i]= k represent triad type =0, has DIREDGE[k] as null edge
		{ 	{0},	//existing edge for type 1:  	TRIEDGE[1][0][i]= k represent triad type =1, has DIREDGE[k] as existing edge
				{1, 2, 3, 4, 5}},
		{ 	{0, 2},	//existing edge for type 2
				{1, 3, 4, 5}},
		{ 	{0, 1},	//existing edge for type 3
				{ 2, 3, 4, 5}},
		{ 	{2, 4},	//existing edge for type 4
				{0, 1, 3, 5}},
		{ 	{1, 2}, // existing edge for type 5
				{ 0, 3, 4, 5 } },
		{ 	{0, 3, 5}, // existing edge for type 6
				{1, 2, 4 } },
		{ 	{2, 3, 5}, // existing edge for type 7
				{ 0, 1, 4} },
		{ 	{2, 3, 4}, // existing edge for type 8
				{ 0, 1, 5 } },
		{ 	{0, 3, 4}, // existing edge for type 9
				{1, 2, 5 } },
		{ 	{0, 2, 3, 5}, // existing edge for type 10
				{1,4} },
		{ 	{0, 1, 3, 5}, // existing edge for type 11
				{ 2, 4} },
		{ 	{2, 3, 4, 5}, // existing edge for type 12
				{ 0, 1} },
		{ 	{1, 2, 3, 5}, // existing edge for type 13
				{0,4 } },
		{ 	{1, 2, 3, 4, 5 }, // existing edge for type 14
				{ 0} },
		{ 	{0, 1, 2, 3, 4, 5}, // existing edge for type 15
				{  } },
	};/*
	/**
	 * calculate the probability given a 3-node graph structure given the io degree of the nodes in a configuration model 
	 * @param type
	 * @param inDeg
	 * @param outDeg
	 * @param numEdge
	 * @return
	 */
	/*
	private static double calProb3NodeStructure(int type, int[] inDeg, int[] outDeg, int numEdge){
		for(int i = 0; i< 3 ; i++){
			if(inDeg[i] < triadInOutDegs[type][i][0] || outDeg[i] < triadInOutDegs[type][i][1])
				//if in out degree of nID[i] do not satisfy the i-th noe in triad
				return 0;
		}
		double res = 1;
		//calculate existing edges
		for(int i: TRIEDGE[type][0]){	// i is DIREDGE idx
			res *=  Math.min(1, 1.0 * outDeg[DIREDGE[i][0]] * inDeg[DIREDGE[i][1]] / numEdge);
			outDeg[DIREDGE[i][0]]--;
			inDeg[DIREDGE[i][1]]--;
			numEdge--;
		}
		//calculate null edges
		for(int i: TRIEDGE[type][1]){
			res *= Math.max(numEdge - outDeg[DIREDGE[i][0]] * inDeg[DIREDGE[i][1]], 0.0)/numEdge;
		}
		return res;
	}
	
	private static double getProb3NodeStructureWithCacheIODeg(int type, int[] nIDs, int[] in, int[] out, int[] cacheIn, int[] cacheOut, int numEdge){
		for(int i=0; i< 3; i++){
			cacheIn[i] = in[nIDs[i]];
			cacheOut[i] = out[nIDs[i]];
		}
		return calProb3NodeStructure(type, cacheIn, cacheOut, numEdge);
	}*/
}
