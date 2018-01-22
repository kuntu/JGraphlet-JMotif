package motifs;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import mathFunctions.*;

public class TriadTransition {
	public static int[] tmpCensus;
	public static int[] tmpEditDistanceDistr;
	public static int[] curEditDistanceDistr;
	public static int[][] allEditDistanceFreq;
	public static int[][] randomGraphTriadDistr;
	public static ArrayList<Double> KLDTriadCensus;
	public static ArrayList<Double> KLDTriadEditDisDistr;
	
	public int[][][] edgeGraph;
	public String dataName;
	public String outputDir;
	public int time;
	public int[][] triadFreq;
	public int[][][] transitionCnt;
	public HashMap<Set<Integer>, List<MotifChainNode>> motifChains;
	HashMap<Set<Integer>, StatSubGraphChange> allStatSubGraph;
	public int numNodes;
	
	public TriadTransition(int[][][] eg, int numNode){
		edgeGraph = eg;
		numNodes = numNode;
		if(edgeGraph!=null){
			time = edgeGraph.length;
			triadFreq = new int[time][16];
			transitionCnt = new int[time][16][16];
			motifChains = new HashMap<Set<Integer>, List<MotifChainNode>>();
			dataName = "";
			outputDir = "";
		}
	}
	
	@SuppressWarnings("finally")
	public static TriadTransition getDynNetFromEdgeFile(String fileName){
		TriadTransition tt = null;
		try{
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			System.out.println("\t read:"+fileName);
			String line = null;
			Pattern p = Pattern.compile("-?\\d+");
			Matcher m = null;
			int numN = 0;
			int t = 0;
			ArrayList<Integer> ls = new ArrayList<Integer>();
			//read basic information: the first line contains 1. number of nodes; 2. number of time steps; 
			if((line = br.readLine())!= null){
				m = p.matcher(line);
				if(m.find()) numN = Integer.parseInt(m.group());
				if(m.find()) t = Integer.parseInt(m.group());
			}else{
				br.close();
				return null;
			}
			int[][][] edges = new int[t][][];
			//read edges
			t = 0;
			while((line = br.readLine()) !=null){
				m = p.matcher(line);
				while(m.find()){
					ls.add(Integer.parseInt(m.group()));
				}
				edges[t] = new int[ls.size()/2][2];
				for(int i=0; i<edges[t].length; i++){
					edges[t][i][0] = ls.get(2*i);
					edges[t][i][1] = ls.get(2*i+1);
				}
				t++;
				ls.clear();
			}
			while(t<edges.length) edges[t++] = new int[0][2];
			
			tt = new TriadTransition(edges, numN);
			br.close();
		}catch(Exception e){
			e.printStackTrace();
		}finally{
			return tt;
		}
	}
	//----------------- 
	//    output files
	//-----------------
	
	public void outputTriadCensus(){
		try{
			String triadCensusFile = outputDir + dataName + "TriadFreqCensus.csv";
			BufferedWriter bw = new BufferedWriter(new FileWriter(triadCensusFile));
			//System.out.println("\twrite triad census to file: \n\t"+ triadCensusFile);
			StringBuilder sb = new StringBuilder();
			for(int[] freq: triadFreq){
				sb.setLength(0);
				for(int f: freq) sb.append(f+"\t");
				sb.setLength(sb.length()-1);
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void outputRandGraphSameAsymDyadTriadFreq(){
		try{
			String randGraphTriadFile = outputDir + dataName + "randomGraphSameAsymDyadFreqs.csv";
			BufferedWriter bw = new BufferedWriter(new FileWriter(randGraphTriadFile));
			//System.out.println("\tcompute random graph with the same dyads...\n\twrite to file: " + randGraphTriadFile);
			double[] randGraphTriadDistr = null;
			StringBuilder sb = new StringBuilder();
			for(int[][] graph: edgeGraph){
				sb.setLength(0);
				MotifGraph mg = new MotifGraph(numNodes, graph, true);
				randGraphTriadDistr = mg.getTriadCountFromRandomGraphWSameAsymDyads();
				for(double d: randGraphTriadDistr){
					sb.append(d + "\t");
				}
				sb.setLength(sb.length() -1);
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void outputRandGraphSameDyadTriadFreq(){
		try{
			String randGraphTriadFile = outputDir + dataName + "randomGraphSameDyadFreqs.csv";
			BufferedWriter bw = new BufferedWriter(new FileWriter(randGraphTriadFile));
			//System.out.println("\tcompute random graph with the same dyads...\n\twrite to file: " + randGraphTriadFile);
			double[] randGraphTriadDistr = null;
			StringBuilder sb = new StringBuilder();
			for(int[][] graph: edgeGraph){
				sb.setLength(0);
				MotifGraph mg = new MotifGraph(numNodes, graph, true);
				randGraphTriadDistr = mg.getTriadCountFromRandomGraphWSameDyads();
				for(double d: randGraphTriadDistr){
					sb.append(d + "\t");
				}
				sb.setLength(sb.length() -1);
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void outputRandGraphSameDyadTriadDistr(){
		try{
			String randGraphTriadFile = outputDir + dataName + "randomGraphSameDyadDistrs.csv";
			BufferedWriter bw = new BufferedWriter(new FileWriter(randGraphTriadFile));
			//System.out.println("\tcompute random graph with the same dyads...\n\twrite to file: " + randGraphTriadFile);
			double[] randGraphTriadDistr = null;
			StringBuilder sb = new StringBuilder();
			for(int[][] graph: edgeGraph){
				sb.setLength(0);
				MotifGraph mg = new MotifGraph(numNodes, graph, true);
				randGraphTriadDistr = mg.getTriadDistrFromRandomGraphWSameDyads();
				for(double d: randGraphTriadDistr){
					sb.append(d + "\t");
				}
				sb.setLength(sb.length() -1);
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * output the frequency of transition between triads
	 */
	public void outputAllTransitionMatrices(){
		try{
			String tranMtrFile = outputDir + dataName+"triadTransitionCountMatrix.csv";
			BufferedWriter bw = new BufferedWriter(new FileWriter(tranMtrFile));
			for(int i=1; i< transitionCnt.length; i++){
				for(int[] tran: transitionCnt[i]){
					for(int n: tran) bw.write(n+"\t");
				}
				bw.newLine();
			}
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void outputSplitedGraphs(int[][] mappings, int[] subgraphSize){
		String[] filenames = new String[subgraphSize.length];
		BufferedWriter[] bws = new BufferedWriter[subgraphSize.length];
		try{
			for(int i=0; i<filenames.length; i++){
				filenames[i] = outputDir +"subGraph" + i+".txt";
				if(subgraphSize[i]>10){
					bws[i] = new BufferedWriter(new FileWriter(filenames[i]));
					bws[i].write(subgraphSize[i] +"\t"+ time);
					bws[i].newLine();
				}
			}
			for(int[][] graph: edgeGraph){
				for(int[] edge: graph){
					if(bws[mappings[edge[0]][0]]==null) continue;
					bws[mappings[edge[0]][0]].write(mappings[edge[0]][1] +"\t" );
					bws[mappings[edge[1]][0]].write(mappings[edge[1]][1] +"\t");
				}
				for(int i=0; i<bws.length; i++){
					if(bws[i]==null) continue;
					bws[i].newLine();
				}
			}
			for(int i=0; i<filenames.length; i++){
				if(bws[i]==null) continue;
				bws[i].close();
			}
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * output a file containing a vector of the KL-divergence of triad census throughout time
	 * 
	 */
	public void outputKLDivergenceOfTriadCensus(){
		try{
			String tranMtrFile = outputDir + dataName+"triadCensusKLDivergence.csv";
			BufferedWriter bw = new BufferedWriter(new FileWriter(tranMtrFile));
			for(double kld: KLDTriadCensus){
				bw.write(kld+"\t");
			}
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void outputKLDivergenceOfEditDistance(){
		try{
			String file = outputDir + dataName+"triadEditDistanceKLDivergence.csv";
			BufferedWriter bw = new BufferedWriter(new FileWriter(file));
			for(double kld: KLDTriadEditDisDistr){
				bw.write(kld+"\t");
			}
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void outputAllEditDistanceDistribution(){
		try{
			String file = outputDir + dataName+"triadAllEditDistanceFreq.csv";
			BufferedWriter bw = new BufferedWriter(new FileWriter(file));
			for(int[] dist: allEditDistanceFreq){
				for(int i: dist)
					bw.write(i+"\t");
				bw.newLine();
			}
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void outputMotifChainAll3NodeSubGraphs(){
		try{
			String subGraphFile = outputDir + dataName + "3NodeSubgraphMotifChain.csv";
			BufferedWriter br = new BufferedWriter(new FileWriter(subGraphFile));
			StatSubGraphChange ssgc = null;
			StringBuilder sb = new StringBuilder();
			for(Set<Integer> s: allStatSubGraph.keySet()){
				ssgc = allStatSubGraph.get(s);
				for(int i: s) {
					sb.append(i);
					sb.append('\t');
				}
				sb.append(ssgc.nonNullTime);
				for(int i: ssgc.editDisDistr){
					sb.append('\t');
					sb.append(i);
				}
				br.write(sb.toString());
				br.newLine();
				sb.setLength(0);
			}
			br.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public StatSubGraphChange getStatOfSubgraph(List<MotifChainNode> ar){
		StatSubGraphChange res = new StatSubGraphChange();
		StatSubGraphChange.nonNullBeg = -1;
		StatSubGraphChange.nonNullEnd = -1;
		MotifChainNode pre = null;
		for(MotifChainNode mcn: ar){
			res.calChange(pre, mcn);
			pre = mcn; 
		}
		if(pre!=null&& StatSubGraphChange.nonNullEnd<StatSubGraphChange.nonNullBeg && pre.motifType !=0) 
			res.nonNullTime += time - StatSubGraphChange.nonNullBeg;
		return res;
	}
	
	public HashMap<Set<Integer>, StatSubGraphChange> getAllSubGraphChanges(){
		HashMap<Set<Integer>, StatSubGraphChange> res = new HashMap<Set<Integer>, StatSubGraphChange>();
		StatSubGraphChange.totalTime = time;
		for(Set<Integer> key: motifChains.keySet()){
			res.put(key, getStatOfSubgraph(motifChains.get(key)));
		}
		this.allStatSubGraph = res;
		return res;
	}
	
	//-------- static method for computing MotifGraphs
	public static void initialTmpVar(){

		if(KLDTriadCensus==null) KLDTriadCensus = new ArrayList<Double>();
		else KLDTriadCensus.clear();
		if(KLDTriadEditDisDistr==null) KLDTriadEditDisDistr = new ArrayList<Double>();
		else KLDTriadEditDisDistr.clear();
	}
	
	/**
	 * Count Triads for an edge network; store subgraphs that have at lease one edge.
	 * @param graph
	 * @return
	 */
	public static HashMap<Set<Integer>, List<MotifChainNode>> countTriads(MotifGraph graph, int[] census, int time){
		HashMap<Set<Integer>, List<MotifChainNode>> res = new HashMap<Set<Integer>, List<MotifChainNode>>();
		if(census == null || census.length!=16){
			System.out.println("census length should be 16 for triads");
			return res;
		}
		int[] subgraph = new int[3];
		int[] iodegrees = new int[3];
		int triadHashKey = -1;
		int motifType = -1;
		Set<Integer> subGraphKey = null;
		ArrayList<MotifChainNode> triadChain = null;
		for(long e: graph.edges){
			MotifGraph.edgeKeyToSubgraphArray(e, subgraph, 0, 1);	//assign two nodes to a 3-node subgraph from edage e
			for(subgraph[2] = 0; subgraph[2] < graph.numNode; subgraph[2]++){	// process all 3-node subgraphs containing edge e  
				if(subgraph[2] == subgraph[0] || subgraph[2] == subgraph[1]) continue;
				subGraphKey = MotifGraph.getSubgraphkey(subgraph);
				if(res.containsKey(subGraphKey)) continue;	// the subGraph has already been processed
				//detect triad type for subgraph by obtain all edges and calculate the in/out degree for the three nodes
				triadHashKey = MotifGraph.subGraphToTriadHashKey(subgraph, graph, iodegrees);
				if(!MotifGraph.triadCodeIdx.containsKey(triadHashKey)){
					System.out.println("error! no such hashKey found!");
					return res;
				}
				// For a subgraph that has at least one edge, record its specific motif type
				triadChain = new ArrayList<MotifChainNode>();
				res.put(subGraphKey, triadChain);
				motifType = MotifGraph.triadCodeIdx.get(triadHashKey);
				triadChain.add(new MotifChainNode(time, motifType));
				//count the corresponding motif type
				census[motifType]++;
			}
		}
		//calculate frequency for 003 triad
		motifType = graph.numNode;
		motifType = motifType * (motifType-1) * (motifType-2) / 6;
		for(int i=1; i<census.length; i++) motifType -= census[i];
		census[0] = motifType;
		return res;
	}
	
	/**
	 * update triad information of current nextwork snapshot with info from previous snapshots
	 * @param preGraph
	 * @param curGraph
	 * @param census
	 * @param curTime
	 * @param allMotifChain
	 * @param subgraph
	 * @param iodegrees
	 * @param tranMtr
	 * @return
	 */
	public static HashMap<Set<Integer>, List<MotifChainNode>> updateTriadsWNextNetSnapshot(MotifGraph preGraph, 
			MotifGraph curGraph, int[] census, int curTime, HashMap<Set<Integer>, List<MotifChainNode>> allMotifChain, 
			int[] subgraph, int[] iodegrees, int[][] tranMtr){
		if(census == null || census.length!=16){
			System.out.println("census length should be 16 for triads");
			return allMotifChain;
		}
		if(tranMtr==null){
			System.out.println("transition matrix need to be initialized");
			return allMotifChain;
		}
		tmpCensus = MathFun.cloneCensus(tmpCensus, census);
		if(curEditDistanceDistr==null) curEditDistanceDistr = new int[7];
		else Arrays.fill(curEditDistanceDistr,0);
		
		//HashMap<Set<Integer>, List<MotifChainNode>> res = new HashMap<Set<Integer>, List<MotifChainNode>>();
		HashMap<Set<Integer>, MotifChainNode> res = new HashMap<Set<Integer>, MotifChainNode>();
		if(subgraph==null) subgraph = new int[3];
		if(iodegrees==null) iodegrees = new int[3];
		int triadHashKey = -1;
		int motifType = -1;
		HashSet<Long> changingEdges = null;
		changingEdges = MotifGraph.getChangingEdges(preGraph, curGraph, changingEdges);	// obtain edges appearing/disappearing in current snapshots
		Set<Integer> subGraphKey = null;
		List<MotifChainNode> triadChain = null;
		for( int i=0; i< 16; i++) tranMtr[i][i] = census[i];
		for(long e: changingEdges){
			MotifGraph.edgeKeyToSubgraphArray(e, subgraph, 0, 1);	//assign two nodes to a 3-node subgraph from edage e
			for(subgraph[2] = 0; subgraph[2] < curGraph.numNode; subgraph[2]++){	// process all 3-node subgraphs containing edge e  
				if(subgraph[2] == subgraph[0] || subgraph[2] == subgraph[1]) continue;
				subGraphKey = MotifGraph.getSubgraphkey(subgraph);
				if(res.containsKey(subGraphKey)) continue;	// the subGraph has already been processed
				//detect triad type for subgraph by obtain all edges and calculate the in/out degree for the three nodes
				triadHashKey = MotifGraph.subGraphToTriadHashKey(subgraph, curGraph, iodegrees);
				if(!MotifGraph.triadCodeIdx.containsKey(triadHashKey)){
					System.out.println("error! no such hashKey found!");
					return allMotifChain;
				}
				// For a subgraph that has at least one edge, record its specific motif type				
				motifType = MotifGraph.triadCodeIdx.get(triadHashKey);
				res.put(subGraphKey, new MotifChainNode(curTime, motifType));
				//count the corresponding motif type
				census[motifType]++;				
			}
		}
		// calculate motif transition
		for(Set<Integer> subG: res.keySet()){
			if(allMotifChain.containsKey(subG)){
				triadChain = allMotifChain.get(subG);
				motifType = triadChain.get(triadChain.size()-1).motifType;
				census[motifType]--;
				tranMtr[motifType][motifType]--;
				triadChain.add(res.get(subG));
				tranMtr[motifType][res.get(subG).motifType]++;
				curEditDistanceDistr[MotifGraph.TRIAD_EDIT_DISTANCE[motifType][res.get(subG).motifType]]++;
			}else{	// subG previously has been the triad "003" and was not recorded
				triadChain = new ArrayList<MotifChainNode>();
				triadChain.add(res.get(subG));
				allMotifChain.put(subG, triadChain);
				census[0]--;
				tranMtr[0][0]--;
				tranMtr[0][res.get(subG).motifType]++;
				curEditDistanceDistr[MotifGraph.TRIAD_EDIT_DISTANCE[0][res.get(subG).motifType] ]++;
			}
		}
		motifType = curGraph.numNode;
		motifType = motifType * (motifType-1) * (motifType-2) / 6;
		for(int i=1; i<curEditDistanceDistr.length; i++) motifType -= curEditDistanceDistr[i];
		curEditDistanceDistr[0] = motifType;
		
		// compute the KL-divergence of triadcensus between previus and current step
		KLDTriadCensus.add(MathFun.KLDiversionFromFreqWithBayesPrior(tmpCensus, census));
		if(tmpEditDistanceDistr==null) tmpEditDistanceDistr = new int[7];
		else KLDTriadEditDisDistr.add(MathFun.KLDiversionFromFreqWithBayesPrior(tmpEditDistanceDistr, curEditDistanceDistr));
		tmpEditDistanceDistr = MathFun.cloneCensus(tmpEditDistanceDistr, curEditDistanceDistr);
		MathFun.cloneCensus(allEditDistanceFreq[curTime], curEditDistanceDistr);
		//calculate frequency for 003 triad
//		motifType = curGraph.numNode;
//		motifType = motifType * (motifType-1) * (motifType-2) / 6;
//		for(int i=1; i<census.length; i++) motifType -= census[i];
//		census[0] = motifType;
		return allMotifChain;
	}
	
	public HashMap<Set<Integer>, Integer> updateTriadFreqAtT(int t, MotifGraph pre, MotifGraph cur, HashMap<Set<Integer>, Integer> preSubgraphType){
		HashMap<Set<Integer>, Integer> curSubgraphType = new HashMap<Set<Integer>, Integer>();
		int[] subgraph = new int[3];
		int[] iodegrees = new int[3];
		int triadHashKey = -1;
		int motifType = -1;
		int preType = -1;
		Set<Integer> subGraphKey = null;
		int[] edgeNums = new int[]{cur.edges.size(), pre.edges.size(), 0};
		//assume triad freq is calculated at t-1;
		for(long e: cur.edges){	// deal with new edges
			if(pre.edges.contains(e)){
				pre.edges.remove(e);
			}else{
				MotifGraph.edgeKeyToSubgraphArray(e, subgraph, 0, 1);
				for(subgraph[2] = 0; subgraph[2] < cur.numNode; subgraph[2]++){	// process all 3-node subgraphs containing edge e  
					if(subgraph[2] == subgraph[0] || subgraph[2] == subgraph[1]) continue;
					subGraphKey = MotifGraph.getSubgraphkey(subgraph);
					if(curSubgraphType.containsKey(subGraphKey)) continue;	// the subGraph has already been processed
					//detect triad type for subgraph by obtain all edges and calculate the in/out degree for the three nodes
					triadHashKey = MotifGraph.subGraphToTriadHashKey(subgraph, cur, iodegrees);
					// For a subgraph that has at least one edge, record its specific motif type
					motifType = MotifGraph.triadCodeIdx.get(triadHashKey);
					triadFreq[t][motifType]++;
					curSubgraphType.put(subGraphKey, motifType);
					if(preSubgraphType.containsKey(subGraphKey)) {
						preType = preSubgraphType.get(subGraphKey);
						preSubgraphType.remove(subGraphKey);
					}else preType = 0;
					transitionCnt[t][preType][motifType]++;
					transitionCnt[t][preType][preType]--;
				}
			}
		}
		HashSet<Set<Integer>> nullSubgraph = new HashSet<Set<Integer>>();
		edgeNums[2] = pre.edges.size();
		edgeNums[1] -= pre.edges.size();
			//System.out.println("\t save step: " + (-edgeNums[2]+edgeNums[1]));
		for(long e: pre.edges){//deal with disappearing edges
			MotifGraph.edgeKeyToSubgraphArray(e, subgraph, 0, 1);
			for(subgraph[2] = 0; subgraph[2] < cur.numNode; subgraph[2]++){	// process all 3-node subgraphs containing edge e  
				if(subgraph[2] == subgraph[0] || subgraph[2] == subgraph[1]) continue;
				subGraphKey = MotifGraph.getSubgraphkey(subgraph);
				if(curSubgraphType.containsKey(subGraphKey) || nullSubgraph.contains(subGraphKey)) continue;	// the subGraph has already been processed
				//detect triad type for subgraph by obtain all edges and calculate the in/out degree for the three nodes
				triadHashKey = MotifGraph.subGraphToTriadHashKey(subgraph, cur, iodegrees);
				preType = preSubgraphType.get(subGraphKey);
				motifType = MotifGraph.triadCodeIdx.get(triadHashKey);
				if(motifType != 0){
					curSubgraphType.put(subGraphKey, motifType);
				}else nullSubgraph.add(subGraphKey);
				preSubgraphType.remove(subGraphKey);
				//count the corresponding motif type
				triadFreq[t][motifType]++;
				transitionCnt[t][preType][motifType]++;
				transitionCnt[t][preType][preType]--;
			}
		}
		if(t>0){
			for(int i=0; i<triadFreq[t].length; i++){
				transitionCnt[t][i][i] += triadFreq[t-1][i];
			}
		}else{
			transitionCnt[t][0][0] += (int) (1L* cur.numNode * (cur.numNode-1) /2 * (cur.numNode-2) /3) ;
		}
		//calculate frequency for 003 triad
		motifType = cur.numNode;
		motifType = (int) (1L* cur.numNode * (cur.numNode-1) /2 * (cur.numNode-2) /3);
		for(int i=1; i<triadFreq[t].length; i++) motifType -= triadFreq[t][i];
		triadFreq[t][0] = motifType;
		curSubgraphType.putAll(preSubgraphType);
		return curSubgraphType;
	}
	public HashMap<Set<Integer>, Integer> countTriadFreqAtT(int t, MotifGraph graph, HashMap<Set<Integer>, Integer> preSubgraphType){
		HashMap<Set<Integer>, Integer> curSubgraphType = new HashMap<Set<Integer>, Integer>();
		int[] subgraph = new int[3];
		int[] iodegrees = new int[3];
		int triadHashKey = -1;
		int motifType = -1;
		int preType = -1;
		Set<Integer> subGraphKey = null;
		for(long e: graph.edges){
			MotifGraph.edgeKeyToSubgraphArray(e, subgraph, 0, 1);	//assign two nodes to a 3-node subgraph from edage e
			for(subgraph[2] = 0; subgraph[2] < graph.numNode; subgraph[2]++){	// process all 3-node subgraphs containing edge e  
				if(subgraph[2] == subgraph[0] || subgraph[2] == subgraph[1]) continue;
				subGraphKey = MotifGraph.getSubgraphkey(subgraph);
				if(curSubgraphType.containsKey(subGraphKey)) continue;	// the subGraph has already been processed
				//detect triad type for subgraph by obtain all edges and calculate the in/out degree for the three nodes
				triadHashKey = MotifGraph.subGraphToTriadHashKey(subgraph, graph, iodegrees);
				// For a subgraph that has at least one edge, record its specific motif type
				motifType = MotifGraph.triadCodeIdx.get(triadHashKey);
				curSubgraphType.put(subGraphKey, motifType);
				//count the corresponding motif type
				triadFreq[t][motifType]++;
				if(preSubgraphType.containsKey(subGraphKey)) {
					preType = preSubgraphType.get(subGraphKey);
					preSubgraphType.remove(subGraphKey);
				}else preType = 0;
				transitionCnt[t][preType][motifType]++;
				transitionCnt[t][preType][preType]--;
			}
		}
		//
		for(int val:preSubgraphType.values()) {
			transitionCnt[t][val][0]++;
			transitionCnt[t][val][val]--;
		}
		if(t>0){
			for(int i=0; i<triadFreq[t].length; i++){
				transitionCnt[t][i][i] += triadFreq[t-1][i];
			}
		}else{
			transitionCnt[t][0][0] += (int) (1L* graph.numNode * (graph.numNode-1) /2 * (graph.numNode-2) /3);
		}
		//calculate frequency for 003 triad
		motifType = graph.numNode;
		motifType = motifType * (motifType-1) * (motifType-2) / 6;
		for(int i=1; i<triadFreq[t].length; i++) motifType -= triadFreq[t][i];
		triadFreq[t][0] = motifType;
		return curSubgraphType;
	}
	
	//----------end of static methods
	
	public static void main(String[] args) {
		System.out.println("Working Directory = " + System.getProperty("user.dir"));
		//for(String s: args) System.out.println(s);
		String filename = "test.txt";
		String dataName ="";
		if(args!=null && args.length>0) {
			filename = args[0];
			if(args.length>1) dataName = args[1];
		}
		TriadTransition tt = TriadTransition.getDynNetFromEdgeFile(filename);
		tt.dataName = dataName;
		initialTmpVar();
		TriadTransition.allEditDistanceFreq = new int[tt.time][7];
		MotifGraph.initializeTriadCode();
		MotifGraph pre = new MotifGraph(tt.numNodes, tt.edgeGraph[0], true);
		MotifGraph cur = null;
		
		tt.motifChains = countTriads(pre, tt.triadFreq[0], 0);;
		int[] subgraph = new int[3];
		int[] iodegrees = new int[3];
		for(int t=1; t<tt.time; t++){
			MathFun.cloneCensus(tt.triadFreq[t], tt.triadFreq[t-1]);
			cur = new MotifGraph(tt.numNodes, tt.edgeGraph[t], true);
			updateTriadsWNextNetSnapshot(pre, cur, tt.triadFreq[t], t, tt.motifChains, subgraph,
					iodegrees, tt.transitionCnt[t]);
			pre = cur;
		}
		//------------------
		// output files
		//-------------------
		tt.outputAllTransitionMatrices();
		tt.outputKLDivergenceOfTriadCensus();
		tt.outputKLDivergenceOfEditDistance();
		tt.outputAllEditDistanceDistribution();
//		tt.getAllSubGraphChanges();
//		tt.outputMotifChainAll3NodeSubGraphs();
		System.out.println("done");
		
		return ;
	}

}
