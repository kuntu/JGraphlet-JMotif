import java.io.*;
import java.util.*;

import graphs.*;
import mathFunctions.*;
import motifs.*;
import randomgraph.*;


public class ExperimentTempNetTriad {
	public static void speedTest(String configFile){
		try{
			BufferedReader br = new BufferedReader(new FileReader(configFile));
			String line = null;
			String dataFile = br.readLine();
			TriadTransition tt = TriadTransition.getDynNetFromEdgeFile(dataFile);
			tt.dataName = br.readLine();
			System.out.println("dataset: " + tt.dataName);
			tt.outputDir = br.readLine()+"/"+tt.dataName+"/";
			File outDir = new File(tt.outputDir);
			if(!outDir.exists()) {
				outDir.mkdirs();
				System.out.println("\tcreate: "+ tt.outputDir);
			}
			HashMap<Set<Integer>, Integer> curSubgraphType = new HashMap<Set<Integer>, Integer>();
			MotifGraph[] graphs = new MotifGraph[tt.time];
			for(int i=0; i<tt.time; i++){
				graphs[i] = new MotifGraph(tt.numNodes, tt.edgeGraph[i], true);
			}
			long begTime = System.nanoTime();
			for(int t = 0; t<tt.edgeGraph.length; t++){
				curSubgraphType = tt.countTriadFreqAtT(t, graphs[t], curSubgraphType);
//				for(int i: tt.triadFreq[t]) System.out.print(i+" ");
//				System.out.print(": ");
//				for(int[] trans: tt.transitionCnt[t]){
//					System.out.print("\t");
//					for(int tr: trans) System.out.print(tr+" ");
//				}
//				System.out.println();
			}
			long elapseTime = System.nanoTime() - begTime;
			System.out.println("\tbase line: \t\t" + elapseTime);
			//fast algorithm
			curSubgraphType = new HashMap<Set<Integer>, Integer>();
			tt.transitionCnt = new int[tt.time][16][16];
			tt.triadFreq = new int[tt.time][16];
			begTime = System.nanoTime();
			curSubgraphType = tt.countTriadFreqAtT(0, graphs[0], curSubgraphType);
			for(int t = 1; t<tt.edgeGraph.length; t++){
				curSubgraphType = tt.updateTriadFreqAtT(t, graphs[t-1], graphs[t], curSubgraphType);
//				for(int i: tt.triadFreq[t]) System.out.print(i+" ");
//				System.out.print(": ");
//				for(int[] trans: tt.transitionCnt[t]){
//					System.out.print("\t");
//					for(int tr: trans) System.out.print(tr+" ");
//				}
//				System.out.println();
			}
			elapseTime = System.nanoTime()-begTime;
			System.out.println("\tfast algorithm\t\t" + elapseTime);
			
			//from Pinghui
			GraphNodePairVal[] newGraphs = new GraphNodePairVal[tt.time];
			for(int t = 0; t < tt.time; t++){
				newGraphs[t] = new GraphNodePairVal(tt.numNodes, true, tt.edgeGraph[t]);
			}
			begTime = System.nanoTime();
			for(GraphNodePairVal g: newGraphs){
				g.getMotifFreq(3);
			}
			elapseTime = System.nanoTime() - begTime;
			System.out.println("\tPinghui algorithm\t" + elapseTime);
			br.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static long[][] getMotifDistributionOfTemporalNetFromFile(int motifSize, String fileName){
		GraphNodePairVal g = null;
		int numMotif = 0;
		if(motifSize == 3) numMotif = 16;
		else if(motifSize == 4) numMotif = 199;
		try{
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			String line = br.readLine();
			String[] data = line.split("[ \t]");
			int num = Integer.parseInt(data[0]);
			int T = Integer.parseInt(data[1]);
			int[][] edgeArray = null;
			long[][] motifDis = new long[T][numMotif];
			for(int t = 0; t < T; t++){
				line = br.readLine();
				if(line == null ||line.length()<2){
					motifDis[t] = new long[numMotif];
					motifDis[t][0] = num;
					motifDis[t][0] = motifDis[t][0] * (motifDis[t][0] -1 )/ 2 * (motifDis[t][0] -2) / 3;
					continue;
				}
				data = line.split("[ \t]");
				edgeArray = new int[data.length/2][2];
				for(int i=0; i < edgeArray.length; i++){
					edgeArray[i][0] = Integer.parseInt(data[2*i]);
					edgeArray[i][1] = Integer.parseInt(data[2*i+1]);
				}
				g = new GraphNodePairVal(num, true, edgeArray);
				motifDis[t] = g.getMotifFreq(motifSize);
			}
			br.close();
			return motifDis;
		}catch(Exception e){
			e.printStackTrace();
			return new long[0][0];
		}
	}
	
	/**
	 * run experiment with config file
	 * @param configFile: config file name
	 */
	public static void runExperiment(String configFile){
		//System.out.println("Working Directory = " + System.getProperty("user.dir"));
		try{
			BufferedReader br = new BufferedReader(new FileReader(configFile));
			//read file setting
			String dataFile = br.readLine();
			TriadTransition tt = TriadTransition.getDynNetFromEdgeFile(dataFile);
			tt.dataName = br.readLine();
			tt.outputDir = br.readLine()+"/"+tt.dataName+"/";
			File outDir = new File(tt.outputDir);
			if(!outDir.exists()) {
				outDir.mkdirs();
				System.out.println("\tcreate: "+ tt.outputDir);
			}
			MotifGraph.initializeTriadCode();
			TriadTransition.initialTmpVar();
			TriadTransition.allEditDistanceFreq = new int[tt.time][7];
			LinkedList<String> lines = new LinkedList<String>();
			String line = null;
			while((line = br.readLine())!= null){
				lines.add(line);
			}
			computeTasks(lines, tt);
			br.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	public static void computeTasks(List<String> tasks, TriadTransition tt){
		TriadTransition.tmpCensus = new int[16];
		MotifGraph pre = new MotifGraph(tt.numNodes, tt.edgeGraph[0], true);
		MotifGraph cur = null;
		// motif transitions
		tt.motifChains = TriadTransition.countTriads(pre, tt.triadFreq[0], 0);;
		int[] subgraph = new int[3];
		int[] iodegrees = new int[3];
		for(int t=1; t<tt.time; t++){
			MathFun.cloneCensus(tt.triadFreq[t], tt.triadFreq[t-1]);
			cur = new MotifGraph(tt.numNodes, tt.edgeGraph[t], true);
			TriadTransition.updateTriadsWNextNetSnapshot(pre, cur, tt.triadFreq[t], t, tt.motifChains, subgraph,
					iodegrees, tt.transitionCnt[t]);
			pre = cur;
		}
		for(String task:tasks){
			//
			if(task.equalsIgnoreCase("triadTransition")){
				System.out.println(task);
				tt.outputAllTransitionMatrices();
			}else if(task.equalsIgnoreCase("triadDistribution")){
				System.out.println("triadDistribution");
				tt.outputTriadCensus();
			}else if(task.equalsIgnoreCase("randomGraphSameDyadDistr")){
				System.out.println("randomGraphSameDyadDistr");
				//suppose a random graph with the same number of dyads, same number of nodes as each snapshot of temporal network, compute the probability distribution of this random network
				tt.outputRandGraphSameDyadTriadDistr();
			}else if(task.equalsIgnoreCase("allTriadTransitionChain")){
				System.out.println("allTriadTransitionChain");
				//for all triads with at least one edge during the process, record their transition between network snapshots
				if(tt.motifChains != null){
					tt.getAllSubGraphChanges();
					tt.outputMotifChainAll3NodeSubGraphs();
				}
			}else if(task.equalsIgnoreCase("randomGraphSameDyad")){
				System.out.println("randomGraphSameDyadFreq");
				//suppose a random graph with the same number of dyads, same number of nodes as each snapshot of temporal network, compute the probability distribution of this random network
				tt.outputRandGraphSameDyadTriadFreq();
			}else if(task.equalsIgnoreCase("randomGraphSameAsymDyadFreq")){
				System.out.println("randomGraphSameAsymDyadFreq");
				//suppose a random graph with the same number of dyads, same number of nodes as each snapshot of temporal network, compute the probability distribution of this random network
				tt.outputRandGraphSameAsymDyadTriadFreq();
			}else if(task.equalsIgnoreCase("splitStaticComponent")){
				System.out.println("splitStaticComponent");
				//suppose a random graph with the same number of dyads, same number of nodes as each snapshot of temporal network, compute the probability distribution of this random network
//				NetworkConnectedComponents ncc = new NetworkConnectedComponents();
//				ncc.getComponentsFromEdgeListsFile(tt.numNodes, tt.time, tt.edgeGraph, true);
//				ncc.nodeIDMapping(true);
//				tt.outputSplitedGraphs(ncc.mapping, ncc.componentSize);
			}
		}
	}
	
	public static void defaultOperation(){
			//javaData = [inDir,'/',dataset{ds},'/',networks{nn},sampleName,granName{tmp}]; EdgeLists.txt
//		String[] dataset={"emailEU", "collegeMsg", "mathoverflow", "askubuntu"};
//			//network name
//		String[][] networkName ={{ "emailEUDept2", "emailEUDept4", "emailEUDept1",  "emailEUDept3","emailEUAll"},
//			{"collegeMsg"},
//			{"mathoverflowA2Q", "mathoverflow-C2Q", "mathoverflowC2A",  "mathoverflowAll"},
//			{"askubuntuAll","askubuntuA2Q","askubuntuC2Q","askubuntuC2A"}};
		String[] dataset={"emailEU"};
		//network name
		String[][] networkName ={{ "emailEUDept2", "emailEUDept4", "emailEUDept1",  "emailEUDept3","emailEUAll"}};
		String inDir = "/nfs/pantanal/scratch1/kuntu/workspace/netMotif/transitions/experimentConfigs/";
		String[] granName = {"Hour","Day","Week","Biweek","FourWeek"};
		String configFile = null;
		for(int i=0; i<dataset.length; i++){
//		for(int i=dataset.length-1; i>0; i--){
			for(String nn: networkName[i]){
				for(int j=1; j< granName.length; j++){
					String gn = granName[j];
					configFile = inDir  +  nn + gn + ".cfg";
					runExperiment(configFile);
				}
			}
		}
	}
	
	public static void runSpeedTestsForDataSets(){
		String[] dataset={"emailEU"};
		//network name
		String[][] networkName ={{ "emailEUDept2", "emailEUDept4", "emailEUDept1",  "emailEUDept3","emailEUAll"}};
		String inDir = "/nfs/pantanal/scratch1/kuntu/workspace/netMotif/transitions/experimentConfigs/";
		String[] granName = {"Hour","Day","Week","Biweek","FourWeek"};
		String configFile = null;
		for(int i=0; i<dataset.length; i++){
//		for(int i=dataset.length-1; i>0; i--){
			for(String nn: networkName[i]){
				for(int j=1; j< granName.length; j++){
					String gn = granName[j];
					configFile = inDir  +  nn + gn + ".cfg";
					speedTest(configFile);
				}
			}
		}
	}
	
	public static void splitDataSets(){
		String[] dataset={"emailEU", "collegeMsg", "mathoverflow", "askubuntu"};
		//network name
		String[][] networkName ={{ "emailEUDept2", "emailEUDept4", "emailEUDept1",  "emailEUDept3","emailEUAll"},
				{"collegeMsg"},
				{"mathoverflowA2Q", "mathoverflowC2Q", "mathoverflowC2A",  "mathoverflowAll"},
				{"askubuntuAll","askubuntuA2Q","askubuntuC2Q","askubuntuC2A"}};
		String[] granName = {"Hour", "Day","Week","Biweek","FourWeek"};
		String configFile = null;
		String inDir = "/nfs/pantanal/scratch1/kuntu/workspace/netMotif/transitions/experimentConfigs/";
		for(int i=0; i<dataset.length; i++){
			for(String nn: networkName[i]){
				for(int j=1; j< granName.length; j++){	// not for hourd
					String gn = granName[j];
					configFile = inDir  +  nn + gn + ".cfg";
					//splitEdgeListNetwork(configFile);
				}
			}
		}
	}
	
//	public static void splitEdgeListNetwork(String configFile){
//		try{
//			BufferedReader br = new BufferedReader(new FileReader(configFile));
//			//read file setting
//			String dataFile = br.readLine();
//			TriadTransition tt = TriadTransition.getDynNetFromEdgeFile(dataFile);
//			tt.dataName = br.readLine();
//			tt.outputDir = br.readLine()+"/"+tt.dataName+"/";	//need to change output dir
//			tt.outputDir = "/nfs/pantanal/scratch1/kuntu/data/dynNet/splitted/"+tt.dataName + "/";
//			File outDir = new File(tt.outputDir);
//			if(!outDir.exists()) {
//				outDir.mkdirs();
//				System.out.println("\tcreate: "+ tt.outputDir);
//			}
//			
//						NetworkConnectedComponents ncc = new NetworkConnectedComponents();
//						ncc.getComponentsFromEdgeListsFile(tt.numNodes, tt.time, tt.edgeGraph, true);
//						ncc.nodeIDMapping(true);
//						tt.outputSplitedGraphs(ncc.mapping, ncc.componentSize);
//			br.close();
//		}catch(Exception e){
//			e.printStackTrace();
//		}
//	}
	
	public static void motifTransitionForSetsOfData(String config){
		try{
			BufferedReader br = new BufferedReader(new FileReader(config));
			String line = null;

			int motifsize = Integer.parseInt(br.readLine());
			String inDir = br.readLine();
			String outDir = br.readLine();
			line = br.readLine();
			String[] dataNames = line.split(";");	//row 3 is data names seperated by ;

			line = br.readLine();
			String[][] networkNames = new String[dataNames.length][];	//row 4 
			String[] tmp = line.split(";");
			for(int i=0; i<tmp.length; i++){
				networkNames[i] = tmp[i].split(",");
			}
			//row 5 is granularity
			line = br.readLine();
			String[] granularities = line.split(",");
			File folder = null;
			
			long[][] motifDists = null;
			br.close();
			String tmpStr = null, tmpfile = null;
			for(String[] networks: networkNames){
				for(String net: networks){
					for(String gran: granularities){
						folder = new File(inDir+"/" + net+gran);
						//System.out.println("input folder: "+ inDir+"/" + net+gran);
						for(File f: folder.listFiles()){
							if(f.isFile()){
								motifDists = getMotifDistributionOfTemporalNetFromFile(motifsize, folder.getAbsolutePath() +"/"+ f.getName());
								tmpStr = outDir+"/" + net+gran;
								File tmpFolder = new File(tmpStr);
								if(!tmpFolder.exists()) tmpFolder.mkdirs();
								System.out.println("Output Dir: " + tmpFolder.getAbsolutePath());
								tmpfile = tmpStr +"/"+f.getName();
								tmpfile = tmpfile.replace(".txt", ".csv");
								BufferedWriter bw = new BufferedWriter(new FileWriter(tmpfile));
								for(long[] snapshot: motifDists){
									for(int i=0; i< snapshot.length-1; i++){
										bw.append(snapshot[i] +"\t");
									}
									bw.append(snapshot[snapshot.length-1]+"\n");
								}
								bw.close();
							}
						}
					}
				}
			}
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	private static void addAllFiles(File folder, List<String> ls){
		if(folder.isFile()) ls.add(folder.getName());
		else{
			String folderName = folder.getName();
			for(File f: folder.listFiles()){
				if(f.isFile()) ls.add(folderName +"/"+ f.getName());
				else if(f.isDirectory()) addAllFiles(f, ls);
			}
		}
		
	}
	
	
	public static void main(String[] args) {
		MotifGraph.initializeTriadCode();
		TriadTransition.initialTmpVar();
		String configFile = null;
		//for(String s: args) System.out.println(s);
		if(args.length>0){
			configFile= args[0];
			if(args[0].equals("speedTest")){
				if(args.length>1){
					speedTest(args[1]);
				}else{
					runSpeedTestsForDataSets();
				}
			}else if(args[0].equals("splitDataSets")){
				splitDataSets();
			}else if(args[0].equals("motifTransitions")){
				if(args.length>1 && args[1] != null){
					motifTransitionForSetsOfData(args[1]);
				}
			}else
				runExperiment(configFile);
		}else {
			System.out.println("config file is not specified, run operation on defaul data sets");
//			defaultOperation();
//			 configFile= "config.txt";
//			runExperiment(configFile);
//			 speedTest(configFile);
			long[][] triadDis = getMotifDistributionOfTemporalNetFromFile(3, "test.txt");
			for(long[] dis: triadDis){
				for(long i: dis) System.out.print(i+"\t");
				System.out.println();
			}
		}
	}

}
