import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

import motifs.RandomGraphMotif;
import graphs.*;
import randomgraph.*;

public class StaticGraphExperiment {
	
	public static String[][] getExperimentsFromCfgFile(String cfg){
		String[][] res = null;
		ArrayList<String[]> ls = new ArrayList<String[]>();
		try{
			BufferedReader br = new BufferedReader(new FileReader(cfg));
			String line = null;
			while((line = br.readLine())!= null){
				if(line.startsWith("//")|| line.startsWith("#")|| line.length() == 0) continue;
				ls.add(line.split("\\s?,?\\s+"));
			}
			res = new String[ls.size()][];
			int idx = 0;
			for(String[] l: ls) res[idx++] = l;
			br.close();
		}catch(Exception e){
			e.printStackTrace();
			res = new String[0][1];
		}
		return res;
	}
	
	/**
	 * 
	 * @param graph
	 * @param command
	 * @param outDir	output folder
	 * @param outFile dataName
	 */
	public static void executeCommand(GraphOfEdgeArray graph, String[] command, String outDir, String outFile, int fIdx){
		//for(int i=0; i< command.length ; i++) System.out.print(command[i] + " ");
		File folder = new File(outDir);
		if(!folder.exists()) folder.mkdirs();
		String fType = ".txt";
		if(fIdx>0) fType += fIdx;
		if(command[0].equals("conversion")){	// converstion
			if(command.length> 1 && command[1].equalsIgnoreCase("nodeList")){
				int[][] output = new int[1][];
				output[0] = GraphIO.edgesToNodeList(graph.edges);
				System.out.println("\toutput:"+ outDir +"/"+outFile+"NodeList"+fType);
				GraphIO.outputMatrix(outDir +"/"+outFile+"NodeList"+fType, output);
			}
		}else if(command.length >= 3&&command[1].equalsIgnoreCase("generateRandomGraph")){
			RandomGraphModel rg = null;	
			if(command[0].equalsIgnoreCase("MANPairModel")){
				rg = GraphFactory.getRandomGraphW_MAN_Pair(graph);
			}else if(command[0].equalsIgnoreCase("numNodeEdge")){
				//random graph with same number of nodes and edges
				rg = GraphFactory.getRandomGraphWNumNodeEdge(graph);
			}else if(command[0].equalsIgnoreCase("reciprocalInOutDegree")){
				rg = GraphFactory.getRandomGraphReciprocalAndInOutDegreeFromGraphOfEdgeArray(graph);
			}else if(command[0].equalsIgnoreCase("joinInOutDegree")){
				rg = GraphFactory.getRandomGraphWSameJointIODegree(graph);
			}else if(command[0].equalsIgnoreCase("cfgModel")){
				rg = GraphFactory.getRandomGraphLoopAndMultiEdgeWSameJointIODegree(graph);
			}
			int repeat = Integer.parseInt(command[2]);
			System.out.printf("\n\t[Operation]: generate %d random graph %s \n", repeat, command[0]);
			//GraphIO.outputMatrix(outDir+"/"+outFile+"RandomGraph_Info"+fType, rg.getGraphInfo());
			int[][] outputM = null;
			for(int i =0; i< repeat; i++){
				outputM = rg.generateRandomGraph().edges;
				System.out.println("\t***"+outDir+"/"+outFile+"_RandomGraph"+ i+""+fType);
				GraphIO.outputMatrix(outDir+"/"+outFile+"_"+command[0]+"_RandomGraph"+ i+""+fType, outputM);
			}
		}else if(command[0].equalsIgnoreCase("motifCensus") && command.length >1){
			System.out.println("\n\t");
			//count motif census
			int mSize = Integer.parseInt(command[1]);
			if(mSize == 3){
				double time = System.currentTimeMillis() ;
				long[] freq = graph.getMotifFreq(mSize);
				time  = System.currentTimeMillis() - time;
				double[][] outputM = new double[1][16];
				for(int i=0; i< 16; i++) outputM[0][i] = freq[i];
				System.out.println("[output]: "+outDir+"/"+outFile + "TriadCencus"+fType);
				GraphIO.outputMatrix(outDir+"/"+outFile + "TriadCencus"+fType, outputM);
				outputM[0] = new double[]{time};
				GraphIO.outputMatrix(outDir+"/"+outFile + "TriadCencus_time"+fType, outputM);
			}else if(mSize==-3){
				long[] freq = graph.getMotifFreq(mSize);
				double[][] outputM = new double[1][freq.length];
				for(int i=0; i< freq.length; i++) outputM[0][i] = freq[i];
				System.out.println("[output]: "+outDir+"/"+outFile + "TriadCencus"+fType);
				GraphIO.outputMatrix(outDir+"/"+outFile + "TriadCencus"+fType, outputM);
			}
		}else if(command[0].equalsIgnoreCase("motifNodeCensus") && command.length >1){
			System.out.println("\n\t");
			//count motif census
			int mSize = Integer.parseInt(command[1]);
			if(mSize == 3){
				double time = System.currentTimeMillis() ;
				long[][] freq = graph.getNodeMotifFreq(mSize);
				time  = System.currentTimeMillis() - time;
				double[][] outputM = new double[freq.length][17];
				for(int i=0; i< outputM.length; i++) {
					for(int j = 0; j < outputM[0].length; j++){
						outputM[i][j] = freq[i][j];
					}
				}
				System.out.println("[output]: "+outDir+"/"+outFile + "Node_TriadCencus"+fType);
				GraphIO.outputMatrix(outDir+"/"+outFile + "Node_TriadCencus"+fType, outputM);
				double[][] timeM = new double[1][1];
				timeM[0][0] = time;
				GraphIO.outputMatrix(outDir+"/"+outFile + "Node_TriadCencus_time"+fType, timeM);
			}
		}else if(command[0].equalsIgnoreCase("joinInOutDegree")){
			
			if(command[1].equalsIgnoreCase("ExpectedTriadFreqOptTime") && command.length >= 4){
				System.out.println("\n\t[Operation]:Analogical computation random graph and count motif frequency");
				int motifSize = Integer.parseInt(command[2]);
				int opt = 0;
				if(command.length >=5) opt = Integer.parseInt(command[4]);
				int repeat = Integer.parseInt(command[3]);
				double[] runTimes = new double[repeat];
				double[][] outputM = null;
				long startime = 0;
				RandomGraphJointInOutDegree rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
				//output snapshot info:
//				outputM = rgjiod.getGraphInfo();
//				GraphIO.outputMatrix(outDir+"/"+outFile + "Info"+fType, outputM);
				//compute 
				outputM = new double[1][];
				for(int i=0; i<repeat; i++){
					startime = System.currentTimeMillis();
					outputM[0] = rgjiod.getMotifFreq(motifSize, opt);
					runTimes[i] = (double) (System.currentTimeMillis() - startime);
				}
				GraphIO.outputMatrix(outDir+"/"+outFile + command[1] +"opt"+ opt+ ""+fType, outputM);
				outputM[0] = runTimes;
				GraphIO.outputMatrix(outDir +"/" + outFile + command[1] +"opt"+ opt+ "_time"+fType, outputM);
			}else if(command[1].equalsIgnoreCase("sampleGraphMotifFreq")&& command.length >= 4){
				System.out.println("\n\t[Operation]:Simulation to generate multiple random graph and count motif frequency");
				RandomGraphMotif rgjiod = null;
//				rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
				if(command.length>=5 && command[4].equals("allowLoopAndMultiEdges")) {
					rgjiod = GraphFactory.getRandomGraphLoopAndMultiEdgeWSameJointIODegree(graph);
					outFile +="loopMultiedge";
				}else rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
				int repeat = Integer.parseInt(command[3]);
				int motifSize = Integer.parseInt(command[2]);
				long startime = System.currentTimeMillis();
				double[][] outputM = rgjiod.getMotifFreqFromSampledGraphs(motifSize, repeat);
				double time = System.currentTimeMillis() - startime;
				GraphIO.outputMatrix(outDir+"/"+outFile+"sampleGraphMotifFreq"+fType, outputM);
				outputM = new double[][]{{(double)repeat, time}};
				GraphIO.outputMatrix(outDir+"/"+outFile+"sampleGraphMotifFreq_time"+fType, outputM);
			}else if(command[1].equalsIgnoreCase("getRandomGraphInfo") ){
				System.out.println("\n\t[Operation]:Output Graph Info");
				String last  = "Info"+fType;
				if(command.length>=3 && command[2].equalsIgnoreCase("removeNullNode")){
					System.out.println("\n reduce network size");
					graph = graph.removeNullNodes();
					last = "NoNullNode" + last;
				}
				double[][] outputM = null;
				RandomGraphJointInOutDegree rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
				//output snapshot info:
				outputM = rgjiod.getGraphInfo();
				GraphIO.outputMatrix(outDir+"/"+outFile + last, outputM);
			}else if(command[1].equalsIgnoreCase("outputConnectComponents") && command.length > 3){
				System.out.println("\t\t[Operation]: obtain Larges Components");
				try{
					File f = new File(command[2]);
					if(!f.exists()) f.mkdirs();
					int compSize = Integer.parseInt(command[3]);
					int[][][] tEdges = GraphPropertiesToolBox.obtainConnectComponentEdges(graph.edges, compSize);
					for(int i = 0; i<tEdges.length; i++){
						BufferedWriter bw = new BufferedWriter(new FileWriter(f+"/"+ outFile+"Comp"+i+""+fType));
						System.out.println("output file: "+ f+"/"+ outFile+"Comp"+i+""+fType);
						bw.write(tEdges[i][0][0] +"\t" + tEdges[i][0][1]);
						for(int j = 1; j< tEdges[i].length; j++){
							bw.newLine();
							bw.write(tEdges[i][j][0] +"\t" + tEdges[i][j][1]);
						}
						bw.close();
					}
				}catch(Exception e){
					e.printStackTrace();
				}
				
			}else if(command[1].equalsIgnoreCase("ExpectedTriadFreq") && command.length >= 3){
				System.out.println("\n\t[Operation]:Analogical computation random graph with joint in/out degree and count motif frequency");
				int motifSize = Integer.parseInt(command[2]);
				int opt = 0;
				int repeat = Integer.parseInt(command[3]);
				double[] runTimes = new double[repeat];
				double[][] outputM = null;
				long startime = 0;
				RandomGraphJointInOutDegree rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
				outputM = new double[1][];
				for(int i=0; i<repeat; i++){
					startime = System.currentTimeMillis();
					outputM[0] = rgjiod.getMotifFreq(motifSize, opt);
					runTimes[i] = (double) (System.currentTimeMillis() - startime);
				}
				GraphIO.outputMatrix(outDir+"/"+outFile + command[1] +"opt"+ opt+ ""+fType, outputM);
				outputM[0] = runTimes;
				GraphIO.outputMatrix(outDir +"/" + outFile + command[1] +"opt"+ opt+ "_time"+fType, outputM);
			}
		}else if(command[0].equalsIgnoreCase("reciprocalInOutDegree")){
			if(command[1].equalsIgnoreCase("sampleGraphMotifFreq")&& command.length >= 4){
				System.out.println("\n\t[Operation]:Simulation to generate multiple random graph and count motif frequency");
				RandomGraphMotif rgriod = null;
				boolean loop = false;
				if(command.length>=5 && command[4].equals("allowLoopAndMultiEdges")) {
					rgriod = GraphFactory.getRandomGraphReciprocalAndInOutDegreeFromGraphOfEdgeArray(graph);
					loop =true;
				}else rgriod = GraphFactory.getRandomGraphReciprocalAndInOutDegreeFromGraphOfEdgeArray(graph);
				int repeat = Integer.parseInt(command[3]);
				int motifSize = Integer.parseInt(command[2]);
				long startime = System.currentTimeMillis();
				double[][] outputM = rgriod.getMotifFreqFromSampledGraphs(motifSize, repeat);
				double time = System.currentTimeMillis() - startime;
				String tmp = "sampleReciprocalInOutDegreeGraphMotifFreq";
				if(loop) tmp += "loopMultiedges";
				GraphIO.outputMatrix(outDir+"/"+outFile+tmp+""+fType, outputM);
				outputM = new double[][]{{(double)repeat, time}};
				GraphIO.outputMatrix(outDir+"/"+outFile+tmp+""+fType, outputM);
			}
		}else if(command[0].equalsIgnoreCase("numNodeEdge")){
			if(command[1].equalsIgnoreCase("sampleGraphMotifFreq")&& command.length >= 4){
				System.out.println("\n\t[Operation]:Simulation to generate multiple random graph with same numNode and numEdge and count motif frequency");
				RandomGraphMotif rgne = null;	//random graph with same number of nodes and edges
				rgne = GraphFactory.getRandomGraphWNumNodeEdge(graph);
				int repeat = Integer.parseInt(command[3]);
				int motifSize = Integer.parseInt(command[2]);
				long startime = System.currentTimeMillis();
				double[][] outputM = rgne.getMotifFreqFromSampledGraphs(motifSize, repeat);
				double time = System.currentTimeMillis() - startime;
				GraphIO.outputMatrix(outDir+"/"+outFile+"sampleWNumNodeEdgeGraphMotifFreq"+fType, outputM);
				outputM = new double[][]{{(double)repeat, time}};
				GraphIO.outputMatrix(outDir+"/"+outFile+"sampleWNumNodeEdgeGraphMotifFreq_time"+fType, outputM);
			}else if(command[1].equalsIgnoreCase("ExpectedTriadFreq") && command.length >= 4){
				System.out.println("\n\t[Operation]:Analogical computation random graph w same # nodes and # edges and count motif frequency");
				int motifSize = Integer.parseInt(command[2]);
				int repeat = Integer.parseInt(command[3]);
				double[] runTimes = new double[repeat];
				double[][] outputM = null;
				long startime = 0;
				RandomGraphMotif rgne = GraphFactory.getRandomGraphWNumNodeEdge(graph);
				//compute 
				outputM = new double[1][];
				for(int i=0; i<repeat; i++){
					startime = System.currentTimeMillis();
					outputM[0] = rgne.getMotifFreq(motifSize);
					runTimes[i] = (double) (System.currentTimeMillis() - startime);
				}
				GraphIO.outputMatrix(outDir+"/"+outFile + command[1] +"WNumNodeEdge"+ ""+fType, outputM);
				outputM[0] = runTimes;
				GraphIO.outputMatrix(outDir +"/" + outFile + command[1] +"WNumNodeEdg"+ "_time"+fType, outputM);
			}
		}else if(command[0].equalsIgnoreCase("MANPairModel")){
			if(command[1].equalsIgnoreCase("sampleGraphMotifFreq")&& command.length >= 4){
				System.out.println("\n\t[Operation]:Simulation to generate multiple random graph with same MAN pair and count motif frequency");
				RandomGraphMotif rgMAN = null;	//random graph with same number of nodes and edges
				rgMAN = GraphFactory.getRandomGraphW_MAN_Pair(graph);
				int repeat = Integer.parseInt(command[3]);
				int motifSize = Integer.parseInt(command[2]);
				long startime = System.currentTimeMillis();
				double[][] outputM = rgMAN.getMotifFreqFromSampledGraphs(motifSize, repeat);
				double time = System.currentTimeMillis() - startime;
				GraphIO.outputMatrix(outDir+"/"+outFile+"sample_MANPair_MotifFreq"+fType, outputM);
				outputM = new double[][]{{(double)repeat, time}};
				GraphIO.outputMatrix(outDir+"/"+outFile+"sample_MANPair_MotifFreq_time"+fType, outputM);
			}else if(command[1].equalsIgnoreCase("ExpectedTriadFreq") && command.length >= 4){
				System.out.println("\n\t[Operation]:Analogical computation random graph w Specified MAN node-Pair and count motif frequency");
				int motifSize = Integer.parseInt(command[2]);
				int repeat = Integer.parseInt(command[3]);
				double[] runTimes = new double[repeat];
				double[][] outputM = null;
				long startime = 0;
				RandomGraphMotif rgMAN = GraphFactory.getRandomGraphW_MAN_Pair(graph);
				//compute 
				outputM = new double[1][];
				for(int i=0; i<repeat; i++){
					startime = System.currentTimeMillis();
					outputM[0] = rgMAN.getMotifFreq(motifSize);
					runTimes[i] = (double) (System.currentTimeMillis() - startime);
				}
				GraphIO.outputMatrix(outDir+"/"+outFile + command[1] +"WMANPair"+ ""+fType, outputM);
				outputM[0] = runTimes;
				GraphIO.outputMatrix(outDir +"/" + outFile + command[1] +"WMANPair"+ "_time"+fType, outputM);
			}else if(command[1].equalsIgnoreCase("getRandomGraphInfo")){
				System.out.println("\n\t[Operation]:Output Graph Info");
				String last  = "_MANPair_Info"+fType;
				if(command.length>=3 && command[2].equalsIgnoreCase("removeNullNode")){
					System.out.println("\n reduce network size");
					graph = graph.removeNullNodes();
					last = "NoNullNode" + last;
				}
				double[][] outputM = null;
				RandomGraphNumMANPair rgMAN = GraphFactory.getRandomGraphW_MAN_Pair(graph);
				//output snapshot info:
				outputM = rgMAN.getGraphInfo();
				GraphIO.outputMatrix(outDir+"/"+outFile + last, outputM);
			}
		}
	}
	
	public static void main(String[] args) {
		if(args.length < 2){
			System.out.println("need two configuration files as input:\n\t1. configfile for data file paths\n\t2. experiments to execute");
			return ;
		}
		int fidx = -1;
		if(args.length>2 && args[2].equalsIgnoreCase("-i")){
			fidx = 1;
		}
		String dataFileCfg = args[0];
		String expCfg = args[1];
		String[] files = GraphIO.getDataFileNamesFromFileList(dataFileCfg);
		GraphOfEdgeArray graph = null;
		//GraphOfEdgeArray[] tGraphs = getGraphsFromCfgFiles(files);
		String[][] commands = getExperimentsFromCfgFile(expCfg);
		if(commands[0].length < 2 || commands[1].length < 2 || !commands[0][0].equalsIgnoreCase("expName")
				|| !commands[1][0].equalsIgnoreCase("outDir")){
			System.out.println("experiment configuration file format error:"
					+ "\n\tThe first line should be in the format: expName experiment_name"
					+ "\n\tThe second line should be in the format: outDir output_folder_for_result");
			return ;
		}
		String expName = commands[0][1]; // first line is experiment Name
		String outDir = commands[1][1];	// second line is output folder
		for(int i=0; i<files.length; i++){
			System.out.printf("Process file[%d]: %s\n", i, files[i]);
			graph = GraphFactory.makeEdgeGraphFromFileOfEdgeList(files[i]);
			//graph = graph.removeNullNodes();
			for(int j = 2; j<commands.length; j++){
				executeCommand( graph, commands[j],outDir+"/"+expName, 
						files[i].substring(files[i].lastIndexOf('/')+1, files[i].lastIndexOf('.')), (i+1) * fidx);
			}
		}
	}

}
