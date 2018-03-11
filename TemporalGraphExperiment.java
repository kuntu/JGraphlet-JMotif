import java.io.*;
import java.util.*;

import motifs.RandomGraphMotif;
import graphs.*;
import randomgraph.RandomGraphJointInOutDegree;
import randomgraph.RandomGraphModel;


public class TemporalGraphExperiment {
	
	public static TemporalGraphWEdgeArray[] getTempGraphsFromCfgFiles(String[] files){
		TemporalGraphWEdgeArray[] tGraphs = new TemporalGraphWEdgeArray[files.length];
		for(int i = 0; i< files.length; i++){
			tGraphs[i] = GraphFactory.getTemproalGraphWEdgeArrayFromFile(files[i]);
		}
		return tGraphs;
	}
	
	public static String[][] getExperimentsFromCfgFile(String cfg){
		String[][] res = null;
		ArrayList<String[]> ls = new ArrayList<String[]>();
		try{
			BufferedReader br = new BufferedReader(new FileReader(cfg));
			String line = null;
			while((line = br.readLine())!= null){
				if(line.startsWith("#")||line.startsWith("//")||line.length() == 0) continue;
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

	public static void main(String[] args) {
		if(args.length < 2){
			System.out.println("need two configuration files as input:\n\t1. configfile for data file paths\n\t2. experiments to execute");
			return ;
		}
		String dataFileCfg = args[0];
		String expCfg = args[1];
		String[] files = GraphIO.getDataFileNamesFromFileList(dataFileCfg);
		//TemporalGraphWEdgeArray[] tGraphs = getTempGraphsFromCfgFiles(files);
		TemporalGraphWEdgeArray tGraph = null;
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
			tGraph = GraphFactory.getTemproalGraphWEdgeArrayFromFile(files[i]);
			for(int j = 2; j<commands.length; j++){
				executeCommand( tGraph, commands[j],outDir+"/"+expName, 
						files[i].substring(files[i].lastIndexOf('/')+1, files[i].lastIndexOf('.')));
			}
		}
	}
	
	public static void executeCommand(TemporalGraphWEdgeArray g, String[] command, String outDir, String outFile0){
		// command[0] represent if operation if for a snapshot, command[1] is the slide number of temproal graph. command[3] is the operation
		if(command.length >= 4&&command[0].equals("snapshot")){
			int snapshotId =Integer.parseInt(command[1]);
			if(snapshotId>=g.time) return;
			String outFile = null;
			File folder = new File(outDir);
			if(!folder.exists()) folder.mkdirs();
			int low =0, high = g.time;
			if(snapshotId >-1){
				low = snapshotId;
				high = snapshotId +1;
			}
			for(snapshotId = low; snapshotId < high; snapshotId++){
				GraphOfEdgeArray graph = g.getSnapshot(snapshotId);
				graph = graph.removeNullNodes();
				outFile = outFile0 + command[0] + snapshotId;
				RandomGraphMotif rgm = null;
				if(command[2].equalsIgnoreCase("joinInOutDegree")) {
					if(command.length>6 && command[command.length-1].equalsIgnoreCase("allowLoopAndMultiEdges")){
							rgm = GraphFactory.getRandomGraphLoopAndMultiEdgeWSameJointIODegree(graph);
							outFile +="loopMultiedge";
					}else rgm = GraphFactory.getRandomGraphWSameJointIODegree(graph);
				}
				else if(command[2].equalsIgnoreCase("reciprocalInOutDegree")) rgm = GraphFactory.getRandomGraphReciprocalAndInOutDegreeFromGraphOfEdgeArray(graph);
				else if(command[2].equalsIgnoreCase("numNodeEdge")) rgm = GraphFactory.getRandomGraphWNumNodeEdge(graph);
				else if(command[2].equalsIgnoreCase("MANPairModel")) rgm = GraphFactory.getRandomGraphW_MAN_Pair(graph);
				outFile = outFile0 + command[0] + snapshotId + command[2];
				if(command.length>= 5 && command[3].equalsIgnoreCase("generateRandomGraph")){
					RandomGraphModel rg = (RandomGraphModel) rgm;
					int repeat = Integer.parseInt(command[4]);
					System.out.printf("\n\t[Operation]: generate %d random graph %s \n", repeat, command[2]);
					int[][] outputM = null;
					for(int i =0; i< repeat; i++){
						outputM = rg.generateRandomGraph().edges;
						GraphIO.outputMatrix(outDir+"/"+outFile+"_RandomGraph"+ i+".txt", outputM);
					}
				}else if(command[3].equalsIgnoreCase("motifCensus")){
					System.out.println("\n\t");
					//count motif census
					int mSize = Integer.parseInt(command[1]);
					if(mSize == 3){
						double time = System.currentTimeMillis() ;
						long[] freq = graph.getMotifFreq(mSize);
						time  = System.currentTimeMillis() - time;
						double[][] outputM = new double[1][16];
						for(int i=0; i< 16; i++) outputM[0][i] = freq[i];
						System.out.println("[output]: "+outDir+"/"+outFile + "TriadCencus.txt");
						GraphIO.outputMatrix(outDir+"/"+outFile + "TriadCencus.txt", outputM);
						outputM[0] = new double[]{time};
						if(command.length > 4 && command[4].equalsIgnoreCase("computeTime"))
							GraphIO.outputMatrix(outDir+"/"+outFile + "TriadCencus_time.txt", outputM);
					}
				}else if(command[3].equalsIgnoreCase("ExpectedTriadFreq")){
					double[] expectMotifFreq = rgm.getMotifFreq(Integer.parseInt(command[4]));
					double[][] outputM = new double[1][];
					outputM[0] = expectMotifFreq;
					GraphIO.outputMatrix(outDir+"/"+outFile + command[3] + ".txt", outputM);
				}else if(command[3].equalsIgnoreCase("SampleGraphMotifFreq")){
					//output simulation result
					int repeat = Integer.parseInt(command[5]);
					long startime = System.currentTimeMillis();
					double[][] outputM = rgm.getMotifFreqFromSampledGraphs(Integer.parseInt(command[4]), repeat);
					double time = System.currentTimeMillis() - startime;
					GraphIO.outputMatrix(outDir+"/"+outFile+"sampleGraphMotifFreq.txt", outputM);
					outputM = new double[][]{{(double)repeat, time}};
					GraphIO.outputMatrix(outDir+"/"+outFile+"sampleGraphMotifFreq_time.txt", outputM);
				}else if(command[3].equalsIgnoreCase("ExpectedTriadFreqOptTime")){
					int repeat = (int) Math.max(1, Integer.parseInt(command[4]));	//num of repeating experiment
					int opt = 0;
					if(command.length>5) opt = Integer.parseInt(command[5]);	//choose method to compute connection probability
					double[][] outputM = null;
					double[] runTimes = new double[repeat];
					long startime = 0;
					//compute expected triad freq
					outputM = new double[1][];
					for(int i=0; i<repeat; i++){
						startime = System.currentTimeMillis();
						outputM[0] = rgm.getMotifFreq(Integer.parseInt(command[3]));
						runTimes[i] = (double) (System.currentTimeMillis() - startime);
					}
					GraphIO.outputMatrix(outDir+"/"+outFile + command[3] +"opt"+ opt+ ".txt", outputM);
					outputM = new double[1][];
					outputM[0] = runTimes;
					GraphIO.outputMatrix(outDir +"/" + outFile + command[3] +"opt"+ opt+ "_time.txt", outputM);
				}else if(command[3].equalsIgnoreCase("getRandomGraphInfo")){
					String last = "Info.txt";
					//output snapshot info:
					RandomGraphModel rg = (RandomGraphModel) rgm;
					double[][] outputM = rg.getGraphInfo();
					GraphIO.outputMatrix(outDir+"/"+outFile + last, outputM);
				}else if(command[3].equalsIgnoreCase("outputConnectComponents") && command.length > 4){
					try{
						File f = new File(command[4]);	//directory for outputed connected components
						if(!f.exists()) f.mkdirs();
						int compSize = Integer.parseInt(command[4]);
						int[][][] tEdges = GraphPropertiesToolBox.obtainConnectComponentEdges(graph.edges, compSize);
						for(int i = 0; i<tEdges.length; i++){
							BufferedWriter bw = new BufferedWriter(new FileWriter(f+"/"+ outFile+"Comp"+i+".txt"));
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
					
				}
			}
		}
	}

}
