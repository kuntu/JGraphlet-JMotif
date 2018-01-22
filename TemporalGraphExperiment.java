import java.io.*;
import java.util.*;

import motifs.RandomGraphMotif;
import graphs.*;
import randomgraph.RandomGraphJointInOutDegree;


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
				if(line.startsWith("//")||line.length() == 0) continue;
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
		// command[0] represent if operation if for a snapshot, command[1] is the slide number of temproal graph. command[2] is the operation
		if(command[0].equals("snapshot")){
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
				outFile = outFile0 + command[0] + snapshotId;
				if(command[2].equalsIgnoreCase("randomGraphIOSeqCompareMotifFreqAlg")){
					RandomGraphJointInOutDegree rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
					long startime = System.currentTimeMillis();
					double[] expectMotifFreq = rgjiod.getMotifFreq(Integer.parseInt(command[3]));
					System.out.println(outFile+"expect: "+(System.currentTimeMillis() - startime));
					double[][] outputM = new double[1][];
					outputM[0] = expectMotifFreq;
					GraphIO.outputMatrix(outDir+"/"+outFile+"ExpectMotifFreq.txt", outputM);
					startime = System.currentTimeMillis();
					double[][] sampleGraphMotifFreq = rgjiod.getMotifFreqFromSampledGraphs(Integer.parseInt(command[3]), Integer.parseInt(command[4]));
					System.out.println(outFile+"sample: "+(System.currentTimeMillis() - startime));
					outputM = sampleGraphMotifFreq;
					GraphIO.outputMatrix(outDir+"/"+outFile+"sampleGraphMotifFreq.txt", outputM);
					//output snapshot info:
					outputM = rgjiod.getGraphInfo();
					GraphIO.outputMatrix(outDir+"/"+outFile + "Info.txt", outputM);
				}else if(command[2].equalsIgnoreCase("ExpectedTriadFreq")){
					RandomGraphJointInOutDegree rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
					double[] expectMotifFreq = rgjiod.getMotifFreq(Integer.parseInt(command[3]));
					double[][] outputM = new double[1][];
					outputM[0] = expectMotifFreq;
					GraphIO.outputMatrix(outDir+"/"+outFile + command[2] + ".txt", outputM);
					//output snapshot info:
					outputM = rgjiod.getGraphInfo();
					GraphIO.outputMatrix(outDir+"/"+outFile + "Info.txt", outputM);
				}else if(command[2].equalsIgnoreCase("SampleGraphMotifFreq")){
					//output simulation result
					if(command.length>=6 && command[5].equalsIgnoreCase("removeNullNode")){
						graph = graph.removeNullNodes();
					}
					RandomGraphMotif rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
					int repeat = Integer.parseInt(command[4]);
					long startime = System.currentTimeMillis();
					double[][] outputM = rgjiod.getMotifFreqFromSampledGraphs(Integer.parseInt(command[3]), repeat);
					double time = System.currentTimeMillis() - startime;
					GraphIO.outputMatrix(outDir+"/"+outFile+"sampleGraphMotifFreq.txt", outputM);
					outputM = new double[][]{{(double)repeat, time}};
					GraphIO.outputMatrix(outDir+"/"+outFile+"sampleGraphMotifFreq_time.txt", outputM);
				}else if(command[2].equalsIgnoreCase("ExpectedTriadFreqOptTime")){
					if(command.length<6) return;
					if(command.length>=7 && command[6].equalsIgnoreCase("removeNullNode")){
						graph = graph.removeNullNodes();
					}
					int repeat = (int) Math.max(1, Integer.parseInt(command[5]));	//num of repeating experiment
					int opt = Integer.parseInt(command[4]);	//choose method to compute connection probability
					RandomGraphJointInOutDegree rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
					double[][] outputM = null;
					double[] runTimes = new double[repeat];
					long startime = 0;
					//compute expected triad freq
					outputM = new double[1][];
					for(int i=0; i<repeat; i++){
						startime = System.currentTimeMillis();
						outputM[0] = rgjiod.getMotifFreq(Integer.parseInt(command[3]), opt);
						runTimes[i] = (double) (System.currentTimeMillis() - startime);
					}
					GraphIO.outputMatrix(outDir+"/"+outFile + command[2] +"opt"+ opt+ ".txt", outputM);
					outputM = new double[1][];
					outputM[0] = runTimes;
					GraphIO.outputMatrix(outDir +"/" + outFile + command[2] +"opt"+ opt+ "_time.txt", outputM);
				}else if(command[2].equalsIgnoreCase("getRandomGraphInfo")){
					if(command.length>=4 && command[3].equalsIgnoreCase("removeNullNode")){
						System.out.println("\n reduce network size");
						graph = graph.removeNullNodes();
					}
					RandomGraphJointInOutDegree rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
					//output snapshot info:
					double[][] outputM = rgjiod.getGraphInfo();
					GraphIO.outputMatrix(outDir+"/"+outFile + "Info.txt", outputM);
				}
			}
		}
	}

}
