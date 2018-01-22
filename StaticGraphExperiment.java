import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
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
	
	public static void executeCommand(GraphOfEdgeArray graph, String[] command, String outDir, String outFile){
		if(command[0].equalsIgnoreCase("static")){
			File folder = new File(outDir);
			if(!folder.exists()) folder.mkdirs();
			if(command[1].equalsIgnoreCase("ExpectedTriadFreqOptTime") && command.length >= 5){
				System.out.println("\n\t[Operation]:Analogical computation random graph and count motif frequency");
				int motifSize = Integer.parseInt(command[2]);
				int opt = Integer.parseInt(command[3]);
				int repeat = Integer.parseInt(command[4]);
				double[] runTimes = new double[repeat];
				double[][] outputM = null;
				long startime = 0;
				RandomGraphJointInOutDegree rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
				//output snapshot info:
//				outputM = rgjiod.getGraphInfo();
//				GraphIO.outputMatrix(outDir+"/"+outFile + "Info.txt", outputM);
				//compute 
				outputM = new double[1][];
				for(int i=0; i<repeat; i++){
					startime = System.currentTimeMillis();
					outputM[0] = rgjiod.getMotifFreq(motifSize, opt);
					runTimes[i] = (double) (System.currentTimeMillis() - startime);
				}
				GraphIO.outputMatrix(outDir+"/"+outFile + command[1] +"opt"+ opt+ ".txt", outputM);
				outputM[0] = runTimes;
				GraphIO.outputMatrix(outDir +"/" + outFile + command[1] +"opt"+ opt+ "_time.txt", outputM);
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
				GraphIO.outputMatrix(outDir+"/"+outFile+"sampleGraphMotifFreq.txt", outputM);
				outputM = new double[][]{{(double)repeat, time}};
				GraphIO.outputMatrix(outDir+"/"+outFile+"sampleGraphMotifFreq_time.txt", outputM);
			}else if(command[1].equalsIgnoreCase("getRandomGraphInfo") ){
				System.out.println("\n\t[Operation]:Output Graph Info");
				double[][] outputM = null;
				RandomGraphJointInOutDegree rgjiod = GraphFactory.getRandomGraphWSameJointIODegree(graph);
				//output snapshot info:
				outputM = rgjiod.getGraphInfo();
				GraphIO.outputMatrix(outDir+"/"+outFile + "Info.txt", outputM);
			}else if(command[1].equalsIgnoreCase("motifCensus") && command.length >2){
				System.out.println("\n\t");
				//count motif census
				int mSize = Integer.parseInt(command[2]);
				if(mSize == 3){
					double time = System.currentTimeMillis() ;
					long[] freq = graph.getMotifFreq(mSize);
					time  = System.currentTimeMillis() - time;
					double[][] outputM = new double[1][16];
					for(int i=0; i< 16; i++) outputM[0][i] = freq[i];
					System.out.println("[output]: "+outDir+"/"+outFile + "TriadCencus.txt");
					GraphIO.outputMatrix(outDir+"/"+outFile + "TriadCencus.txt", outputM);
					outputM[0] = new double[]{time};
					GraphIO.outputMatrix(outDir+"/"+outFile + "TriadCencus_time.txt", outputM);
				}
			}
		}
	}
	
	public static void main(String[] args) {
		if(args.length < 2){
			System.out.println("need two configuration files as input:\n\t1. configfile for data file paths\n\t2. experiments to execute");
			return ;
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
			for(int j = 2; j<commands.length; j++){
				executeCommand( graph, commands[j],outDir+"/"+expName, 
						files[i].substring(files[i].lastIndexOf('/')+1, files[i].lastIndexOf('.')));
			}
		}
	}

}
