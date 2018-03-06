import graphs.GraphFactory;
import graphs.GraphIO;
import graphs.GraphOfEdgeArray;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

import motifs.RandomGraphMotif;
import randomgraph.RandomGraphJointInOutDegree;


public class DataPreprocessing {

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
	
	public static void executeCommand(String[] command, String inFile, String outDir, String outFile){
		if(command[0].equalsIgnoreCase("Preprocessing")){
			File folder = new File(outDir);
			if(!folder.exists()) folder.mkdirs();
			if(command[1].equalsIgnoreCase("convertMatrixToEdgeList")){
				System.out.println("\n\t[Operation]:conver matrix to edge list");
				int[][] m = GraphIO.convertAdjacentMatrixToEdgeList(GraphIO.getMatrixFromFile(inFile));
				outFile = outDir + "/"  + inFile.substring(inFile.lastIndexOf('/')+1, inFile.lastIndexOf("_adj")) + "EdgeList.txt";
				GraphIO.outputMatrix(outFile, m);
			}else if(command[1].equalsIgnoreCase("ID2Integer")){
				boolean isEgo = false;
				if(command.length > 2 && command[2].equalsIgnoreCase("ego")) isEgo = true;
				GraphIO.convertStrIDToIntInFile(inFile, isEgo);
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
		int[][] intMtr = null;
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
			//intMtr = GraphIO.getMatrixFromFile(files[i]);
			for(int j = 2; j<commands.length; j++){
				executeCommand( commands[j], files[i], outDir+"/"+expName, null);
						//files[i].substring(files[i].lastIndexOf('/')+1, files[i].lastIndexOf('.')));
			}
		}
	}

}
