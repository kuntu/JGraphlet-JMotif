import graphs.GraphIO;


public class ExperimentPipeline {
	public static String[][][] getExperimentSettingsFromArgs(String[] args){
		String[][][] res = new String[2][][];
		String[] datafiles = GraphIO.getDataFileNamesFromFileList(args[0]);
		if(datafiles == null) return null;
		else{
			res[0] = new String[][]{datafiles};
		}
		String[][] commands = GraphIO.getExperimentCommandsFromCfgFile(args[1]);
		if(commands == null||commands.length==0|| commands[0].length < 2 || commands[1].length < 2 || !commands[0][0].equalsIgnoreCase("expName")
				|| !commands[1][0].equalsIgnoreCase("outDir")){
			System.out.println("experiment configuration file format error:"
					+ "\n\tThe first line should be in the format: expName experiment_name"
					+ "\n\tThe second line should be in the format: outDir output_folder_for_result");
			return null;
		}else res[1] = commands;
		return res;
	}
}
