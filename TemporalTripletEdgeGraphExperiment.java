import java.io.File;

import randomgraph.TemporalRandomGraph;
import randomgraph.TemporalTimeStampRandomGraph;
import graphs.GraphFactory;
import graphs.TemporalTripletEdgeGraph;


public class TemporalTripletEdgeGraphExperiment {
	public static void main(String[] args){
		if(args.length<2){
			System.out.println("input the filename that store data_file_list and the filename that store operations");
			return;
		}
		String[][][] tmp = ExperimentPipeline.getExperimentSettingsFromArgs(args);
		if(tmp == null) return;
		String[] datafiles = tmp[0][0];
		String[][] commands = tmp[1];
		//
		String expName = commands[0][1];
		String outDir = commands[1][1];
		File folder = new File(outDir+"/"+expName);
		if(!folder.exists()) folder.mkdirs();
		TemporalTripletEdgeGraph tg = null;
		for(int i = 0; i< datafiles.length; ++i){
			tg = GraphFactory.makeTemporalTripletGraphFromFile(datafiles[i]);
			for(int j = 2; j<commands.length; ++j){
				if(commands[j].length==0) continue;
				executeCommands(commands[j], tg.tripletEdges, expName, outDir, datafiles[i].substring(datafiles[i].lastIndexOf('/')+1, datafiles[i].lastIndexOf('.')));
			}
		}
	}
	
	private static void executeCommands(String[] cmd, long[][] edges, String expName, String outDir, String dataFilePrefix){
		if(cmd.length > 2 && cmd[0].equalsIgnoreCase("generateRandomGraph")){
			TemporalRandomGraph tg = null;
			if(cmd[1].equalsIgnoreCase("timestamp")){	//generate random graph by shuffle the time stamps of the tempoal triplet edges
				tg = new TemporalTimeStampRandomGraph(edges);
			}else{
				
			}
			if(tg == null){
				System.out.println("no support for random graphs type: " + cmd[1]);
				return;
			}
			int num = Integer.parseInt(cmd[2]);
			
			tg.saveRandomGraphs(outDir+"/"+expName+"/"+dataFilePrefix, num);
		}else{
			
		}
	}
}
