package graphs;
import java.io.*;
import java.util.*;


public class GraphIO {
	
	/**
	 * obtain file names of data files from a configuration file.
	 * @param dataCfgFile
	 * @return
	 */
	public static String[] getDataFileNamesFromConfig(String dataCfgFile){
		String[] files = null;
		try {
			BufferedReader br = new BufferedReader(new FileReader(dataCfgFile));
			String line = null;
			String[] lineComp = null;
			ArrayList<String[]> nameComponents = new ArrayList<String[]>();
			String inDir = "";
			int numFile = 1;
			
			while((line =br.readLine())!=null){
				if(line.startsWith("//")||line.length()==0) continue;	// "//" means comment text, skip
				lineComp = line.split("\\s?,?\\s+");
				if(lineComp[0].equals("inDir") && lineComp.length>1) inDir = lineComp[1];
				else{
					//process data set file name components
					nameComponents.add(lineComp);
					numFile *= lineComp.length;
				}
			}
			// create filenames
			int idx = 0;
			files = new String[numFile];
			Arrays.fill(files,  inDir+"/");
			int outRepeat = 1, inRepeat = numFile;
			for(String[] comps: nameComponents){
				idx = 0;
				inRepeat /= comps.length;
				for(int i = 0; i < outRepeat; i++){
					for(String s: comps){
						for(int j=0; j< inRepeat; j++){
							files[idx++] += s;
						}
					}
				}
				outRepeat *= comps.length;	//increase outter loop for comps.length.
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return files;
	}
	
	public static String[] getDataFileNamesFromFileList(String fileList){
		String[] files = null;
		try{
			BufferedReader br = new BufferedReader(new FileReader(fileList));
			List<String> ls = new LinkedList<String>();
			String line = null;
			while((line = br.readLine())!=null) {
				if(line.length()==0 || line.startsWith("//")) continue;
				ls.add(line);
			}
			files = new String[ls.size()];
			int idx = 0;
			for(String s: ls) files[idx++] = s;
			br.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		return files;
	}
	

	public static void outputMatrix(String fileName, double[][] m){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
			StringBuilder sb = new StringBuilder();
			for(int r= 0; r < m.length; r++){
				for(double d: m[r]) sb.append(d+"\t");
				sb.setLength(sb.length() - 1);
				bw.append(sb.toString());
				sb.setLength(0);
				if(r<m.length -1) bw.newLine();
			}
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
}
