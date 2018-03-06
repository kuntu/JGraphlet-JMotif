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
				if(line.startsWith("//")||line.startsWith("#")|| line.isEmpty()) continue;	// "//" means comment text, skip
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
				if(line.length()==0 || line.startsWith("#")||line.startsWith("//")) continue;
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
	
	public static int[][] getMatrixFromFile(String fileName){
		int[][] m = null;
		BufferedReader br = null;
		try{
			br = new BufferedReader(new FileReader(fileName));
			String line = null;
			String[] data = null;
			int N = -1, node = 0;
			while((line = br.readLine()) != null){
				if(line.startsWith("#")|| line.startsWith("//")|| line.isEmpty() ) continue;
				data = line.split("\\s?[\\s,]+");
				if(N == -1){
					N = data.length;
					m = new int[N][N];
				}
				for(int i=0; i<N; i++){
					m[node][i] = Integer.parseInt(data[i]);
				}
				++node;
			}
			br.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		return m;
	}
	
	public static int[][] convertAdjacentMatrixToEdgeList(int[][] m){
		int numEdge = 0;
		for(int[] row: m){
			for(int i: row){
				if(i==1) numEdge++;
			}
		}
		int node = 1;
		int[][] edges = new int[numEdge][2];
		numEdge = 0;
		for(int i = 0; i< m.length; i++){
			node = i + 1;
			for(int j = 0; j < m[i].length; j++){
				if(m[i][j] >0 && i != j){
					edges[numEdge][0] = node;
					edges[numEdge][1] = j+1;
					++numEdge;
				}
			}
		}
		return edges;
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
	
	public static void outputIntegerMatrix(String fileName, int[][] m){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
			StringBuilder sb = new StringBuilder();
			for(int r= 0; r < m.length; r++){
				for(int d: m[r]) sb.append(d+"\t");
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
	
	public static void convertAdjMatFileToEdgeList(String fileName){
		int[][] m = getMatrixFromFile(fileName);
		int[][] edges = convertAdjacentMatrixToEdgeList(m);
		fileName = fileName.substring(0, fileName.lastIndexOf("_adj")) + "EdgeList.txt";
		outputIntegerMatrix(fileName, edges);
	}
	
	public static void convertStrIDToIntInFile(String fileName, boolean isEgo){
		try{
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			int idx = fileName.lastIndexOf('.');
			String outfile = fileName.substring(0, idx) + (isEgo?"Ego":"") + "Formated.txt";
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			String line = null;
			String[] data = null;
			HashMap<String, Integer> map = new HashMap<String, Integer>();
			int s = 0, t = 0;
			while((line = br.readLine()) != null){
				if(line.startsWith("//") || line.startsWith("#") || line.isEmpty()) continue;
				data = line.split("\\s?[,\\s]+");
				if(map.containsKey(data[0])) s = map.get(data[0]);
				else{
					s = map.size() + 1;
					map.put(data[0], s);
				}
				if(map.containsKey(data[1])) t = map.get(data[1]);
				else{
					t = map.size() + 1;
					map.put(data[1], t);
				}
				bw.write(s + "\t" + t+"\n");
			}
			if(isEgo){
				s = map.size() + 1;
				for(int i: map.values())
					bw.write(s + "\t" + i+"\n");
			}
			br.close();
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
}
