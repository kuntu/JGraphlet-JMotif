import java.io.*;
import java.util.*;


public class RandGraphGenerator {
	
	public int[][] generateRandGraph(int size, int avgDeg){
		int numEdge = size * avgDeg;
		HashSet<Long> hs = new HashSet<Long>();
		Random rnd = new Random();
		int s = 0, t =0;
		long e = 0;
		int[][] edges = new int[numEdge][2];
		for(int i=0; i< numEdge; ++i){
			s = rnd.nextInt(size) + 1;
			t = rnd.nextInt(size) + 1;
			if(s == t) --i;
			else{
				e = s;
				e = (e<<32) + t;
				if(!hs.contains(e)){
					edges[i] = new int[]{s, t};
					hs.add(e);
				}else{
					--i;
				}
			}
		}
		return edges;
	}
	
	public void saveTemporalRandomGraph(String fileName, int t, int size, int avgDeg){
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(fileName)));
			bw.write(size+"\t"+t);
			int[][] edges = null;
			StringBuilder sb = new StringBuilder();
			for(int i = 0; i<t; ++i){
				edges = generateRandGraph(size, avgDeg);
				sb.setLength(0);
				for(int[] edge: edges){
					sb.append(edge[0] + " " + edge[1] +" ");
				}
				sb.setLength(sb.length()-1);
				bw.newLine();
				bw.write(sb.toString());
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		RandGraphGenerator rgg = new RandGraphGenerator();
		rgg.saveTemporalRandomGraph("./dataSets/testData/dynGraph.txt", 4, 50, 5);
	}

}
