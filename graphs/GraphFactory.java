package graphs;
import java.io.*;
import java.util.*;

import randomgraph.*;


public class GraphFactory {
	
	/**
	 * generate a directed graph with int[][] as edges representation from a network data file
	 * the first line of the network data file may contains the size of the graph, in the format:
	 * 		size xxx
	 * where xxx is an integer indicate the number of nodes
	 * The rest lines contains 2 number, represent the ID of the node, in the format:
	 * 		source_node_ID		target_node_ID
	 * where the node IDs begins with 1
	 * @param file (String): the full path of file name to read in.
	 * @return
	 */
	public static GraphOfEdgeArray makeEdgeGraphFromFileOfEdgeList(String file){
		GraphOfEdgeArray graph = null;
		try{
			BufferedReader  br = new BufferedReader(new FileReader(file));
			int size = -1, minID = Integer.MAX_VALUE;
			String line = null;
			String[] edgeStr = null;
			int[] edge = null;
			LinkedList<int[]> edgeList = new LinkedList<int[]>();
			while((line = br.readLine())!=null){
				if(line.startsWith("//") || line.startsWith("#") || line.length() == 0) continue;
				edgeStr = line.split("\\s?[,\\s]+");
				if(edgeStr[0].equalsIgnoreCase("size")){
					size = Integer.parseInt(edgeStr[1]);
				}else{
					edge = new int[2];
					edge[0] = Integer.parseInt(edgeStr[0]);
					edge[1] = Integer.parseInt(edgeStr[1]);
					if(edge[0]!=edge[1]) edgeList.add(edge);
					minID = Math.min(minID, edge[0]);
					minID = Math.min(minID, edge[1]);
				}
			}
			int[][] edges = new int[edgeList.size()][];
			int idx = 0;
			for(int[] e: edgeList) edges[idx++] = e;
			edges = removeLoopAndMultiEdges(edges);
			String[] nameInt = null;
			if(size == -1){
				HashSet<Integer> hs = new HashSet<Integer>();
				for(int[] e: edges){
					hs.add(e[0]);
					hs.add(e[1]);
					size = Math.max(size, e[1]);
					size = Math.max(size, e[0]);
				}
				if(size >= 2* hs.size() || minID == 0){
					size = hs.size();
					System.out.println("\n[warning]: Node ID is too large for graph size or contains '0'. \n\tGraph size is set to the number of unique node IDs");
					HashMap<Integer, Integer> hm = new HashMap<Integer, Integer>();
					for(int i: hs) hm.put(i, hm.size()+1);
					for(int i=0; i< edges.length; i++){
						edge = edges[i];
						edge[0] = hm.get(edge[0]);
						edge[1] = hm.get(edge[1]);
					}
				}
			}
			graph = new GraphOfEdgeArray(edges, true, size);
			graph.nodeNames = nameInt;
			br.close();
		}catch(Exception e){
			e.printStackTrace();
			System.out.println("Error when creating graph");
			graph = new GraphOfEdgeArray(new int[0][2], true, 0);
		}
		return graph;
	}
	
	/**
	 * generate a directed graph with int[][] as edges representation from a network data file
	 * the first line of the network data file may contains the size of the graph, in the format:
	 * 		size xxx
	 * where xxx is an integer indicate the number of nodes
	 * The rest lines contains 2 number, represent the ID of the node, in the format:
	 * 		source_node_ID		target_node_ID
	 * where the node IDs begins with 1
	 * This function removes duplicate edges.
	 * @param file (String): the full path of file name to read in.
	 * @return
	 */
	public static GraphOfEdgeArray makeEdgeGraphFromFileOfEdgeListWithDuiplicate(String file){
		GraphOfEdgeArray graph = null;
		try{
			BufferedReader  br = new BufferedReader(new FileReader(file));
			int size = -1;
			String line = null;
			String[] edgeStr = null;
			int[] edge = null;
			HashSet<Long> edgeSet = new HashSet<Long>();
			long s =0, t = 0;
			while((line = br.readLine())!=null){
				if(line.startsWith("//") || line.startsWith("#") || line.length() == 0) continue;
				edgeStr = line.split("\\s?,?\\s+");
				if(edgeStr[0].equalsIgnoreCase("size")){
					size = Integer.parseInt(edgeStr[1]);
				}else{
					s = Long.parseLong(edgeStr[0]);
					t = Long.parseLong(edgeStr[1]);
					if(s== t) continue;
					s = (s<<32) + t;
					edgeSet.add(s);
				}
			}
			int[][] edges = new int[edgeSet.size()][2];
			int idx = 0;
			for(long l: edgeSet){
				edges[idx][1] = (int) l;
				edges[idx][0] = (int) (l>>32);
				++idx;
			}
			String[] nameInt = null;
			if(size == -1){
				HashSet<Integer> hs = new HashSet<Integer>();
				for(int[] e: edges){
					hs.add(e[0]);
					hs.add(e[1]);
					size = Math.max(size, e[1]);
					size = Math.max(size, e[0]);
				}
				if(size >= 2* hs.size()){
					System.out.println("\t[warning]: Node ID is too large for graph size for computation efficiency. \n\tGraph size is set to the number of unique node IDs");
					size = hs.size();
					nameInt = new String[size];
					size = 0;
					for(int i: hs) nameInt[size++] = String.valueOf(i);
					Arrays.sort(nameInt);
					HashMap<String, Integer> hm = new HashMap<String, Integer>();
					for(int i=0; i<nameInt.length; i++){
						hm.put(nameInt[i], i+1);
					}
					for(int i=0; i< edges.length; i++){
						edge = edges[i];
						edge[0] = hm.get(String.valueOf(edge[0]));
						edge[1] = hm.get(String.valueOf(edge[1]));
					}
				}
			}
			graph = new GraphOfEdgeArray(edges, true, size);
			graph.nodeNames = nameInt;
			br.close();
		}catch(Exception e){
			e.printStackTrace();
			System.out.println("Error when creating graph");
			graph = new GraphOfEdgeArray(new int[0][2], true, 0);
		}
		return graph;
	}
	
	
	/**
	 * 
	 * @param file
	 * @return
	 */
	public static TemporalGraphWEdgeArray getTemproalGraphWEdgeArrayFromFile(String file){
		TemporalGraphWEdgeArray g = null;
		try{
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while(line!=null && (line.startsWith("//") || line.startsWith("#") || line.length() == 0) ){	// read the firt valid line to get num of nodes and number of time snapshots
				line = br.readLine();
			}
			String[] data = line.split("\\s?,?\\s+");
			int size = Integer.parseInt(data[0]);
			int time = Integer.parseInt(data[1]), t =0;
			int numOfEdge = 0, numSelfLoop = 0;
			int[][][] tEdges = new int[time][][];
			int[][] edges = null, tmpEdges = null;
			while(t< time && (line = br.readLine())!= null){
				if(line.startsWith("//")|| line.length()==0) continue;
				data = line.split("\\s?,?\\s+");
				numOfEdge = data.length /2;
				numSelfLoop = 0;
				edges = new int[numOfEdge][2];
				for(int i = 0; i< numOfEdge; i++){
					edges[numSelfLoop][0] = Integer.parseInt(data[2*i]);
					edges[numSelfLoop][1] = Integer.parseInt(data[2*i + 1]);
					if(edges[numSelfLoop][0] == edges[numSelfLoop][1]) numSelfLoop -= 1;//remove self-loop
					++numSelfLoop;
				}
				if(numSelfLoop < numOfEdge){
					tmpEdges = new int[numSelfLoop][2];
					for(int i = 0; i< tmpEdges.length; i++){
						tmpEdges[i][0] = edges[i][0];
						tmpEdges[i][1] = edges[i][1]; 
					}
				}else{
					tmpEdges = edges;
				}
				tEdges[t] = tmpEdges;
				t++;
			}
			numOfEdge = 0;
			while(t<time){
				tEdges[t++] = new int[0][2];
				numOfEdge++;	// re-use to record number of empty snapshots
			}
			g = new TemporalGraphWEdgeArray(size, time, tEdges);
			g.numNullSnapshot = numOfEdge;	// re-use to record number of empty snapshots
			br.close();
		}catch(Exception e){
			e.printStackTrace();
			System.out.println("Create a empty temporal network.");
			g = new TemporalGraphWEdgeArray(0, 0, new int[0][0][2]);
		}
		return g;
	}
	
	public static RandomGraphJointInOutDegree getRandomGraphWSameJointIODegree(GraphOfEdgeArray graph){
		RandomGraphJointInOutDegree g = null;
		int[][] ioDegSeq = graph.getDegreeSeq();
		g = new RandomGraphJointInOutDegree(ioDegSeq);
		return g;
	}
	public static RandomGraphJoinInOutDegreeRemoveLoopMultiEdge getRandomGraphLoopAndMultiEdgeWSameJointIODegree(GraphOfEdgeArray graph){
		RandomGraphJoinInOutDegreeRemoveLoopMultiEdge g = null;
		int[][] ioDegSeq = graph.getDegreeSeq();
		g = new RandomGraphJoinInOutDegreeRemoveLoopMultiEdge(ioDegSeq);
		return g;
	}
	
	public static int[][] removeLoopAndMultiEdges(int[][] edges){
		HashSet<Long> hs = new HashSet<Long>();
		long edgeHash = 0;
		for(int[] edge: edges){
			if(edge[0] == edge[1]) continue;
			edgeHash = edge[0];
			edgeHash = (edgeHash<<32) | edge[1];
			hs.add(edgeHash);
		}
		if(hs.size() == edges.length) return edges; 
		edges = new int[hs.size()][2];
		int idx = 0;
		edgeHash = (1L << 32);
		for(long e: hs){
			edges[idx][1] = (int) (e % edgeHash);
			edges[idx][0] = (int) (e>>>32);
			idx++;
		}
		return edges;
	}
}
