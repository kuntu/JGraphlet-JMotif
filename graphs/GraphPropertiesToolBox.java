package graphs;
import java.util.HashMap;
import java.util.Map;

public class GraphPropertiesToolBox {
	/**
	 * obtain in/out degree sequences from edge array, the result is represented as a 2 by (N+1) matrix m. m[0][i] is the in degree of node i, where i\in [1, N], i=0 is consider as tmpvarialbe and has no meaning
	 * @param size num of nodes in the graph
	 * @param edges, array ofdirected edges
	 * @return a array int[2][N+1], inOut[0][i] is  in degree of node i, (node 0 is not counted), inOut[1][i] is the out degree of node i. 
	 */
	public static int[][] getInOutDegreeFromEdges(int size, int[][] edges){
		int[][] inOut = new int[2][size+1];
		for(int[] edge:edges){
			inOut[0][edge[1]]++;
			inOut[1][edge[0]]++;
		}
		return inOut;
	}
	
	/**
	 * return the frequencies of joint in/out degree
	 * @param inOutDeg int[2][N+1], inOutDeg[0][i] represent the in degree of node i and inOutDeg[1][i] is the out degree of node i
	 * 
	 * @return int[3][M] array, M is the number of unique in/out pair, inOutJointFreq[0] is inDegree, [1] is out and [2] is
	 * 		the frequency of (inOutJointFreq[0], inOutJointFreq[1]) pair
	 */
	public static int[][] getJointInOutFreqFromInOutDegree(int[][] inOutDeg){
		HashMap<Long, Integer> hm = new HashMap<Long, Integer>();
		long e = 0;
		for(int i=1; i<inOutDeg[0].length; i++){
			e =(long) inOutDeg[0][i];
			e = (e<<32) + inOutDeg[1][i];
			if(!hm.containsKey(e)) hm.put(e, 1);
			else hm.put(e, hm.get(e)+1);
		}
		int[][] inOutJointFreq = new int[3][hm.size()];
		int idx= 0;
		for(Map.Entry<Long, Integer> en: hm.entrySet()){
			e = en.getKey();
			inOutJointFreq[1][idx] = (int) e;
			inOutJointFreq[0][idx] = (int) (e >> 32);
			inOutJointFreq[2][idx++] = en.getValue();
		}
		return inOutJointFreq;
	}
	
	public double getExpectedProductInOut(int[][] jointInOutDegFreq){
		double res = 0;
		double size = 0;
		for(int i = 0; i< jointInOutDegFreq.length; i++){
			size += jointInOutDegFreq[i][0]*jointInOutDegFreq[i][2];
		}
		return res;
	}
}
