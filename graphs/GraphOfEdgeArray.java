package graphs;

import java.util.HashMap;

import motifs.Motif;

public class GraphOfEdgeArray extends BasicGraph implements GraphProperty, Motif{
	public String[] nodeNames;
	public int[][] edges;
	public GraphOfEdgeArray(int[][] e, boolean direct, int s){
		edges = e;
		this.directed = direct;
		size =s;
	}
	@Override
	public int[][] getDegreeSeq() {
		return GraphPropertiesToolBox.getInOutDegreeFromEdges(size, edges);
	}
	@Override
	public int[][] getDegreeFreq() {
		return GraphPropertiesToolBox.getJointInOutFreqFromInOutDegree(GraphPropertiesToolBox.getInOutDegreeFromEdges(size, edges));
	}
	@Override
	public long[] getMotifFreq(int motifSize) {
		GraphNodePairVal gnpv = new GraphNodePairVal(size, directed, edges);
		return gnpv.getMotifFreq(motifSize);
	}
	
	public GraphOfEdgeArray removeNullNodes(){
		HashMap<Integer, Integer> hm = new HashMap<Integer, Integer>();
		int[][] newEdges = new int[edges.length][2];
		for(int i = 0; i<edges.length; i++){
			for(int j=0; j< 2; j++){
				if(!hm.containsKey(edges[i][j])){
					newEdges[i][j] = hm.size() + 1;
					hm.put(edges[i][j], newEdges[i][j]);
				}else{
					newEdges[i][j] = hm.get(edges[i][j]);
				}
			}
		}
		return new GraphOfEdgeArray(newEdges, directed, hm.size());
	}
}
