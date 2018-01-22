package graphs;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import motifs.*;

public class TemporalGraphNodePairVal extends TemporalNetwork implements DynamicMotif {
	
	public GraphNodePairVal[] graphs;
	public int[] edgeNum;
	public GraphNodePairVal[] preMissingEdges;
	public TemporalGraphNodePairVal(int size, int[][][] edgeGraphs, boolean dir){
		this.size = size;
		this.time = edgeGraphs.length;
		edgeNum = new int[time];
		for(int t=0; t<time; t++){
			graphs[t] = new GraphNodePairVal(size, dir, edgeGraphs[t]);
			edgeNum[t] = edgeGraphs[t].length;
		}
		calMissingEdgesInNextTimeStep();
	}
	
	private void calMissingEdgesInNextTimeStep(){
		int[][] tmp = new int[][]{};
		HashMap<Integer, Integer> preEdgeType = null;
		HashMap<Integer, Integer> curEdgeType = null;
		HashMap<Integer, Integer> hm = null;
		for(int t=0; t<time-1; t++){	// prepare for edges that missing in the next time step
			preMissingEdges[t+1] = new GraphNodePairVal(size, graphs[t+1].directed, tmp);
			for(int n: graphs[t].nodes.keySet()){
				preEdgeType = graphs[t].nodes.get(n);
				curEdgeType = graphs[t+1].nodes.get(n);
				if(curEdgeType == null) {
					preMissingEdges[t+1].nodes.put(n, preEdgeType);
					continue;
				}
				hm = new HashMap<Integer, Integer>();
				for(int u: preEdgeType.keySet()){
					if(!curEdgeType.containsKey(u)){
						hm.put(u, preEdgeType.get(u));
					}
				}
				if(!hm.isEmpty()) preMissingEdges[t+1].nodes.put(n, hm);
			}
		}
	}
	
	@Override
	public long[][] getMotifFreq(int numNode) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int[][] getMotfiTransition(int numNode) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public int[][] getAllMotifFreq(int subGraphSize){
		int[][] res = null;
		if(subGraphSize ==3){
			res = new int[time][16];
			
		}
		
		return res;
	}
	public void getAllTriadFreq(long[][] allFreqs){
		long[] tmp = graphs[0].getMotifFreq(3);
		for(int i=0;i < tmp.length; i++) allFreqs[0][i] = tmp[i];
		HashMap<Integer, Integer> neig = null;
		HashMap<Integer, Integer> preNeig = null;
		int nuVal=-1, preNUVal=-1, nvVal = -1, preNVVal = -1, uvVal = -1, preUVVal=-1;
		int n = 0, motifType = -1;
		HashSet<Integer> hs = new HashSet<Integer>();
		for(int t = 1; t < time; t++){
			for(int i=1; i<16; i++) allFreqs[t][i] = allFreqs[t-1][i];
			//new appear edges and changing edges at time t
			for(Entry<Integer, HashMap<Integer, Integer>> entry: graphs[t].nodes.entrySet()){
				n = entry.getKey();
				neig = entry.getValue();
				preNeig = graphs[t-1].nodes.get(n);
				hs.clear();
				hs.addAll(neig.keySet());
				for(int u: hs){
					nuVal = graphs[t].getPairVal(n, u);
					preNUVal = graphs[t-1].getPairVal(n, u);
					if(nuVal == preNUVal) continue;	//unchanged edge, skip
					if(preNUVal==0 ){	// appearing edge e(n,u)
						for(int v: hs){
							uvVal = graphs[t].getPairVal(u, v);
							if(u==v|| (n>u && uvVal > 0)) continue;
							else{	//u!=v && (n<u || !e(u,v))
								nvVal = graphs[t].getPairVal(n, v);
								preNVVal = graphs[t-1].getPairVal(n, v);
								preUVVal = graphs[t-1].getPairVal(u, v);
								if(v< n && uvVal >0 && ( preNVVal!=nvVal|| preUVVal!=uvVal )) continue;
								//current motiftype
								motifType = GraphNodePairVal.getTriadTypeFromStructureVal(nuVal, nvVal, uvVal);
								allFreqs[t][motifType]++;
								//remove previous motiftype
								if(preNVVal>0 && preUVVal>0) motifType = GraphNodePairVal.getTriadTypeFromStructureVal(nuVal, nvVal, uvVal);
								else if(preNVVal==0 && preUVVal == 0) continue;
								else motifType = (preNVVal+preUVVal==3)? 3:2;
								allFreqs[t][motifType]--;
							}
						}
						// count dyads
						if(n<u){
							motifType = nvVal==3? 3:2;
							
							allFreqs[t][motifType] += graphs[t].size - hs.size() - graphs[t].nodes.get(u).size();
							for(int v:graphs[t].nodes.get(u).keySet()){
								if(hs.contains(v)) allFreqs[t][motifType]++;
							}
						}
					}else{	// changing edge
						nuVal = neig.get(u);
						if(nuVal == preNeig.get(u)) continue;
						
					}
						
					//check missing edges
						
				}
			}
			//missing edges at time t
		}
	}
	private void updateTriadFreq(int[][] allFreqs, int t){
		
	}

	@Override
	public double[] getAvgMotifFreq(int motifSize) {
		// TODO Auto-generated method stub
		return null;
	}
}
