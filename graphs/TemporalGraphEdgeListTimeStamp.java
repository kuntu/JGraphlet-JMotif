package graphs;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import motifs.Motif;

public class TemporalGraphEdgeListTimeStamp extends BasicGraph implements  GraphProperty, Motif {
	long[][] timeStampEdges;
	public TemporalGraphEdgeListTimeStamp(long[][] edgesOfTimeStamp) {
		timeStampEdges = edgesOfTimeStamp;
		
	}
	@Override
	public long[] getMotifFreq(int motifSize) {
		return null;
	}

	@Override
	public int[][] getDegreeSeq() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int[][] getDegreeFreq() {
		// TODO Auto-generated method stub
		return null;
	}
	
	public List<long[]> trackNodePairThreeEdgeMotif(long timeWin){
		//long[] 0 -3: count; idx; preEdgeDir; firstEdgeTimeStamp, SecEdgeTimeStamp
		//e1=e2=e3:0 (ei is the direction of edge i
		//e1=e2!=e3: 2
		//e1!=e2, e1==e2: 1
		//e1!=e2, e2!= e3: 3
		HashMap<Long, long[]> nodepairMotifMap = new HashMap<Long, long[]>();
		int idx = 0;
		long key = 0, dir = 0;
		long[] nodePair = null;
		for(long[] edge: timeStampEdges){
			key = GraphPropertiesToolBox.edgeToKey(edge[0], edge[1]);
			nodePair = nodepairMotifMap.get(key);
			if(nodePair == null){
				nodePair = new long[8];
				nodePair[6] = - timeWin - 1;	//timestamp for e1
				nodePair[7] = nodePair[6];		//timestamp for e2
				nodePair[5] = 1;	// 1 for edge direction small ID -> large ID; -1 for the other
				nodepairMotifMap.put(key, nodePair);
			}
			dir = (edge[0]<edge[1]? 1: -1);
			idx >>= 1;
			if(dir * nodePair[5] < 0) nodePair[4] += 2;
			nodePair[5] = dir;
			if(edge[2]-nodePair[6] < timeWin){
				idx = (int) nodePair[4];
				++nodePair[idx];
			}
			nodePair[6] = nodePair[7];
			nodePair[7] = edge[2];
			if(nodePair[6] == nodePair[7]) nodePair[6] = - timeWin - 1;
		}
//		res = new long[nodepairMotifMap.size()][6];
		List<long[]> ls = new LinkedList<long[]>();
		idx = 0;
		long[] res = new long[6];
		for(Map.Entry<Long, long[]> en: nodepairMotifMap.entrySet()){
			nodePair = en.getValue();
			idx = 2;
			key = 0;
			for(int i = 0; i < 4; ++i){
				res[i+2] = nodePair[i];
				key += nodePair[i];
			}
			if(key == 0) continue;	// only track node pairs that form a three edge temporal motif
			nodePair = GraphPropertiesToolBox.keyToEdgeLong(en.getKey());
			key = 0;
			res[0] = nodePair[0];
			res[1] = nodePair[1];
			ls.add(res);
			res = new long[6];
		}
		return ls;
	}
	@Override
	public long[][] getNodeMotifFreq(int motifSize) {
		// TODO Auto-generated method stub
		return null;
	}
}
