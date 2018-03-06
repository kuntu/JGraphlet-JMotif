package networkSampling;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

/**
	 * sample target nodes in a directed edges given their in-degree??
	 * Note: node ID begin from 1 - N, ID 0 represent null node
	 * in degree of null node 0 is also considered as 0
	 * @author kuntu
	 *
	 */
public class RandomNodeSampler {
	private ArrayList<Integer> ls;	// node IDs whose in degree is larger than 0
	private HashMap<Integer, Integer> id2Idx;
	public int[] freq;	// samp;
	private Random rnd;
	private int lastIdx;
	public RandomNodeSampler(int[] f){
		ls = new ArrayList<Integer>();
		id2Idx = new HashMap<Integer, Integer>();
		freq = f;
		for(int i = 1 ; i< freq.length; i++){
			if(freq[i]>0){
				ls.add(i);
				id2Idx.put(i, ls.size()-1);
			}
		}
		lastIdx = ls.size() - 1;
		rnd = new Random();
	}
	
	/**
	 * For source node j, sample randomly a number of target nodes that are connected by j. no loop is allowed. 
	 * Freq[i] represents the number of directed edges that can point to target node i. 
	 * @param sID
	 * @param num
	 * @return
	 */
	public int[] sampleTargetNodesForSourceNode(int sID, int num){
		int len = lastIdx+ 1;
		if(id2Idx.containsKey(sID) && id2Idx.get(sID) <= lastIdx){
			len--;
			swapNodeToIdx(sID, lastIdx);
		}
		if(len< num){
			num = len;
		}
		int[] res = new int[num];
		if(sID<0 || sID >=freq.length) sID = 0;	//sID<0 or >= networksize + 1 = freq.length
		int idx = 0, ii = 0;
			for(int i =0; i< num && i<len ; i++){
				idx = i + rnd.nextInt(len - i);	// swap value between i and idx, nodes are chosen by 1/ls.size(), not by indegree distribution
				res[ii] = ls.get(idx);
				id2Idx.put(res[ii], i);
				id2Idx.put(ls.get(i), idx);
				ls.set(idx, ls.get(i));	// replace 
				ls.set(i, res[ii]);
				if(freq[res[ii]]==0){
					System.out.println("error");
				}
				if(--freq[res[ii]] == 0){
					//remove res[i] from ls because there no more edges can point to it
					swapNodeToIdx(res[ii], lastIdx);
					lastIdx--;
					i--;
					num--;
					len--;
					if(len == lastIdx){
						swapNodeToIdx(sID, lastIdx);
					}
					if(len>lastIdx+1){
						System.out.println("error len");
					}
				}
				ii++;
			}
		return res;
	}
	private void swapNodeToIdx(int nID, int idx){
		int nIdx = id2Idx.get(nID);
		int targetID = ls.get(idx);
		ls.set(nIdx, targetID);
		ls.set(idx, nID);
		id2Idx.put(nID, idx);
		id2Idx.put(targetID, nIdx);
	}
	
	
}
