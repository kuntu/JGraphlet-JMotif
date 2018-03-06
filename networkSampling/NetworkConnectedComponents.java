package networkSampling;
import java.io.*;
import java.util.*;
import java.util.Map.Entry;


public class NetworkConnectedComponents {
	private HashMap<Integer, HashSet<Integer>> components;
	private HashMap<Integer, Integer> par;
	private HashMap<Integer, Integer> rank;
	public int[][] mapping;
	public int[] componentSize;
	
	public HashMap<Integer, HashSet<Integer>> getComponentsFromFile(String fileName){
		try{
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			par = new HashMap<Integer, Integer>();
			rank = new HashMap<Integer, Integer>();
			String line = null;
			String[] data = null;
			while((line = br.readLine())!=null){
				data = line.split(" ");
				union((int) Integer.parseInt(data[0]), (int) Integer.parseInt(data[1]) );
			}
			br.close();
			
		}catch(Exception e){
			e.printStackTrace();
		}
		setComponentsFromUnionFind();
		return components;
	}
	
	public void setComponentsFromUnionFind(){
		components = new HashMap<Integer, HashSet<Integer>>();
		int setID = 0;
		HashSet<Integer> neig = null;
		for(Entry<Integer, Integer> en: par.entrySet()){
			setID = find(en.getValue());
			if(!components.containsKey(setID)){
				neig = new HashSet<Integer>();
				components.put(setID, neig);					
			}else neig = components.get(setID);
			neig.add(en.getKey());
		}
	}
	
	private void union(int a, int b){
		int ap = find(a);
		int bp = find(b);
		if(ap != bp){
			if(!rank.containsKey(ap)){
				if(!rank.containsKey(bp)){
					par.put(bp, ap);
					rank.put(ap, 2);
				}else{
					par.put(ap, bp);
					rank.put(bp,rank.get(bp)+1);
				}
			}else{
				if(!rank.containsKey(bp)){
					par.put(bp, ap);
					rank.put(ap, rank.get(ap)+1);
				}else if (rank.get(ap)>=rank.get(bp)){
					par.put(bp, ap);
					rank.put(ap, rank.get(ap)+rank.get(bp));
				}else{
					par.put(ap, bp);
					rank.put(bp,rank.get(bp)+rank.get(ap));
				}
			}
		}
	}
	private Integer find(int a){
		if(!par.containsKey(a)) {
			par.put(a, a);
			return a;
		}
		int p = par.get(a);
		int g = par.get(p);
		while(p!= g){
			par.put(a, g);
			a = p;
			p = g;
			g= par.get(p);
		}
		return p;
	}
	
	//edge graph
	public HashMap<Integer, HashSet<Integer>> getComponentsFromEdgeListsFile(int size, int time, int[][][] edgelists, boolean skipZero){
		int len = skipZero ? size+1:size;
		int[] pars = new int[len];
		int[] ranks= new int[len];
		for(int i=0;i<len; i++){
			pars[i] =i;
			ranks[i] = 1;
		}
		for(int[][] graph: edgelists){
			for(int[] edge: graph){
				union(edge[0], edge[1], pars, ranks);
			}
		}
		setComponentFromUnionFind(pars, skipZero);
		return components;
	}
	
	public void setComponentFromUnionFind(int[] pars, boolean skipZero){
		components = new HashMap<Integer, HashSet<Integer>>();
		int idx = skipZero? 1:0;
		int setID = 0;
		HashSet<Integer> hs = new HashSet<Integer>();
		while(idx<pars.length){
			setID = find(idx, pars);
			hs = components.get(setID);
			if(hs == null){
				hs = new HashSet<Integer>();
				components.put(setID, hs);
			}
			hs.add(idx);
			idx++;
		}
	}
	
	private void union(int a, int b, int[] pars, int[] ranks){
		a = find(a, pars);
		b = find(b, pars);
		if(a == b) return;
		if(ranks[a]>= ranks[b]){
			pars[b] = a;
			ranks[a] += ranks[b];
		}else{
			pars[a] = b;
			ranks[b] += ranks[a];
		}
	}
	private int find(int a, int[] pars){
		int p = pars[a];
		int g = pars[p];
		while(p != pars[p]){
			pars[a] = pars[p];
			a = p;
			p = pars[p];
		}
		return p;
	}
	
	public int[][] nodeIDMapping( boolean skipZero){
		int len = 0;
		for(HashSet<Integer> hs: components.values()){
			len += hs.size();
		}
		if(skipZero) len++;
		mapping = new int[len][2];
		componentSize = new int[components.size()];
		int setID = 0, newNID= skipZero? 1:0;
		for(HashSet<Integer> hs: components.values()){
			componentSize[setID] = hs.size();
			newNID= skipZero? 1:0;
			for(int nID: hs){
				mapping[nID][0] = setID;
				mapping[nID][1] = newNID++;
			}
			setID++;
		}
		return mapping;
	} 
	
	public static void main(String[] args){
		NetworkConnectedComponents ncc = new NetworkConnectedComponents();
		HashMap<Integer, HashSet<Integer>> components = ncc.getComponentsFromFile("3colNetwork.txt");
		for(int key: components.keySet()){
			for(int n: components.get(key)) System.out.print(n+", ");
			System.out.println();
		}
		
		
	}
}
