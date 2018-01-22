package motifs;


public class StatSubGraphChange {
	public static int totalTime;
	public static int nonNullBeg;
	public static int nonNullEnd;
	public int[] editDisDistr;
	public int nonNullTime;
	public StatSubGraphChange(){
		editDisDistr = new int[7];
	}
	
	public void calChange(MotifChainNode pre, MotifChainNode cur){
		if(cur.timeStep==0) {
			nonNullBeg = 0;
			return ;
		}
		int editDis = 0;
		//pre is null, subgraph is null motif
		if(pre==null || pre.motifType ==0){
			editDis = MotifGraph.TRIAD_EDIT_DISTANCE[0][cur.motifType];
			nonNullBeg = cur.timeStep;
		}else{
			editDis = MotifGraph.TRIAD_EDIT_DISTANCE[pre.motifType][cur.motifType];
			editDisDistr[0] +=  cur.timeStep -1 -pre.timeStep;	//if subgraph remain stable in time window [pre.time, cur.time-1], edit distatnce at that period is 0
		}
		editDisDistr[editDis]++;
		if(cur.motifType == 0){	//calculate the time when subgraph is not anull motif
			nonNullTime += cur.timeStep - nonNullBeg;
			nonNullBeg = -1;
		}
	}
	
	public String toString(){
		String res = nonNullTime +":";
		int sumdis = 0;
		for(int i=1; i< editDisDistr.length; i++) sumdis+= i*editDisDistr[i];
		return res+sumdis;
	}
}
