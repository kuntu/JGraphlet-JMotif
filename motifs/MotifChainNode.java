package motifs;

public class MotifChainNode implements Comparable<MotifChainNode> {
	
	public int timeStep;
	public int motifType;
	public MotifChainNode(int t, int type){
		timeStep = t;
		motifType = type;
	}
	
	public String toString(){
		return timeStep +":"+ motifType+" ";
	}
	
	@Override
	public int compareTo(MotifChainNode b){
		return timeStep - b.timeStep;
	}
}
