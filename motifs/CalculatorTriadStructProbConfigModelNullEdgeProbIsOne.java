package motifs;

/**
 * In a configuration model, the probability that u does not connect to v, 
 * 			P(!(u->v)|I_u,O_u, I_v, O_v, numEdge) = probNull = 1		 
 * 			P(u->v|I_u,O_u, I_v, O_v, numEdge)= \frac{O_uI_v}{numEdge } 
 * @author captor
 *
 */
public class CalculatorTriadStructProbConfigModelNullEdgeProbIsOne extends TriadStructProbability
		implements StructProbFromDegreeEdgeCalculator {

	@Override
	public double computeProbForStruct(double[][] deg, int numEdge, int t) {
		double res = 1, tmp = 1;
		for(int i=0; i<3; i++){
			for(int j=0; j<2; j++){
				if(deg[i][j] < triadInOutDegs[t][i][j]) return 0;
			}
		}
		int v = 0;
		for(int u =0; u<3; u++){
			for(int i = 0; i<2; i++){
				v = DyadNode2Node[u][i];
				tmp = deg[u][1]/numEdge * deg[v][0];
				if(DyadEdgeIndicator[t][u][i]){	//u->v: O_u * I_v /|E|
					if(tmp <1) res *= tmp;
					deg[u][1]--;
					deg[v][0]--;
				}else{
					if(tmp >=1) return 0;
				}
				if(res == 0) return 0;
			}
		}
		return res;
	}

}
