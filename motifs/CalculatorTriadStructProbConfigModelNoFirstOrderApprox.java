package motifs;

/**
 * In a configuration model, the probability that u does not connect to v, P(!(u->v)|I_u,O_u, I_v, O_v, numEdge):
 * 		 probNull = \prod_{j=0}^{O_u-1}(1 - \frac{I_v}{numEdge - I_u - j}) 
 * @author captor
 *
 */
public class CalculatorTriadStructProbConfigModelNoFirstOrderApprox extends TriadStructProbability
		implements StructProbFromDegreeEdgeCalculator {

	@Override
	public double computeProbForStruct(double[][] deg, int numEdge, int t) {
		double res = 0;
		for(int i=0; i<3; i++){
			for(int j=0; j<2; j++){
				if(deg[i][j] < triadInOutDegs[t][i][j]) return 0;
			}
		}
		res = 1.0;
		//node u
		double probNull = 1;
		int u =0, v=0;
		for(u =0; u<3; u++){
			for(int i=0; i<2; i++){
				v = DyadNode2Node[u][i];
				//compute probNull
				probNull = 1;
				if(deg[u][1]> deg[v][0]){	// compute with probability that v's incoming edge come from nodes other than u: prod_j=1^{I_v}(1-O_u/(|E|-O_v-j))
					for(int j=0; j<deg[v][0]; j++){
						probNull *= 1- deg[u][1] /(numEdge - deg[v][0] - j);
					}
				}else{
					for(int j=0; j<deg[u][1]; j++){
						probNull *= 1- deg[v][0] /(numEdge - deg[u][0] - j);
					}
				}
				if(DyadEdgeIndicator[t][u][i]){	// if u->v compute in[v]/|E|
					res *= (1 - probNull);
					deg[u][1]--;
					deg[v][0]--;
					numEdge--;
				}else{
					res *= probNull;
				}
				if(res == 0) return 0;
			}
		}
		return res;
	}

}
