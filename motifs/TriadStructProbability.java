package motifs;

public class TriadStructProbability {
	protected static final int[][][] triadInOutDegs = new int[][][]{	// 16 triads: 
		
		{{0,0}, {0,0},{0,0}},	//1 null 003
		{{0,1}, {1,0},{0,0}},	//2 012
		{{1,1}, {1,1},{0,0}},	//3 102
		{{0,2}, {1,0},{1,0}},	//4 021D
		{{2,0}, {0,1}, {0,1}},	//5 021U
		{{1,1}, {0,1}, {1,0}},	//6 021c
		{{0,1}, {2,1},{1,1}},	//7 111D
		{{1,0}, {1,2},{1,1}},	//8 111U
		{{2,0}, {0,2}, {1,1}},	//9 030T
		{{1,1}, {1,1},{1,1}},	//10 030C
		{{1,1}, {2,2},{1,1}},	//11 201
		{{0,2}, {2,1},{2,1}},	//12 120D
		{{2,0}, {1,2}, {1,2}},	//13 120U
		{{1,1}, {1,2},{2,1}},	//14 120C
		{{2,1}, {1,2}, {2,2}},	//15 210
		{{2,2}, {2,2},{2,2}}	//16 300
	};
	protected static final int[][] triadNodePermutation = new int[][]{
		{0},	// 1. only 1 possible combination {0, 1,2} for type 1{ (0,0) (0,0) (0,0)}, which is triadPermutation[0], (0,0) (0,0) (0,0)
		{0, 1, 2, 3, 4, 5}, // 2. 6 possible combination, see triadPermutation 
		{0, 3, 4},	// 3. 3 possible combination, see triadPermutation 
		{0, 3, 4},	// 4. 3 6 possible combination, see triadPermutation 
		{0, 3, 4},	// 5. 3 possible combination
		{0, 1, 2, 3, 4, 5},	// 6. 6 possible combination, see triadPermutation
		{0, 1, 2, 3, 4 ,5},	// 7. 6 possible combination, see triadPermutation
		{0, 1, 2, 3, 4, 5},	// 8. 6 possible combination, see triadPermutation
		{0, 1, 2, 3, 4, 5},	// 9. 6 possible combination, see triadPermutation
		{0, 1},	// 10. 2 possible combination, see triadPermutation
		{0, 3, 4},	// 11. only 1 possible combination {0, 1,2}, which is tirad permutation 0 (0,0) (0,0) (0,0)
		{0, 3, 4},	// 12. only 1 possible combination {0, 1,2}, which is tirad permutation 0 (0,0) (0,0) (0,0)
		{0, 3, 4},	// 13. only 1 possible combination {0, 1,2}, which is tirad permutation 0 (0,0) (0,0) (0,0)
		{0, 1, 2, 3, 4, 5},	// 14. 6 possible combination, see triadPermutation
		{0, 1, 2, 3, 4, 5},	// 15. 6 possible combination, see triadPermutation
		{0},	// 16. only 1 possible combination {0, 1,2}, which is tirad permutation 0 (0,0) (0,0) (0,0)
		
	};
	protected static final int[][] triadPermutation = new int[][]{
		{0, 1, 2},
		{0, 2, 1},
		{1, 0, 2},
		{1, 2, 0},
		{2, 0, 1},
		{2, 1 ,0}
	};
	protected static final int[] NumEdgeInTriads = new int[]{
			0, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5 ,6
	};
	protected static final boolean[][][] DyadEdgeIndicator = new boolean[][][]{
	// 0-> 1	->2;	1->0	->2;	 2->0,	->1					0	1		2
		{{false, false}, {false, false}, {false, false}},	//1	{{0,0}, {0,0},{0,0}},	//1 null 003
		{{true, false}, {false, false}, {false, false}},	//2 {{0,1}, {1,0},{0,0}},	//2 012
		{{true, false}, {true, false}, {false, false}},		//3 {{1,1}, {1,1},{0,0}},	//3 102
		{{true, true}, {false, false}, {false, false}},		//4 {{0,2}, {1,0},{1,0}},	//4 021D
		{{false, false}, {true, false}, {true, false}},		//5 {{2,0}, {0,1}, {0,1}},	//5 021U
		{{false, true}, {true, false}, {false, false}},		//6 {{1,1}, {0,1}, {1,0}},	//6 021c
		{{true, false}, {false, true}, {false, true}},		//7 {{0,1}, {2,1},{1,1}},	//7 111D
		{{false, false}, {true, true}, {false, true}},		//8 {{1,0}, {1,2},{1,1}},	//8 111U
		{{false, false}, {true, true}, {true, false}},		//9 {{2,0}, {0,2}, {1,1}},	//9 030T
		{{true, false}, {false, true}, {true, false}},		//10 {{1,1}, {1,1},{1,1}},	//10 030C
		{{true, false}, {true, true}, {false, true}},		//11 {{1,1}, {2, 2},{1,1}},	//11 201
		{{true, true}, {false, true}, {false, true}},		//12 {{0,2}, {2,1},{2,1}},	//12 120D
		{{false, false}, {true, true}, {true, true}},		//13 {{2,0}, {1,2}, {1,2}},	//13 120U
		{{false, true}, {true, true}, {false, true}},		//14 {{1,1}, {1,2},{2,1}},	//14 120C
		{{false, true}, {true, true}, {true, true}},		//15 {{2,1}, {1,2}, {2,2}},	//15 210
		{{true, true}, {true, true}, {true, true}}			//16
	};
	protected static final int[][] DyadNode2Node = new int[][]{
		{1, 2}, {0, 2}, {0, 1}
	};
	private static double[][] Deg = new double[3][2];
	/**
	 * compute the probability that 3 nodes in their in/out degree consist a triad of structure type t
	 * The probability that a has no edge to v or w in the structure is computed as: 
	 * 		probNull = \prod_{j=0)^{O_u-1}(1-frac{I_v+I_w}{|E| - I_u - j})
	 * where O_u is out-degree and I_u is in-degree
	 * @param deg inOutDeg[i][0] is the in degree of the i-th node, inOutDeg[i][1] is the out-degree 
	 * @param numEdge number of edges in graph
	 * @param type structure type
	 * @return probability 
	 */
	public static double approxStructProbNoLoopNoMulti(double[][] deg, int numEdge, int type){
		double res = 0;
		for(int i=0; i<3; i++){
			for(int j=0; j<2; j++){
				if(deg[i][j] < triadInOutDegs[type][i][j]) return res;
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
				for(int j=0; j<deg[u][1]; j++){
					probNull *= 1- deg[v][0] /(numEdge - deg[u][0] - j);
				}
				if(DyadEdgeIndicator[type][u][i]){	// if u->v compute in[v]/|E|
					if(deg[v][0] == 0) return 0;
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
 
	public static double getApproxProbConfigModel(int[] nIds, int[] in, int[] out, int t, int numEdge, StructProbFromDegreeEdgeCalculator opt){
		double res = 0;
		for(int p: triadNodePermutation[t] ){
			for(int n = 0; n<3; n++){//assign in-out degrees to 3 nodes
				Deg[n][0] = in[nIds[triadPermutation[p][n]]];
				Deg[n][1] = out[nIds[triadPermutation[p][n]]];
			}
			res += opt.computeProbForStruct(Deg, numEdge, t);
		}
		return res;
	}
	
	private static double approxStructProbConfigModl(double[][] deg, int numEdge, int t){
		double res = 0;
		for(int i=0; i<3; i++){
			for(int j=0; j<2; j++){
				if(deg[i][j] < triadInOutDegs[t][i][j]) return res;
			}
		}
		int v = 0;
		res = 1;
		for(int u =0; u<3; u++){
			for(int i = 0; i<2; i++){
				v = DyadNode2Node[u][i];
				if(DyadEdgeIndicator[t][u][i]){	//u->v: O_u * I_v /|E|
					res *= deg[u][1]/numEdge * deg[v][0];
					deg[u][1]--;
					deg[v][0]--;
					numEdge--;
				}else res *= 1 - deg[u][1] / numEdge * deg[v][0];
				if(res == 0) return 0;
			}
		}		
		return res;
	}
	public static double getApproxProbConfigModel(int[] nIds, int[] in, int[] out, int t, int numEdge){
		double res = 0;
		for(int p: triadNodePermutation[t] ){
			for(int n = 0; n<3; n++){//assign in-out degrees to 3 nodes
				Deg[n][0] = in[nIds[triadPermutation[p][n]]];
				Deg[n][1] = out[nIds[triadPermutation[p][n]]];
			}
			res += approxStructProbConfigModl(Deg, numEdge, t);
		}
		return res;
	}
	//-------
	/*
	 * <1>
	 if(DyadEdgeIndicator[type][0][0]){
			res *= deg[1][0]/(numEdge - deg[0][0]);
			c++;
			deg[1][0]--;
			numEdge--;
		}
		if(DyadEdgeIndicator[type][0][1]){
			res *= deg[2][0]/(numEdge - deg[0][0]);
			c++;
			deg[2][0]--;
			numEdge--;
		}
		if(c==2) res *= c;
		//other out going edge of u
		sumIn -= c;
		for(int i=(int) c; i<deg[0][1]; i++){
			res *= (numEdge - sumIn)/(numEdge - deg[0][0]);
			numEdge--;
		}
		c = 0;
		
		//node v
		if(DyadEdgeIndicator[type][1][0]){
			res *= deg[0][0]/(numEdge - deg[1][0]);
			c++;
			deg[0][0]--;
			numEdge--;
		}
		if(DyadEdgeIndicator[type][1][1]){
			res *= deg[2][0]/(numEdge - deg[1][0]);
			c++;
			deg[2][0]--;
			numEdge--;
		}
		if(c==2) res *= c;
		//other out going edge of v
		sumIn -= c;
		for(int i=(int) c; i<deg[1][1]; i++){
			res *= (numEdge - sumIn)/(numEdge - deg[1][0]);
			numEdge--;
		}
		c = 0;
		//node w
		if(DyadEdgeIndicator[type][2][0]){
			res *= deg[0][0]/(numEdge - deg[2][0]);
			c++;
			deg[0][0]--;
			numEdge--;
		}
		if(DyadEdgeIndicator[type][2][1]){
			res *= deg[1][0]/(numEdge - deg[2][0]);
			c++;
			deg[1][0]--;
			numEdge--;
		}
		if(c==2) res *= c;
		//other out going edge of v
		sumIn -= c;
		for(int i=(int) c; i<deg[2][1]; i++){
			res *= (numEdge - sumIn)/(numEdge - deg[2][0]);
			numEdge--;
		}
		c = 0;
	 * */
}
