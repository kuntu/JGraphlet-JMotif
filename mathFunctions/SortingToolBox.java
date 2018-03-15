package mathFunctions;

import java.util.Arrays;

public class SortingToolBox {
	public static void sortMatrixColumnInplace(int[][] m, int[] cpIdx, boolean desc){
		if(m.length == 0 || m[0].length == 0) return;
		int r = m.length, c = m[0].length;
		int[][] mtr = new int[c][r];
		for(int i = 0; i<r; ++i){
			for(int j = 0; j < c; ++j){
				mtr[j][i] = m[i][j];
			}
		}
		Arrays.sort(mtr, new ArrayComparator(cpIdx, desc));
		for(int i = 0; i<r; ++i){
			for(int j = 0; j < c; ++j){
				m[i][j] = mtr[j][i];
			}
		}
	}
}
