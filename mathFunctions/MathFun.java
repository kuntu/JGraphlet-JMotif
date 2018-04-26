package mathFunctions;
import java.util.*;


public class MathFun {
	
	
	/**
	 * copy census array to copy 
	 * @param copy
	 * @param org
	 * @return
	 */
	public static int[] cloneCensus(int[] copy, int[] org){
		if(copy==null || copy.length!=org.length) copy = new int[org.length];
		for(int i=0; i<org.length; i++){
			copy[i] = org[i];
		}
		return copy;
	}
	
	/**
	 * compute KL divergence, when there is 0 frequency, apply Baysian prior instead of removing that part
	 * @param preStep
	 * @param curStep
	 * @return
	 */
	public static double KLDiversionFromFreqWithBayesPrior(int[] preStep, int[] curStep){
		double res = 0;
		double sum1 = 0, sum2 = 0;
		double p1 = 0;
		double prior = 0;
		boolean bs = false;
		if(preStep.length!=curStep.length){
			System.out.println("Distribution size not matched");
			return 1;
		}
		for(int i=0; i<preStep.length; i++){
			sum1 += preStep[i];
			sum2 += curStep[i];
			if(preStep[i]==0 || curStep[i] ==0) bs = true;
		}
		if(bs){
			prior = 1;
			sum1 += preStep.length;
			sum2 += curStep.length;
		}
		for(int i=0; i<preStep.length; i++){
			p1 = (preStep[i]+prior) / sum1;
			res += p1*(Math.log(p1) - Math.log((curStep[i]+prior)/sum2));
		}
		return res;
	}
	
	public static double KLDiversionFromNonZeroFreq(int[] preStep, int[] curStep){
		double res = 0;
		double sum1 = 0, sum2 = 0;
		double p1 = 0;
		if(preStep.length!=curStep.length){
			System.out.println("Distribution size not matched");
			return 1;
		}
		for(int i=0; i<preStep.length; i++){
			sum1 += preStep[i];
			sum2 += curStep[i];
		}
		for(int i=0; i<preStep.length; i++){
			p1 = preStep[i] / sum1;
			res += p1*(Math.log(p1) - Math.log(curStep[i]/sum2));
		}
		return res;
	}
	
	/**
	 * randomly sample len elements from array and put them at the first len position in the array, in O(len) time. good for array.length is known, 
	 * @param array
	 * @param len
	 */
	public static void durstenfeldShuffle(int[] array, int len){
		Random rnd = new Random();
		int idx = 0;
		int tmp = 0;
		len = Math.min(array.length, len);
		for(int i=0; i<len; i++){
			idx = i+ rnd.nextInt(array.length - i);
			tmp = array[i];
			array[i] = array[idx];
			array[idx] = tmp;
		}
	}
	/**
	 * 
	 * @param m		(int[][]) matrix
	 * @param col
	 * @param len
	 */
	public static void durstenfeldShuffleMatrixColumn(int[][] m, int col, int len){
		len = Math.min(len, m.length);
		Random rnd = new Random();
		int idx = 0, tmp = 0;
		for(int i=0; i<len; i++){
			idx = i+rnd.nextInt(m.length - i);
			tmp = m[idx][col];
			m[idx][col] = m[i][col];
			m[i][col] = tmp;
		}
	}
	
	public static void durstenfeldShuffleMatrixMultiColumn(long[][] m, int[] colIdx, int len){
		long[] tmp = new long[colIdx.length];
		Random rnd = new Random();
		if(len > m.length) len = m.length;
		int idx = 0;
		for(int i = 0; i< len; ++i){
			idx = rnd.nextInt(len - i);
			for(int j = 0; j < tmp.length; j++){
				tmp[j] = m[i][colIdx[j]];
				m[i][colIdx[j]] = m[idx][colIdx[j]];
				m[idx][colIdx[j]] = tmp[j];
			}
		}
	}
	
	/**
	 * reservior sampling with k elements
	 * @param totalNum
	 * @param k
	 * @return
	 */
	public static long[] resorviorSampling(long totalNum, int k){
		if(k > totalNum) k = (int) totalNum;
		long[] res = new long[k];
		int idx = 0;
		int pos =  -1;
		Random rnd = new Random();
		while(idx< k){
			res[idx] = idx;
			idx++;
		}
		while(idx < totalNum){
			pos = rnd.nextInt(idx+1);
			if(pos < k){
				res[pos] = idx;
			}
			++idx;
		}
		return res;
	}
	
	/**
	 * Sample k elements from an array with proportion to their weights
	 * @param weight
	 * @param k
	 * @return
	 */
	public static int[] AChaoWeightedResorviorSampling(double[] weight, int k){
		int[] res = new int[k];
		Random rnd = new Random();
		if(k >= weight.length) k = weight.length;
		int idx = 0;
		double sum = 0;
		while(idx < k) {
			sum += weight[idx]/k;
			res[idx++] = idx+1;
		}
		while(idx<weight.length){
			sum += weight[idx]/k;
			if(rnd.nextDouble()< weight[idx]/sum){
				res[rnd.nextInt(k)] = idx + 1;
			}
			idx++;
		}
		return res;
	}
	
	public static int[] AChaoWeightedResorviorSampling(ArrayList<Integer> elements, ArrayList<Double> weight, int num){
		int[] res = new int[num];
		
		int idx = 0, len = elements.size();
		double sum = 0, w = 0;;
		while(idx < num){
			res[idx] = idx;
			sum += weight.get(idx) / num;
			idx++;
		}
		Random rnd = new Random();
		while(idx < len){
			w = weight.get(idx);
			sum += w / num;
			if(rnd.nextDouble() < w / sum){
				res[rnd.nextInt(num)] = idx;
			}
		}
		Arrays.sort(res);
		for(int i = res.length -1; i>=0; i++){
			if(weight.get(res[i]) == 1){
				weight.remove(res[i]);
			}else weight.set(res[i], weight.get(res[i]) -1);
			res[i] = elements.get(res[i]);
		}
		return res;
	}
	
	/**
	 * sample k elements from 0 ~ n-1 with O(k) space and O(k) time
	 * @param n
	 * @param k
	 * @return
	 */
	public static int[] sampleKIntfromN_withNoReplacement(int n, int k){
		int[] res = null;
		if(k > n){
			res = new int[n];
			for(int i = 0; i< n; i++) res[i] = i;
			return res;
		}
		res = new int[k];
		int tmp = 0, curVal;
		Random rnd = new Random();
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		
		for(int i = 0; i< k; i++){
			if(map.containsKey(i)) curVal = map.get(i); //need to remember to swap value, i could be chosen in previous steps
			else curVal = i;
			res[i] = i + rnd.nextInt(n);
			if(map.containsKey(res[i])){
				tmp = map.get(res[i]);
				map.put(res[i], curVal);
				res[i] = tmp;
			}else map.put(res[i], curVal);
			--n;
		}
		return res;
	}
}
