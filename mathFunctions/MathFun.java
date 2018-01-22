package mathFunctions;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;


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
	
	/**
	 * Sample k elements from an array with proportion to their weights
	 * @param weight
	 * @param num
	 * @return
	 */
	public static int[] AChaoWeightedResorviorSampling(double[] weight, int num){
		int[] res = new int[num];
		Random rnd = new Random();
		if(num >= weight.length) num = weight.length;
		int idx = 0;
		double sum = 0;
		while(idx < num) {
			sum += weight[idx]/num;
			res[idx++] = idx+1;
		}
		while(idx<weight.length){
			sum += weight[idx]/num;
			if(rnd.nextDouble()< weight[idx]/sum){
				res[rnd.nextInt(num)] = idx + 1;
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
}
