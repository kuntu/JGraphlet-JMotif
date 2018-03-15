package mathFunctions;

import java.util.Comparator;

public class ArrayComparator implements Comparator<int[]> {
	public int[] cpIdx;
	public boolean desc;
	public ArrayComparator(int[] idxs, boolean descending){
		cpIdx = idxs;
		desc = descending;
	}
	@Override
	public int compare(int[] a, int[] b) {
		if(a.length < cpIdx.length || b.length < cpIdx.length) return 0;
		for(int i: cpIdx){
			if(a[i] == b[i]) continue;
			else return desc? b[i] - a[i]:a[i] - b[i];
		}
		return 0;
	}

}
