package mathFunctions;

import java.util.Comparator;

public class ArrayLongComparator implements Comparator<long[]> {
	public int[] cpIdx;
	public boolean desc;
	public ArrayLongComparator(int[] idxs, boolean descending){
		cpIdx = idxs;
		desc = descending;
	}
	@Override
	public int compare(long[] a, long[] b) {
		if(a.length < cpIdx.length || b.length < cpIdx.length) return 0;
		for(int i: cpIdx){
			if(a[i] == b[i]) continue;
			else return  (desc? b[i] - a[i]:a[i] - b[i]) >0 ? 1:-1;
		}
		return 0;
	}
}
