package msgf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

public class Histogram<T extends Comparable<T>> extends Hashtable<T,Integer> {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private T minKey = null;
	private T maxKey = null;
	private int size;
	
	public void add(T t)
	{
		if(this.get(t) == null)
			this.put(t, 1);
		else
			this.put(t, this.get(t)+1);
		if(minKey == null || minKey.compareTo(t) > 0)
			minKey = t;
		if(maxKey == null || maxKey.compareTo(t) < 0)
			maxKey = t;
		size++;
	}

	public T minKey()	
	{ 
		return minKey; 
	}
	
	public T maxKey()	
	{ 
		return maxKey; 
	}
	
	public int totalCount()
	{
		return size;
	}
	
	@Override
	public Integer get(Object key)
	{
		Integer num = super.get(key);
		if(num == null)
			return 0;
		else
			return num;
	}
	
	public void printSorted()
	{
		ArrayList<T> keyList = new ArrayList<T>(this.keySet());
		Collections.sort(keyList);
		for(T key : keyList)
			System.out.println(key+"\t"+this.get(key));
	}
	
	public void printSortedRatio()
	{
		int totalCount = totalCount();
		ArrayList<T> keyList = new ArrayList<T>(this.keySet());
		Collections.sort(keyList);
		for(T key : keyList)
		{
			System.out.println(key+"\t"+this.get(key)+"\t"+this.get(key)/(float)totalCount);
//			System.out.print(key+"\t"+this.get(key)+"\t");
//			System.out.format("%.3f\n", this.get(key)/(float)totalCount);
		}
	}
}
