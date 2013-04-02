package edu.ucsd.msjava.msdbsearch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.ucsd.msjava.msutil.Pair;

public class FragmentMassErrors {
	public List<Pair<Float, Float>> getErrorList() {
		return errorList;
	}

	public float getSum() {
		return sum;
	}

	public float getMean() {
		return mean;
	}

	public float getMedian() {
		return median;
	}

	public float getSd() {
		return sd;
	}

	public float getSum7() {
		return sum7;
	}

	public float getMean7() {
		return mean7;
	}

	public float getMedian7() {
		return median7;
	}

	public float getSd7() {
		return sd7;
	}

	private List<Pair<Float, Float>> errorList; // (error, intensity)
	
	// for all peaks
	private float sum;
	private float mean;
	private float median;
	private float sd;
	
	// for top 7 peaks
	private float sum7;
	private float mean7;
	private float median7;
	private float sd7;
	
	public FragmentMassErrors()
	{
		errorList = new ArrayList<Pair<Float, Float>>();
	}
	
	public void add(Pair<Float, Float> error)
	{
		errorList.add(error);
	}
	
	public void computeStats()
	{
		List<Float> allErrors = new ArrayList<Float>();
		List<Float> top7Errors = new ArrayList<Float>();
		
		Collections.sort(errorList, new Pair.PairComparator<Float,Float>(true));	// sort by intensities
		int rank = 0;
		for(Pair<Float,Float> errInfo : errorList)
		{
			float error = errInfo.getFirst();
			allErrors.add(error);
			if(++rank <= 7)
				top7Errors.add(error);
		}
		
		sum = sum(allErrors);
		mean = mean(allErrors);
		median = median(allErrors);
		sd = stdev(allErrors);
		
		sum7 = sum(top7Errors);
		mean7 = mean(top7Errors);
		median7 = median(top7Errors);
		sd7 = stdev(top7Errors);
	}
	
	public static float sum(List<Float> numbers)
	{
		float sum = 0;
		for(float num : numbers)
			sum += num;
		return sum;
	}
	
	public float mean(List<Float> numbers)
	{
		return sum(numbers)/numbers.size();
	}
	
	public float median(List<Float> numbers)
	{
		ArrayList<Float> sorted = new ArrayList<Float>(numbers);
		Collections.sort(sorted);
		int mid = sorted.size()/2;
		if(sorted.size() % 2 == 0)
			return (sorted.get(mid-1) + sorted.get(mid))/2;
		else
			return sorted.get(mid);
	}
	
	public float stdev(List<Float> numbers)
	{
		float sumSq = 0;
		for(float num : numbers)
			sumSq += num*num;
		float mean = mean(numbers);
		
		float var = sumSq - mean*mean;
		return (float)Math.sqrt(var);
	}
}
