package edu.ucsd.msjava.msdbsearch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.ucsd.msjava.msutil.Pair;

public class MassErrorStat {
	private List<Pair<Float, Float>> errorList; // (error, intensity)
	
	// for all peaks (absolute)
//	private float sum;
	private float mean;
//	private float median;
	private float sd;
	
	// for top 7 peaks (absolute)
//	private float sum7;
	private float mean7;
//	private float median7;
	private float sd7;

	// for all peaks (absolute)
//	private float rSum;
	private float rMean;
//	private float rMedian;
	private float rSd;
	
	// for top 7 peaks (absolute)
//	private float rSum7;
	private float rMean7;
//	private float rMedian7;
	private float rSd7;
	
	public MassErrorStat()
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
		
		List<Float> allRErrors = new ArrayList<Float>();
		List<Float> top7RErrors = new ArrayList<Float>();
		
		Collections.sort(errorList, new Pair.PairReverseComparator<Float,Float>(true));	// sort by intensities
		int rank = 0;
		for(Pair<Float,Float> errInfo : errorList)
		{
			float error = errInfo.getFirst();
			float absError = Math.abs(error);
			allErrors.add(absError);
			allRErrors.add(error);
			if(++rank <= 7)
			{
				top7Errors.add(absError);
				top7RErrors.add(error);
			}
		}
		
		mean = mean(allErrors);
		rMean = mean(allRErrors);
		sd = stdev(allErrors);
		rSd = stdev(allRErrors);
		
		mean7 = mean(top7Errors);
		rMean7 = mean(top7RErrors);
		sd7 = stdev(top7Errors);
		rSd7 = stdev(top7RErrors);
	}

	public List<Pair<Float, Float>> getErrorList() {
		return errorList;
	}

	public int size() {
		return errorList.size();
	}
	
//	public float getSum() {
//		return sum;
//	}

	public float getMean() {
		return mean;
	}
	
	public float getRMean() {
		return rMean;
	}

//	public float getMedian() {
//		return median;
//	}

	public float getSd() {
		return sd;
	}

	public float getRSd() {
		return rSd;
	}
	
//	public float getSum7() {
//		return sum7;
//	}

//	public float getRSum7() {
//		return rSum7;
//	}
	
	public float getMean7() {
		return mean7;
	}

	public float getRMean7() {
		return rMean7;
	}
	
	public float getSd7() {
		return sd7;
	}

	public float getRSd7() {
		return rSd7;
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
		double sumSq = 0;
		for(float num : numbers)
			sumSq += num*num;
		float mean = mean(numbers);
		
		float var = (float)sumSq/numbers.size() - mean*mean;
		return (float)Math.sqrt(var);
	}
}
