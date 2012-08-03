package edu.ucsd.msjava.msscorer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ListStat {
	private List<Float> data;
	private List<Float> sortedData;
	
	public ListStat(List<Float> data)
	{
		this.data = data;
	}
	
	public ListStat(float[] data)
	{
		this.data = new ArrayList<Float>();
		for(float val : data)
			this.data.add(val);
	}
	
	public float mean()
	{
		double sum = 0;
		for(Float val : data)
			sum += val;
		return (float)sum/data.size();
	}
	
	public float median()
	{
		return percentile(0.5f);
	}
	
	public float percentile(float percentile)
	{
		if(sortedData == null)
			sortedData = new ArrayList<Float>(data);
		Collections.sort(sortedData);
		float index = (data.size()-1)*percentile;
		int lowIndex = (int)Math.floor(index);
		int highIndex = (int)Math.ceil(index);
		if(lowIndex == highIndex)
			return sortedData.get(lowIndex);
		else
			return (index-lowIndex)*sortedData.get(lowIndex)+(highIndex-index)*sortedData.get(highIndex);
	}
	
	public float stdev()
	{
		float sumSq = 0;
		for(Float val : data)
			sumSq += val*val;
		float mean = mean();
		return sumSq/data.size()-mean*mean;
	}
	
	public float getSignalThreshold(float signalToNoiseRatio)
	{
		float median = median();
		ArrayList<Float> dataWithoutOutliers = new ArrayList<Float>();
		for(Float val : data)
			if(val < median*signalToNoiseRatio)
				dataWithoutOutliers.add(val);
		
		ListStat newList = new ListStat(dataWithoutOutliers);
		float newMedian = newList.median();
		return newMedian*signalToNoiseRatio;
	}
	
	public float[] getOutlieres(float signalToNoiseRatio)
	{
		float newMedian = getSignalThreshold(signalToNoiseRatio);
		ArrayList<Float> outlierList = new ArrayList<Float>();
		for(Float val : data)
			if(val >= newMedian*signalToNoiseRatio)
				outlierList.add(val);
		
		float[] outliers = new float[outlierList.size()];
		int index = -1;
		for(Float val : outlierList)
			outliers[++index] = val;
			
		return outliers;
	}
	            
}
