package msgf;

import java.util.ArrayList;

public class LinearCalibration {
	ArrayList<Float> x;
	ArrayList<Float> y;
	float slope;
	float intercept;
	boolean isUpdated;
	
	public LinearCalibration() 
	{
		isUpdated = false;
		x = new ArrayList<Float>();
		y = new ArrayList<Float>();
	}
	
	public float predict(float x)
	{
		if(!isUpdated)
			update();
		return x*slope+intercept;
	}
	
	public float getSlope()
	{
		if(isUpdated)
			return slope;
		else
		{
			update();
			return slope;
		}
	}

	public float getIntercept()
	{
		if(isUpdated)
			return intercept;
		else
		{
			update();
			return intercept;
		}
	}
		
	public void addData(float x, float y)
	{
		this.x.add(x);
		this.y.add(y);
		isUpdated = false;
	}
	
	private void update() 
	{
		float sumXSq = 0;
		float sumX = 0;
		float sumY = 0;
		float sumXY = 0;
		if(x.size() < 2)
		{
			slope = 1;
			intercept = 0;
			return;
		}
		for(int i=0; i<x.size(); i++)
		{
			sumXSq += x.get(i)*x.get(i);
			sumX += x.get(i);
			sumY += y.get(i);
			sumXY += x.get(i)*y.get(i);
		}
		slope = (x.size()*sumXY-sumX*sumY)/(x.size()*sumXSq-sumX*sumX);
		intercept = (sumY-slope*sumX)/x.size();
		isUpdated = true;
	}
	

}
