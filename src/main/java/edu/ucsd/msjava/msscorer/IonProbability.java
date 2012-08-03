package edu.ucsd.msjava.msscorer;

import java.util.HashSet;
import java.util.Iterator;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.IonType;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Reshape;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.IonType.PrefixIon;

public class IonProbability {
	private Iterator<Spectrum> itr;
	private Reshape filter;
	private IonType[] ions;
	private Tolerance tol;
	private boolean onePerPep = false;
	private int numAllSegments = 1;
	private int targetSegment = 0;
	
	public IonProbability(Iterator<Spectrum> itr, IonType[] ions, Tolerance tol)
	{
		this.itr = itr;
		this.ions = ions;
		this.tol = tol;
	}
	
	public IonProbability segment(int targetSegment, int numAllSegments)
	{
		this.numAllSegments = numAllSegments;
		this.targetSegment = targetSegment;
		return this;
	}
	
	public IonProbability filter(Reshape filter)
	{
		this.filter = filter;
		return this;
	}
	
	public IonProbability onePerPeptide(boolean isOnePerPep)
	{
		this.onePerPep = isOnePerPep;
		return this;
	}
	
	public float[] getIonProb()
	{
		float[] ionProbArr = new float[ions.length];
		int[] numObservedPeaks = new int[ions.length];
		int[] numMissingPeaks = new int[ions.length];
		HashSet<String> pepSet = null;
		if(onePerPep)
			pepSet = new HashSet<String>();
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			if(filter != null)
				spec = filter.apply(spec);
			Peptide pep = spec.getAnnotation();
			if(pep == null)
				continue;
			
			if(onePerPep)
			{
				String pepStr = spec.getAnnotationStr();
				if(pepSet.contains(pepStr))
					continue;
				else
					pepSet.add(pepStr);
			}
			
			int index = -1;
			for(IonType ion : ions)
			{
				index++;
				if(ion instanceof PrefixIon)
				{
					double prm = 0;
					for(int i=0; i<pep.size()-1; i++)
					{
						prm += pep.get(i).getMass();
						float mz = ion.getMz((float)prm);
						if(numAllSegments > 1)
						{
							int segNum = (int)(mz/spec.getParentMass()*numAllSegments);
							if(segNum >= numAllSegments)
								segNum = numAllSegments-1;
							if(segNum != targetSegment)
								continue;
						}
						
						if(spec.getPeakByMass(mz, tol) != null)
							numObservedPeaks[index]++;
						else
							numMissingPeaks[index]++;
					}
				}
				else
				{
					double srm = 0;
					for(int i=0; i<pep.size()-1; i++)
					{
						srm += pep.get(pep.size()-1-i).getMass();
						float mz = ion.getMz((float)srm);
						if(numAllSegments > 1)
						{
							int segNum = (int)(mz/spec.getParentMass()*numAllSegments);
							if(segNum >= numAllSegments)
								segNum = numAllSegments-1;
							if(segNum != targetSegment)
								continue;
						}
						if(spec.getPeakByMass(mz, tol) != null)
						{
							numObservedPeaks[index]++;
//							if(ion.getName().equals("y2-H3PO4"))
//								System.out.println("Debug");
						}
						else
							numMissingPeaks[index]++;
					}
				}
			}
		}
		
		for(int i=0; i<ions.length; i++)
		{
			if(numObservedPeaks[i]+numMissingPeaks[i] == 0)
				ionProbArr[i] = 0;
			else
				ionProbArr[i] = numObservedPeaks[i]/(float)(numObservedPeaks[i]+numMissingPeaks[i]);
		}
		return ionProbArr;
	}
}
