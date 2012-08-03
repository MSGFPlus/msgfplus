package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MgfSpectrumParser;

public class TopDownAnalysis {
	public static void main(String argv[]) throws Exception
	{
		
	}
	
	public static void tagTest() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/TopDown/sangtae_list.mgf";
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			Peptide annotation = spec.getAnnotation();
		}
	}
}
