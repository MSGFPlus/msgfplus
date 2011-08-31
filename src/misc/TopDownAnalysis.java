package misc;

import parser.MgfSpectrumParser;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;

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
