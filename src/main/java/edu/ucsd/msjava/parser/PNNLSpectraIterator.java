package edu.ucsd.msjava.parser;

import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Iterator;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;

public class PNNLSpectraIterator extends SpectraIterator {

	private HashMap<Integer,ActivationMethod> scanNumActMethodMap;
	
	public PNNLSpectraIterator(String fileName) throws FileNotFoundException 
	{
		super(fileName, new PNNLSpectrumParser());
		scanNumActMethodMap = PNNLSpectrumParser.getScanTypeMap(fileName);	
	}

	@Override
	public Spectrum next() 
	{
		if(scanNumActMethodMap == null)
			return super.next();
		
		Spectrum spec = super.next();
		ActivationMethod method = scanNumActMethodMap.get(spec.getScanNum());
		if(method != null)
			spec.setActivationMethod(method);
		return spec;
	}
	
	public static void main(String argv[]) throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Test/Matt/QC_Shew_11_03_200ng_4_23Aug11_Hawk_11-05-04p_dta.txt";
		PNNLSpectraIterator itr = new PNNLSpectraIterator(fileName);
		Iterator<Spectrum> specItr = itr.iterator();
		while(specItr.hasNext())
		{
			Spectrum spec = specItr.next();
			System.out.println(spec.getScanNum()+"\t"+spec.getActivationMethod());
		}
	}
}
