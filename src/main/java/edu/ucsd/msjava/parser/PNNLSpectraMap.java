package edu.ucsd.msjava.parser;

import java.util.HashMap;
import java.util.Iterator;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.SpectraMap;
import edu.ucsd.msjava.msutil.Spectrum;

public class PNNLSpectraMap extends SpectraMap {

	private HashMap<Integer,ActivationMethod> scanNumActMethodMap;
	
	public PNNLSpectraMap(String fileName) 
	{
		super(fileName, new PNNLSpectrumParser());
		scanNumActMethodMap = PNNLSpectrumParser.getScanTypeMap(fileName);	
	}

	@Override
	public synchronized Spectrum getSpectrumBySpecIndex(int specIndex)
	{
		if(scanNumActMethodMap == null)
			return super.getSpectrumBySpecIndex(specIndex);
		else
		{
			Spectrum spec = super.getSpectrumBySpecIndex(specIndex);
			ActivationMethod method = scanNumActMethodMap.get(spec.getScanNum());
			if(method != null)
				spec.setActivationMethod(method);
			return spec;
		}
	}
	
	public static void main(String argv[]) throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Test/Matt/QC_Shew_11_03_200ng_4_23Aug11_Hawk_11-05-04p_dta.txt";
		PNNLSpectraMap map = new PNNLSpectraMap(fileName);
		for(int specIndex : map.getSpecIndexList())
		{
			Spectrum spec = map.getSpectrumBySpecIndex(specIndex);
			System.out.println(spec.getScanNum()+"\t"+spec.getActivationMethod());
		}
	}
	
	
}
