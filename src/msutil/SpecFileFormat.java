package msutil;

import java.io.File;
import java.util.ArrayList;

import jmzparser.MzMLSpectraMap;

import parser.MS2SpectrumParser;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraMap;
import parser.PNNLSpectraMap;
import parser.PklSpectrumParser;
import parser.SpectrumParser;

public class SpecFileFormat extends FileFormat {
	private SpecFileFormat(String[] suffixes)
	{
		super(suffixes);
	}
	
	private SpecFileFormat(String suffix)
	{
		super(suffix);
	}
	
	public static final SpecFileFormat MGF;
	public static final SpecFileFormat MZXML;
	public static final SpecFileFormat MZML;
	public static final SpecFileFormat MS2;
	public static final SpecFileFormat PKL;
	public static final SpecFileFormat MZDATA;
	public static final SpecFileFormat DTA_TXT;
	
	public static SpectrumAccessorBySpecIndex getSpecMap(File specFile)
	{
		SpectrumAccessorBySpecIndex specMap = null;
		SpecFileFormat specFormat = getSpecFileFormat(specFile.getName());
		if(specFormat == null)
			return null;
		
		if(specFormat == MZXML)
			specMap = new MzXMLSpectraMap(specFile.getPath());
		else if(specFormat == MZML)
			specMap = new MzMLSpectraMap(specFile.getPath());
		else if(specFormat == SpecFileFormat.DTA_TXT)
			specMap = new PNNLSpectraMap(specFile.getPath());
		else
		{
			SpectrumParser parser = null;
			if(specFormat == SpecFileFormat.MGF)
				parser = new MgfSpectrumParser();
			else if(specFormat == SpecFileFormat.MS2)
				parser = new MS2SpectrumParser();
			else if(specFormat == SpecFileFormat.PKL)
				parser = new PklSpectrumParser();
			specMap = new SpectraMap(specFile.getPath(), parser);
		}
		
		return specMap;
	}
	
	public static SpecFileFormat getSpecFileFormat(String specFileName)
	{
		String lowerCaseFileName = specFileName.toLowerCase();
		for(SpecFileFormat f : specFileFormatList)
		{
			for(String suffix : f.getSuffixes())
			{
				if(lowerCaseFileName.endsWith(suffix.toLowerCase()))
					return f;
			}
		}
		return null;
	}
	
	private static ArrayList<SpecFileFormat> specFileFormatList;
	static {
		MGF = new SpecFileFormat(".mgf");
		MZXML = new SpecFileFormat(".mzXML");
		MZML = new SpecFileFormat(".mzML");
		MS2 = new SpecFileFormat(".ms2");
		PKL = new SpecFileFormat(".pkl");
		MZDATA = new SpecFileFormat(".mzData");
		DTA_TXT = new SpecFileFormat("_dta.txt");
		
		specFileFormatList = new ArrayList<SpecFileFormat>();
		specFileFormatList.add(MGF);
		specFileFormatList.add(MZXML);
		specFileFormatList.add(MZML);
		specFileFormatList.add(MS2);
		specFileFormatList.add(PKL);
		specFileFormatList.add(MZDATA);
		specFileFormatList.add(DTA_TXT);
	}
}
