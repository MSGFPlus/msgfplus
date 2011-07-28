package misc;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Iterator;

import parser.MS2SpectrumParser;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraIterator;
import parser.PNNLSpectrumParser;
import parser.PklSpectrumParser;
import parser.SpectrumParser;

import msutil.SpecFileFormat;
import msutil.SpectraIterator;
import msutil.Spectrum;

public class ConvertToMgf {
	public static void main(String argv[]) throws Exception
	{
		boolean writeActivationMethod = false;
		if(argv.length != 2 && argv.length != 3)
			printUsageAndExit("Wrong parameters!");
		if(argv.length == 3 && argv[2].equalsIgnoreCase("1"))
			writeActivationMethod = true;
		convert(argv[0], argv[1], writeActivationMethod);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println(message);
		System.out.println("Usage: java misc.SpecFileConverter SourceFileName TargetFileName(*.mgf) [0/1] (0: don't write ActivationMethod, 1: write ActivationMethod) ");
		System.exit(-1);
	}
	
	public static void convert(String source, String target, boolean writeActivationMethod) throws Exception
	{
		String specFileName = source;
		SpecFileFormat sourceFileFormat = null;
		int posDot = specFileName.lastIndexOf('.');
		if(posDot >= 0)
		{
			String extension = specFileName.substring(posDot);
			if(extension.equalsIgnoreCase(".mzXML"))
				sourceFileFormat = SpecFileFormat.MZXML;
			else if(extension.equalsIgnoreCase(".mgf"))
				sourceFileFormat = SpecFileFormat.MGF;
			else if(extension.equalsIgnoreCase(".ms2"))
				sourceFileFormat = SpecFileFormat.MS2;
			else if(extension.equalsIgnoreCase(".pkl"))
				sourceFileFormat = SpecFileFormat.PKL;
		}		
		if(sourceFileFormat == null && specFileName.length() > 8)
		{
			String suffix = specFileName.substring(specFileName.length()-8);
			if(suffix.equalsIgnoreCase("_dta.txt"))
				sourceFileFormat = SpecFileFormat.DTA_TXT;
		}
		
		if(sourceFileFormat == null)
			printUsageAndExit("Unsupported file format: " + source);
		
		if(!target.endsWith(".mgf"))
			printUsageAndExit("Target file must end with .mgf: " + target);
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(target)));
		
    	Iterator<Spectrum> specItr = null;
		
		if(sourceFileFormat == SpecFileFormat.MZXML)
		{
			specItr = new MzXMLSpectraIterator(specFileName);
		}
		else
		{
			SpectrumParser parser = null;
			if(sourceFileFormat == SpecFileFormat.MGF)
				parser = new MgfSpectrumParser();
			else if(sourceFileFormat == SpecFileFormat.DTA_TXT)
				parser = new PNNLSpectrumParser();
			else if(sourceFileFormat == SpecFileFormat.MS2)
				parser = new MS2SpectrumParser();
			else if(sourceFileFormat == SpecFileFormat.PKL)
				parser = new PklSpectrumParser();
			
			try {
				specItr = new SpectraIterator(specFileName, parser);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		int numSpecs = 0;
		while(specItr.hasNext())
		{
			Spectrum spec = specItr.next();
			spec.outputMgf(out, writeActivationMethod);
			numSpecs++;
		}
		out.close();
		
		System.out.println(numSpecs + " spectra converted.");
	}
}
