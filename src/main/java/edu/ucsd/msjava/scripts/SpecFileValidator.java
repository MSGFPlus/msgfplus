package edu.ucsd.msjava.scripts;

import java.io.File;
import java.util.Iterator;

import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MS2SpectrumParser;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;
import edu.ucsd.msjava.parser.PklSpectrumParser;

public class SpecFileValidator {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1)
			printUsageAndExit(null);
		
		File specFile = new File(argv[0]);
		validateSpecFile(specFile);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println("Error: " + message);
		System.out.println("usage: java -jar SpecFileValidator.jar SpecFileValidator fileName");
		System.exit(-1);
	}
	
	public static void validateSpecFile(File specFile)
	{
		if(!specFile.exists())
		{
			System.out.println(specFile.getPath() + " is invalid! (file not exists)");
			System.exit(-1);
		}
		
		String fileName = specFile.getName();
		String ext = fileName.substring(fileName.lastIndexOf('.')+1);
		
		int numSpectra = 0;
		
		try {
			Iterator<Spectrum> itr = null;
			if(ext.equalsIgnoreCase("mzXML"))
			{
				itr = new MzXMLSpectraIterator(specFile.getPath());
			}
			else if(ext.equalsIgnoreCase("mgf"))
			{
				itr = new SpectraIterator(specFile.getPath(), new MgfSpectrumParser());
			}
			else if(ext.equalsIgnoreCase("pkl"))
			{
				itr = new SpectraIterator(specFile.getPath(), new PklSpectrumParser());
			}
			else if(ext.equalsIgnoreCase("ms2"))
			{
				itr = new SpectraIterator(specFile.getPath(), new MS2SpectrumParser());
			}
			
			if(itr == null)
			{
				System.out.println(specFile.getPath() + " is invalid! (not a valid file format)");
				System.exit(-1);
			}
			
			while(itr.hasNext())
			{
				itr.next();
				numSpectra++;
			}
		}
		catch (Exception e)
		{
			System.out.println(specFile.getPath() + " is invalid! (parsing error)");
			System.exit(-1);
		}

		if(numSpectra == 0)
		{
			System.out.println(specFile.getPath() + " is invalid!");
			System.exit(-1);
		}
		
		System.out.println(specFile.getPath() + " is validated.");
		System.out.println("NumberOfSpectra: " + numSpectra);
	}
}
