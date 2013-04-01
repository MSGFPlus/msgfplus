package edu.ucsd.msjava.ims;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Iterator;

import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;

public class MaskSpectra {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit("Illegal parameters!");
		
		File specFile = new File(argv[0]);
		Iterator<Spectrum> itr = new SpectraAccessor(specFile).getSpecItr();
		if(itr == null)
			printUsageAndExit("Illegal spectrum format.");
		
		int windowSize = Integer.parseInt(argv[1]);
		
		String specPath = specFile.getPath();
		File outputFile = new File(specPath.substring(0, specPath.lastIndexOf('.'))+"_Mask"+windowSize+".mgf");
		mask(itr, outputFile, windowSize);
		System.out.println("Done");
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.err.println("Error: " + message);
		System.out.println("Usage: java MaskSpectra SpectrumFile WindowSize");
		System.exit(-1);
	}
	
	public static void mask(Iterator<Spectrum> itr, File outputFile, int windowSize) throws Exception
	{
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			float precurosrMz = spec.getPrecursorPeak().getMz();
			int minMz = Math.round(precurosrMz/windowSize)*windowSize - windowSize/2;
			int maxMz = minMz + windowSize/2;
			
			Spectrum newSpec = spec.getCloneWithoutPeakList();
			for(Peak p : spec)
			{
				if(p.getMz() < minMz || p.getMz() >= maxMz)
					newSpec.add(p);
			}
			newSpec.outputMgf(out);
		}
		out.close();
	}
}
