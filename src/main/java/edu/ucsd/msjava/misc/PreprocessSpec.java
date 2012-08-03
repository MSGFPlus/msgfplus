package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Iterator;


import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.SpecFileFormat;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MS2SpectrumParser;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;
import edu.ucsd.msjava.parser.PNNLSpectrumParser;
import edu.ucsd.msjava.parser.PklSpectrumParser;
import edu.ucsd.msjava.parser.SpectrumParser;

public class PreprocessSpec {
	public static void main(String argv[]) throws Exception
	{
		File source = null;
		File target = null;
		ActivationMethod activationMethod = null;
		InstrumentType instType = InstrumentType.LOW_RESOLUTION_LTQ;
		Enzyme	enzyme = Enzyme.TRYPSIN;
		boolean writeActivationMethod = false;
		
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length)
				printUsageAndExit("Illegal parameters");
			else if(argv[i].equalsIgnoreCase("-s"))
			{
				source = new File(argv[i+1]);
				if(!source.exists())
					printUsageAndExit(argv[i+1] + " doesn't exist!");
			}
			else if(argv[i].equalsIgnoreCase("-t"))
			{
				target = new File(argv[i+1]);
				if(!target.getName().endsWith(".mgf"))
					printUsageAndExit(argv[i+1] + " should end with .mgf!");
			}
			else if(argv[i].equalsIgnoreCase("-m"))	// Fragmentation method
			{
				// (0: written in the spectrum, 1: CID , 2: ETD, 3: HCD)
				if(argv[i+1].equalsIgnoreCase("0"))
				{
					activationMethod = null;
				}
				else if(argv[i+1].equalsIgnoreCase("1"))
				{
					activationMethod = ActivationMethod.CID;
				}
				else if(argv[i+1].equalsIgnoreCase("2"))
				{
					activationMethod = ActivationMethod.ETD;
				}
				else if(argv[i+1].equalsIgnoreCase("3"))
				{
					activationMethod = ActivationMethod.HCD;
				}
				else
					printUsageAndExit("Illegal activation method: " + argv[i+1]);
			}			
			else if(argv[i].equalsIgnoreCase("-m"))	// Fragmentation method
			{
				// (0: written in the spectrum, 1: CID , 2: ETD, 3: HCD)
				if(argv[i+1].equalsIgnoreCase("0"))
				{
					activationMethod = null;
				}
				else if(argv[i+1].equalsIgnoreCase("1"))
				{
					activationMethod = ActivationMethod.CID;
				}
				else if(argv[i+1].equalsIgnoreCase("2"))
				{
					activationMethod = ActivationMethod.ETD;
				}
				else if(argv[i+1].equalsIgnoreCase("3"))
				{
					activationMethod = ActivationMethod.HCD;
				}
				else
					printUsageAndExit("Illegal activation method: " + argv[i+1]);
			}			
			else if(argv[i].equalsIgnoreCase("-inst"))	// Instrument type
			{
				if(argv[i+1].equalsIgnoreCase("0"))
				{
					instType = InstrumentType.LOW_RESOLUTION_LTQ;
				}
				else if(argv[i+1].equalsIgnoreCase("1"))
				{
					instType = InstrumentType.TOF;
				}
				else if(argv[i+1].equalsIgnoreCase("2"))
				{
					instType = InstrumentType.HIGH_RESOLUTION_LTQ;
				}
				else
				{
					printUsageAndExit("Illegal instrument type: " + argv[i+1]);
				}
			}			
			else if(argv[i].equalsIgnoreCase("-e"))	// Enzyme
			{
				// 0: No enzyme, 1: Trypsin, 2: Chymotrypsin, 3: LysC, 4: LysN, 5: GluC, 6: ArgC, 7: AspN
				if(argv[i+1].equalsIgnoreCase("0"))
					enzyme = null;
				else if(argv[i+1].equalsIgnoreCase("1"))
					enzyme = Enzyme.TRYPSIN;
				else if(argv[i+1].equalsIgnoreCase("2"))
					enzyme = Enzyme.CHYMOTRYPSIN;
				else if(argv[i+1].equalsIgnoreCase("3"))
					enzyme = Enzyme.LysC;
				else if(argv[i+1].equalsIgnoreCase("4"))
					enzyme = Enzyme.LysN;
				else if(argv[i+1].equalsIgnoreCase("5"))
					enzyme = Enzyme.GluC;
				else if(argv[i+1].equalsIgnoreCase("6"))
					enzyme = Enzyme.ArgC;
				else if(argv[i+1].equalsIgnoreCase("7"))
					enzyme = Enzyme.AspN;
				else if(argv[i+1].equalsIgnoreCase("8"))
					enzyme = Enzyme.ALP;
				else if(argv[i+1].equalsIgnoreCase("9"))
					enzyme = Enzyme.Peptidomics;
				else
					printUsageAndExit("Illegal enzyme: " + argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-w"))
			{
				if(argv[i+1].equals("0"))
					writeActivationMethod = false;
				else if(argv[i+1].equals("1"))
					writeActivationMethod = true;
			}
			else
			{
				printUsageAndExit("Invalid parameter: " + argv[i]);
			}
		}

		if(source == null || target == null)
			printUsageAndExit("Invalid parameters!");
		convert(source, target, activationMethod, instType, enzyme, writeActivationMethod);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println(message);
		System.out.println("Usage: java ConvertToMgf\n" +
				"\t-s SourceFile or Directory\n" +
				"\t-t TargetFileName (*.mgf)\n" +
				"\t[-m FragmentationMethodID] (0: as written in the spectrum or CID if no info (Default), 1: CID, 2: ETD, 3: HCD\n" +
				"\t[-inst InstrumentID] (0: Low-res LCQ/LTQ (Default for CID and ETD), 1: TOF , 2: High-res LTQ (Default for HCD))\n" +
				"\t[-e EnzymeID] (0: No enzyme, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: endogenous peptides)\n" +
				"\t[-w 0/1] (0: don't write ActivationMethod (default), 1: write ActivationMethod)\n"
				);
		
		System.exit(-1);
	}
	
	public static void convert(File source, File target, ActivationMethod activationMethod, InstrumentType instType, Enzyme enzyme, boolean writeActivationMethod) throws Exception
	{
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(target)));
		
		File[] fileList;
		if(!source.isDirectory())
		{
			fileList = new File[1];
			fileList[0] = source;
		}
		else
			fileList = source.listFiles();
			
		int numFileConverted = 0;
		for(File sourceFile : fileList)
		{
			String specFileName = sourceFile.getName();
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
			if(sourceFileFormat != null)
			{
				convertFile(sourceFile, sourceFileFormat, target, activationMethod, instType, enzyme, out, writeActivationMethod);
				numFileConverted++;
			}
		}
		out.close();
		System.out.println("Converted " + numFileConverted + " files.");
	}
	
	public static void convertFile(File sourceFile, SpecFileFormat sourceFileFormat, File target,  ActivationMethod activationMethod, InstrumentType instType, Enzyme enzyme, PrintStream out, boolean writeActivationMethod) throws Exception
	{
		System.out.print(sourceFile.getName() + ": ");
    	Iterator<Spectrum> specItr = null;
		
		if(sourceFileFormat == SpecFileFormat.MZXML)
		{
			specItr = new MzXMLSpectraIterator(sourceFile.getPath());
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
				specItr = new SpectraIterator(sourceFile.getPath(), parser);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		NewRankScorer scorer = null;
		
		if(activationMethod != null)
			scorer = NewScorerFactory.get(activationMethod, instType, enzyme, null);
		
		int numSpecs = 0;
		while(specItr.hasNext())
		{
			Spectrum spec = specItr.next();
			if(activationMethod == null || activationMethod == ActivationMethod.FUSION)
				scorer = NewScorerFactory.get(spec.getActivationMethod(), instType, enzyme, null);

			scorer.filterPrecursorPeaks(spec);
			if(activationMethod != null && spec.getActivationMethod() != null && spec.getActivationMethod() != activationMethod)
				continue;
			
			spec.outputMgf(out, writeActivationMethod);
			numSpecs++;
		}
		System.out.println(numSpecs + " spectra preprocessed.");
	}
}
