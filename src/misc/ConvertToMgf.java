package misc;

import java.io.BufferedOutputStream;
import java.io.File;
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

import msutil.ActivationMethod;
import msutil.SpecFileFormat;
import msutil.SpectraIterator;
import msutil.Spectrum;

public class ConvertToMgf {
	public static void main(String argv[]) throws Exception
	{
		boolean writeActivationMethod = false;
		ActivationMethod activationMethod = null;
		File source = null;
		File target = null;
		int specIndex = -1;
		
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
			else if(argv[i].equalsIgnoreCase("-w"))
			{
				if(argv[i+1].equals("0"))
					writeActivationMethod = false;
				else if(argv[i+1].equals("1"))
					writeActivationMethod = true;
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
			else if(argv[i].equalsIgnoreCase("-index"))
			{
				specIndex = Integer.parseInt(argv[i+1]);
			}
			else
			{
				printUsageAndExit("Invalid parameter: " + argv[i]);
			}
		}

		if(source == null || target == null)
			printUsageAndExit("Invalid parameters!");
		convert(source, target, writeActivationMethod, activationMethod, specIndex);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println(message);
		System.out.println("Usage: java ConvertToMgf\n" +
				"\t-s SourceFile or Directory\n" +
				"\t-t TargetFileName (*.mgf)\n" +
				"\t[-w 0/1] (0: don't write ActivationMethod (default), 1: write ActivationMethod)\n" +
				"\t[-m FragmentationMethodID] (0: convert all (Default), 1: CID, 2: ETD, 3: HCD)\n" +
				"\t[-index startIndex(,endIndex)] (only write the spectrum of indices with the specified range, startIndex: inclu)");
		System.exit(-1);
	}
	
	public static void convert(File source, File target, boolean writeActivationMethod, ActivationMethod activationMethod, int specIndex) throws Exception
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
				convertFile(sourceFile, sourceFileFormat, target, writeActivationMethod, activationMethod, specIndex, out);
				numFileConverted++;
			}
		}
		out.close();
		System.out.println("Converted " + numFileConverted + " files.");
	}
	
	public static void convertFile(File sourceFile, SpecFileFormat sourceFileFormat, File target, boolean writeActivationMethod, ActivationMethod activationMethod, int specIndex, PrintStream out) throws Exception
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
		
		int numSpecs = 0;
		while(specItr.hasNext())
		{
			Spectrum spec = specItr.next();
			if(specIndex > 0 && spec.getSpecIndex() != specIndex)
				continue;
			if(activationMethod != null && spec.getActivationMethod() != activationMethod)
				continue;
			spec.outputMgf(out, writeActivationMethod);
			numSpecs++;
		}
		System.out.println(numSpecs + " spectra converted.");
	}
}
