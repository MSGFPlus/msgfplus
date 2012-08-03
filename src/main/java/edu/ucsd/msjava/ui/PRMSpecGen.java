package edu.ucsd.msjava.ui;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Iterator;


import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScoredSpectrum;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Pair;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.SpecFileFormat;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.SpectraMap;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumAccessorBySpecIndex;
import edu.ucsd.msjava.parser.MS2SpectrumParser;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;
import edu.ucsd.msjava.parser.MzXMLSpectraMap;
import edu.ucsd.msjava.parser.PNNLSpectrumParser;
import edu.ucsd.msjava.parser.PklSpectrumParser;
import edu.ucsd.msjava.parser.SpectrumParser;
import edu.ucsd.msjava.sequences.Constants;


public class PRMSpecGen {
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();
		if(argv.length < 2 || argv.length % 2 != 0)
			printUsageAndExit("The number of parameters must be even.");

		File 	specFile = null;
		SpecFileFormat specFormat = null;
		File	outputFile = null;
//		Tolerance fragmentMassTolerance = null;
		Enzyme	enzyme = Enzyme.TRYPSIN;
		ActivationMethod activationMethod = null;
		InstrumentType instType = InstrumentType.LOW_RESOLUTION_LTQ;
		
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length)
				printUsageAndExit("Illegal parameters");
			else if(argv[i].equalsIgnoreCase("-s"))
			{
				specFile = new File(argv[i+1]);
				if(!specFile.exists())
				{
					printUsageAndExit(argv[i+1]+" doesn't exist.");
				}
				if(specFile.isDirectory())
				{
					printUsageAndExit(argv[i+1]+" must not be a directory!");
				}
				else
				{
					String specFileName = specFile.getName();
					int posDot = specFileName.lastIndexOf('.');
					if(posDot >= 0)
					{
						String extension = specFileName.substring(posDot);
						if(extension.equalsIgnoreCase(".mzXML"))
							specFormat = SpecFileFormat.MZXML;
						else if(extension.equalsIgnoreCase(".mgf"))
							specFormat = SpecFileFormat.MGF;
						else if(extension.equalsIgnoreCase(".ms2"))
							specFormat = SpecFileFormat.MS2;
						else if(extension.equalsIgnoreCase(".pkl"))
							specFormat = SpecFileFormat.PKL;
					}		
					if(specFormat == null && specFileName.length() > 8)
					{
						String suffix = specFileName.substring(specFileName.length()-8);
						if(suffix.equalsIgnoreCase("_dta.txt"))
							specFormat = SpecFileFormat.DTA_TXT;
					}
				}
				if(specFormat == null)
					printUsageAndExit("Illegal spectrum format: " + argv[i+1]);
			}
//			else if(argv[i].equalsIgnoreCase("-f"))
//			{
//				fragmentMassTolerance = Tolerance.parseToleranceStr(argv[i+1]);
//				if(fragmentMassTolerance == null)
//					printUsageAndExit("Illegal fragment mass tolerance value: " + argv[i+1]);
//			}
//			else if(argv[i].equalsIgnoreCase("-i"))
//			{
//				String[] token = argv[i+1].split(",");
//				ionTypes = new IonType[token.length];
//				for(int j=0; j<token.length; j++)
//				{
//					ionTypes[j] = IonType.getIonType(token[j]);
//					if(ionTypes[j] == null)
//					{
//						printUsageAndExit("Unrecognizable ion type: " + token[i]);
//					}
//				}
//				for(IonType ion : ionTypes)
//					System.out.println(ion.getName()+"\t"+ion.toString());
//				System.exit(-1);
//			}			
			else if(argv[i].equalsIgnoreCase("-o"))
			{
				outputFile = new File(argv[i+1]);
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
				else
					printUsageAndExit("Illegal enzyme: " + argv[i+1]);
			}
			else
			{
				printUsageAndExit("Invalid option: " + argv[i]);
			}
		}
		
		if(specFile == null)
			printUsageAndExit("Spectrum is not specified.");
		
		if(outputFile == null)
			printUsageAndExit("Output file is not specified.");
		
		if(activationMethod == ActivationMethod.HCD)
			instType = InstrumentType.HIGH_RESOLUTION_LTQ;
		
		generatePRMSpectrum(specFile, specFormat, 
	    		outputFile, enzyme, activationMethod, instType);
		System.out.println("Complete.");
		System.out.format("Time: %.3f sec\n", (System.currentTimeMillis()-time)/(float)1000);
	}
	
	public static void printUsageAndExit()
	{
		printUsageAndExit(null);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println("Error: " + message + "\n");
		System.out.println("PRMSpecGen v"+ MSGFDB.VERSION + " (" + MSGFDB.RELEASE_DATE + ")");
		System.out.print("Usage: java -Xmx500M -cp MSGFDB.jar ui.PRMSpecGen\n"
				+ "\t-s SpectrumFile (*.mzXML, *.mgf, *.ms2, *.pkl or *_dta.txt)\n" //, *.mgf, *.pkl, *.ms2)\n"
//				+ "\t-f FragMassTolerance (fragment mass tolerance in ppm or Da. The value must be less than 0.5Da or 100ppm. E.g. 0.5Da or 30ppm)\n"
				+ "\t-o outputFileName (e.g. PRMSpec.mgf)\n"
				+ "\t[-m FragmentationMethodID] (0: as written in the spectrum or CID if no info (Default), 1: CID, 2: ETD, 3: HCD)\n"
				+ "\t[-inst InstrumentID] (0: Low-res LCQ/LTQ (Default for CID and ETD), 1: TOF , 2: High-res LTQ (Default for HCD))\n"
				+ "\t[-e EnzymeID] (0: No enzyme, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N)\n"
				);
		System.exit(-1);
	}
	
    public static void generatePRMSpectrum(
    		File specFile, 
    		SpecFileFormat specFormat, 
    		File outputFile, 
    		Enzyme enzyme, 
    		ActivationMethod activationMethod,  
    		InstrumentType instType
    		)
	{
    	// Output
		PrintStream out = null;
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
    	
		SpectrumAccessorBySpecIndex specMap = null;
    	Iterator<Spectrum> specItr = null;
		
		if(specFormat == SpecFileFormat.MZXML)
		{
			specItr = new MzXMLSpectraIterator(specFile.getPath());
			specMap = new MzXMLSpectraMap(specFile.getPath());
		}
		else
		{
			SpectrumParser parser = null;
			if(specFormat == SpecFileFormat.MGF)
				parser = new MgfSpectrumParser();
			else if(specFormat == SpecFileFormat.DTA_TXT)
				parser = new PNNLSpectrumParser();
			else if(specFormat == SpecFileFormat.MS2)
				parser = new MS2SpectrumParser();
			else if(specFormat == SpecFileFormat.PKL)
				parser = new PklSpectrumParser();
			
			try {
				specItr = new SpectraIterator(specFile.getPath(), parser);
				specMap = new SpectraMap(specFile.getPath(), parser);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		if(specItr == null || specMap == null)
		{
			printUsageAndExit("Error while parsing spectrum file: " + specFile.getPath());
		}

		int totalNumSpecs = specMap.getSpecIndexList().size();
		
		NewRankScorer scorer = null;
		if(activationMethod != null)
			scorer = NewScorerFactory.get(activationMethod, instType, enzyme, null);

		System.out.println("Total number of spectra: " + totalNumSpecs);
		int numSpecs = 0;
		while(specItr.hasNext())
		{
			Spectrum spec = specItr.next();
			numSpecs++;
			if(numSpecs % 1000 == 0)
			{
				System.out.format("Processing spectra... %.4f", (numSpecs*100/(float)totalNumSpecs));
				System.out.println("% done.");
			}
			if(spec.size() < Constants.MIN_NUM_PEAKS_PER_SPECTRUM)
			{
				System.out.println("Spectrum " + spec.getSpecIndex() + " has too few peaks (#Peaks: " + spec.size()+"): ignored.");
				continue;
			}
			if(spec.getCharge() <= 0)
			{
				System.out.println("Spectrum " + spec.getSpecIndex() + " has zero or negative charge: ignored.");
				continue;
			}
			
			if(activationMethod == null || activationMethod == ActivationMethod.FUSION)
				scorer = NewScorerFactory.get(spec.getActivationMethod(), instType, enzyme, Protocol.NOPROTOCOL);
			
			scorer.doNotUseError();
			NewScoredSpectrum<NominalMass> scoredSpec = scorer.getScoredSpectrum(spec);
			int maxNominalMass = NominalMass.toNominalMass(spec.getParentMass());
			
			// PRM spectrum
			out.println("BEGIN IONS");
			out.print("TITLE=PRM_SpecIndex="+spec.getSpecIndex());
		    if(spec.getTitle() != null)
		        out.println(" " + spec.getTitle());
		    else
		    	out.println();
		    if(spec.getAnnotation() != null)
		    	out.println("SEQ=" + spec.getAnnotationStr());
			out.println("PEPMASS=" + spec.getPrecursorPeak().getMz());
			out.println("SCANS=" + spec.getScanNum());
			out.println("CHARGE="+spec.getCharge()+"+");
			for(int m=1; m<maxNominalMass; m++)
			{
				Pair<Float,Float> massScorePair = scoredSpec.getNodeMassAndScore(NominalMass.getMassFromNominalMass(m), true);
				Float peakDerivedMass = massScorePair.getFirst();
				float score = massScorePair.getSecond();
				if(peakDerivedMass == null)
					out.format("%d", m);
				else
					out.format("%f", peakDerivedMass);
				out.format("\t%.3f\n",score);
			}
			out.println("END IONS");

			// SRM spectrum
			out.println("BEGIN IONS");
			out.print("TITLE=SRM_SpecIndex="+spec.getSpecIndex());
		    if(spec.getTitle() != null)
		        out.println(" " + spec.getTitle());
		    else
		    	out.println();
		    if(spec.getAnnotation() != null)
		    	out.println("SEQ=" + spec.getAnnotationStr());
			out.println("PEPMASS=" + spec.getPrecursorPeak().getMz());
			out.println("SCANS=" + spec.getScanNum());
			out.println("CHARGE="+spec.getCharge()+"+");
			for(int m=1; m<maxNominalMass; m++)
			{
				Pair<Float,Float> massScorePair = scoredSpec.getNodeMassAndScore(NominalMass.getMassFromNominalMass(m), false);
				Float peakDerivedMass = massScorePair.getFirst();
				float score = massScorePair.getSecond();
				if(peakDerivedMass == null)
					out.format("%d", m);
				else
					out.format("%f", peakDerivedMass);
				out.format("\t%.3f\n",score);
			}
			out.println("END IONS");
    	}
		out.close();
	}	
}
