package ui;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import parser.MS2SpectrumParser;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraIterator;
import parser.MzXMLSpectraMap;
import parser.PNNLSpectraIterator;
import parser.PNNLSpectraMap;
import parser.PklSpectrumParser;
import parser.SpectrumParser;

import msdbsearch.CompactFastaSequence;
import msdbsearch.CompactSuffixArray;
import msdbsearch.ConcurrentMSGFDB;
import msdbsearch.DBScanner;
import msdbsearch.ReverseDB;
import msdbsearch.ScoredSpectraMap;

import msgf.MSGFDBResultGenerator;
import msgf.Tolerance;
import msscorer.NewScorerFactory.SpecDataType;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.ActivationMethod;
import msutil.InstrumentType;
import msutil.Modification;
import msutil.SpecFileFormat;
import msutil.SpecKey;
import msutil.SpectraIterator;
import msutil.SpectraMap;
import msutil.Spectrum;
import msutil.SpectrumAccessorBySpecIndex;

public class MSGFDB {
	public static final String VERSION = "6955";
	public static final String RELEASE_DATE = "12/11/2011";
	
	public static final String DECOY_DB_EXTENSION = ".revConcat.fasta";
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();
		if(argv.length < 2 || argv.length % 2 != 0)
			printUsageAndExit("The number of parameters must be even.");

		File 	specFile = null;
		SpecFileFormat specFormat = null;
		File 	databaseFile 	= null;
		File	paramFile	= null;
		File 	dbIndexDir 	= null;
		File	outputFile = null;
//		Tolerance parentMassTolerance = null;
		Tolerance leftParentMassTolerance = null;
		Tolerance rightParentMassTolerance = null;
		Boolean isTolerancePPM = null;
		int numAllowedC13 = 1;
		int 	numMatchesPerSpec = 1;
		Enzyme	enzyme = Enzyme.TRYPSIN;
		ActivationMethod activationMethod = null;
		
		InstrumentType instType = InstrumentType.LOW_RESOLUTION_LTQ;
		int numAllowedNonEnzymaticTermini = 1;
//		boolean showTitle = false;
		boolean useTDA = false;
		boolean showFDR = true;
		boolean replicateMergedResults = false;
		boolean doNotDseEdgeScore = false;
		boolean useUniformAAProb = false;
		int minPeptideLength = 6;
		int maxPeptideLength = 40;
		int minCharge = 2;
		int maxCharge = 3;
		int startSpecIndex = 0;
		int endSpecIndex = Integer.MAX_VALUE;
		int numThreads = Runtime.getRuntime().availableProcessors();
		
		AminoAcidSet aaSet = null;
		
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
						if(extension.equalsIgnoreCase(".mzXML") || extension.equalsIgnoreCase(".mzML"))
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
			else if(argv[i].equalsIgnoreCase("-d"))
			{
				databaseFile = new File(argv[i+1]);
				if(!databaseFile.exists())
				{
					printUsageAndExit(argv[i+1]+" doesn't exist.");
				}
				String databaseFileName = databaseFile.getName();
				String dbExt = databaseFileName.substring(databaseFileName.lastIndexOf('.')+1);
				if(!dbExt.equalsIgnoreCase("fasta") && !dbExt.equalsIgnoreCase("fa"))
				{
					printUsageAndExit("Database file must ends with .fa or .fasta!: " + argv[i+1]);
				}
			}
			else if(argv[i].equalsIgnoreCase("-dd"))
			{
				dbIndexDir = new File(argv[i+1]);
				if(!dbIndexDir.exists())
				{
					printUsageAndExit(argv[i+1]+" doesn't exist.");
				}
				if(!dbIndexDir.isDirectory())
				{
					printUsageAndExit(argv[i+1]+" is not a directory.");
				}
			}
			else if(argv[i].equalsIgnoreCase("-param"))
			{
				paramFile = new File(argv[i+1]);
				if(!paramFile.exists())
				{
					printUsageAndExit(argv[i+1]+" doesn't exist.");
				}
			}
			else if(argv[i].equalsIgnoreCase("-t"))
			{
				String[] token = argv[i+1].split(",");
				if(token.length == 1)
				{
					leftParentMassTolerance = rightParentMassTolerance = Tolerance.parseToleranceStr(token[0]);
				}
				else if(token.length == 2)
				{
					leftParentMassTolerance = Tolerance.parseToleranceStr(token[0]);
					rightParentMassTolerance = Tolerance.parseToleranceStr(token[1]);
				}
				if(leftParentMassTolerance == null || rightParentMassTolerance == null)
				{
					printUsageAndExit("Illegal tolerance value: " + argv[i+1]);
				}
				if(leftParentMassTolerance.isTolerancePPM() != rightParentMassTolerance.isTolerancePPM())
				{
					printUsageAndExit("Left and right tolerance units must be the same: " + argv[i+1]);
				}
				if(leftParentMassTolerance.getValue() < 0 || rightParentMassTolerance.getValue() < 0)
				{
					printUsageAndExit("Parent mass tolerance must not be negative: " + argv[i+1]);
				}
			}
			else if(argv[i].equalsIgnoreCase("-u"))	// hidden option for ccms workflow
			{
				if(argv[i+1].equalsIgnoreCase("1"))
					isTolerancePPM = true;
				else if(argv[i+1].equalsIgnoreCase("0"))
					isTolerancePPM = false;
			}
			else if(argv[i].equalsIgnoreCase("-c13"))
			{
				try {
					numAllowedC13 = Integer.parseInt(argv[i+1]);
					if(numAllowedC13 != 0 && numAllowedC13 != 1 && numAllowedC13 != 2)
					{
						printUsageAndExit("Illegal -c13 value (must be 0/1/2): " + argv[i+1]);
					}
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal numMatchesPerSpec: " + argv[i+1]);
				} 
				
			}
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
				else if(argv[i+1].equalsIgnoreCase("4"))
				{
					activationMethod = ActivationMethod.FUSION;
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
			else if(argv[i].equalsIgnoreCase("-mod"))
			{
				File modFile = new File(argv[i+1]);
				if(!modFile.exists())
				{
					printUsageAndExit(modFile + " doesn't exist.");
				}
				String modFileName = modFile.getName();
				if(!modFileName.substring(modFileName.lastIndexOf('.')+1).equalsIgnoreCase("xml"))	// custom mod file
				{
					aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
				}
				else	// mod file for ccms workflow
				{
					aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(modFile.getPath());
				}
			}
			else if(argv[i].equalsIgnoreCase("-n"))
			{
				try {
					numMatchesPerSpec = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal numMatchesPerSpec: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-scan"))
			{
				String[] token = argv[i+1].split(",");
				if(token.length != 1 && token.length != 2)
				{
					printUsageAndExit("Illigal scanNum: " + argv[i+1]);
				}
				try {
					startSpecIndex = Integer.parseInt(token[0]);
					if(token.length == 2)
						endSpecIndex = Integer.parseInt(token[1]);
					else
						endSpecIndex = startSpecIndex+1;
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal scanNum: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-thread"))
			{
				try {
					numThreads = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal number of threads: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-minLength"))
			{
				try {
					minPeptideLength = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal minLength: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-maxLength"))
			{
				try {
					maxPeptideLength = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal maxLength: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-minCharge"))
			{
				try {
					minCharge = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal minCharge: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-maxCharge"))
			{
				try {
					maxCharge = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal maxCharge: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-nnet"))
			{
				try {
					numAllowedNonEnzymaticTermini = Integer.parseInt(argv[i+1]);
					if(numAllowedNonEnzymaticTermini != 0 && numAllowedNonEnzymaticTermini != 1 && numAllowedNonEnzymaticTermini !=2)
						printUsageAndExit("Illegal -nnet (0/1/2 are allowed): " + argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal numMatchesPerSpec: " + argv[i+1]);
				} 
			}
//			else if(argv[i].equalsIgnoreCase("-title"))
//			{
//				if(argv[i+1].equalsIgnoreCase("1"))
//					showTitle = true;
//			}
			else if(argv[i].equalsIgnoreCase("-uniformAAProb"))
			{
				if(argv[i+1].equalsIgnoreCase("1"))
					useUniformAAProb = true;
			}
			else if(argv[i].equalsIgnoreCase("-tda"))
			{
				if(argv[i+1].equalsIgnoreCase("1"))
					useTDA = true;
				else if(argv[i+1].equalsIgnoreCase("0"))
					useTDA = false;
				else
				{
					printUsageAndExit("Illigal -tda parameter: " + argv[i+1]);
				}
			}
			else if(argv[i].equalsIgnoreCase("-showFDR"))
			{
				if(argv[i+1].equalsIgnoreCase("1"))
					showFDR = true;
				else if(argv[i+1].equalsIgnoreCase("0"))
					showFDR = false;
				else
				{
					printUsageAndExit("Illigal -showFDR parameter: " + argv[i+1]);
				}
			}
			else if(argv[i].equalsIgnoreCase("-replicate"))
			{
				if(argv[i+1].equalsIgnoreCase("1"))
					replicateMergedResults = true;
				else if(argv[i+1].equalsIgnoreCase("0"))
					replicateMergedResults = false;
				else
				{
					printUsageAndExit("Illigal -replicate parameter: " + argv[i+1]);
				}
			}
			else
			{
				printUsageAndExit("Invalid parameter: " + argv[i]);
			}
		}
		
		if(specFile == null)
			printUsageAndExit("Spectrum is not specified.");
		if(databaseFile == null)
			printUsageAndExit("Database is not specified.");
		
		if(leftParentMassTolerance == null || rightParentMassTolerance == null)
			printUsageAndExit("Parent mass tolerance is not specified.");

		if(minPeptideLength > maxPeptideLength)
			printUsageAndExit("MinPepLength must not be larger than MaxPepLength!");
			
		if(minCharge > maxCharge)
			printUsageAndExit("MinPrecursorCharge must not be larger than MaxPrecursorCharge!");
		
		if(isTolerancePPM != null)
		{
			leftParentMassTolerance = new Tolerance(leftParentMassTolerance.getValue(), isTolerancePPM);
			rightParentMassTolerance = new Tolerance(rightParentMassTolerance.getValue(), isTolerancePPM);
		}
	
		if(aaSet == null)
			aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		if(activationMethod == ActivationMethod.HCD)
			instType = InstrumentType.HIGH_RESOLUTION_LTQ;
		
		if(rightParentMassTolerance.getToleranceAsDa(1000) >= 0.5f)
			numAllowedC13 = 0;
		
		System.out.println("MS-GFDB v"+ VERSION + " (" + RELEASE_DATE + ")");
		
		runMSGFDB(specFile, specFormat, databaseFile, leftParentMassTolerance, rightParentMassTolerance, numAllowedC13,
	    		outputFile, enzyme, numAllowedNonEnzymaticTermini,
	    		activationMethod, instType, aaSet, numMatchesPerSpec, startSpecIndex, endSpecIndex, useTDA, showFDR,
	    		minPeptideLength, maxPeptideLength, minCharge, maxCharge, numThreads, useUniformAAProb, dbIndexDir, replicateMergedResults, doNotDseEdgeScore);
		System.out.format("MS-GFDB complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis()-time)/(float)1000);
	}
	
	public static void printUsageAndExit()
	{
		printUsageAndExit(null);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println("Error: " + message + "\n");
		System.out.println("MSGFDB v"+ VERSION + " (" + RELEASE_DATE + ")");
		System.out.print("Usage: java -Xmx2000M -jar MSGFDB.jar\n"
				+ "\t-s SpectrumFile (*.mzXML, *.mzML, *.mgf, *.ms2, *.pkl or *_dta.txt)\n" //, *.mgf, *.pkl, *.ms2)\n"
				+ "\t-d Database (*.fasta or *.fa)\n"
				+ "\t-t ParentMassTolerance (e.g. 2.5Da, 30ppm or 0.5Da,2.5Da)\n"
				+ "\t   Use comma to set asymmetric values. E.g. \"-t 0.5Da,2.5Da\" will set 0.5Da to the left (expMass<theoMass) and 2.5Da to the right (expMass>theoMass).\n"
				+ "\t[-o outputFileName] (Default: stdout)\n"
				+ "\t[-thread NumOfThreads] (Number of concurrent threads to be executed, Default: Number of available cores)\n"
				+ "\t[-tda 0/1] (0: don't search decoy database (default), 1: search decoy database to compute FDR)\n"
				+ "\t[-m FragmentationMethodID] (0: as written in the spectrum or CID if no info (Default), 1: CID, 2: ETD, 3: HCD, 4: Merge spectra from the same precursor)\n"
				+ "\t[-inst InstrumentID] (0: Low-res LCQ/LTQ (Default for CID and ETD), 1: TOF , 2: High-res LTQ (Default for HCD))\n"
				+ "\t[-e EnzymeID] (0: No enzyme, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: endogenous peptides)\n"
				+ "\t[-c13 0/1/2] (Number of allowed C13, Default: 1)\n"
				+ "\t[-nnet 0/1/2] (Number of allowed non-enzymatic termini, Default: 1)\n"
				+ "\t[-mod ModificationFileName] (Modification file, Default: standard amino acids with fixed C+57)\n"
				+ "\t[-minLength MinPepLength] (Minimum peptide length to consider, Default: 6)\n"
				+ "\t[-maxLength MaxPepLength] (Maximum peptide length to consider, Default: 40)\n"
				+ "\t[-minCharge MinPrecursorCharge] (Minimum precursor charge to consider if charges are not specified in the spectrum file, Default: 2)\n"
				+ "\t[-maxCharge MaxPrecursorCharge] (Maximum precursor charge to consider if charges are not specified in the spectrum file, Default: 3)\n"
				+ "\t[-n NumMatchesPerSpec] (Number of matches per spectrum to be reported, Default: 1)\n"
				+ "\t[-uniformAAProb 0/1] (0: use amino acid probabilities computed from the input database (Default), 1: use probability 0.05 for all amino acids)\n"
//				+ "\t[-scan scanNum] (scan number to be searched)\n" => hidden option
//				+ "\t[-showFDR 0/1] (0: don't show FDR, 1: show FDR (default)\n" => hidden option
//				+ "\t[-param paramFile]\n"
//				+ "\t[-err 0/1 (0: don't use peak errors (default), 1: use peak errors for scoring]\n"
//				+ "\t[-title 0/1] (0: don't show title (default), 1: show title)\n"
				);
		System.out.println("Example (high-precision): java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 30ppm -c13 1 -nnet 0 -tda 1 -o testMSGFDB.tsv");
		System.out.println("Example (low-precision): java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 0.5Da,2.5Da -nnet 0 -tda 1 -o testMSGFDB.tsv");
		System.exit(-1);
	}
	
    public static void runMSGFDB(
    		File specFile, 
    		SpecFileFormat specFormat, 
    		File databaseFile, 
    		Tolerance leftParentMassTolerance, 
    		Tolerance rightParentMassTolerance, 
    		int numAllowedC13,
    		File outputFile, 
    		Enzyme enzyme, 
    		int numAllowedNonEnzymaticTermini,
    		ActivationMethod activationMethod, 
    		InstrumentType instType,
    		AminoAcidSet aaSet, 
    		int numMatchesPerSpec,
    		int startSpecIndex,
    		int endSpecIndex,
    		boolean useTDA,
    		boolean showFDR,
    		int minPeptideLength,
    		int maxPeptideLength,
    		int minCharge,
    		int maxCharge,
    		int numThreads,
    		boolean useUniformAAProb,
    		File dbIndexDir,
    		boolean replicateMergedResults,
    		boolean doNotDseEdgeScore
    		)
	{
		long time = System.currentTimeMillis();
    	
		System.out.println("Loading database files...");
		if(dbIndexDir != null)
		{
			File newDBFile = new File(dbIndexDir.getPath()+File.separator+databaseFile.getName());
			if(!useTDA)
			{
				if(!newDBFile.exists())
				{
					System.out.println("Creating " + newDBFile.getPath() + ".");
					ReverseDB.copyDB(databaseFile.getPath(), newDBFile.getPath());
				}
			}
			databaseFile = newDBFile;
		}
		if(useTDA)
		{
			String dbFileName = databaseFile.getName();
			String concatDBFileName = dbFileName.substring(0, dbFileName.lastIndexOf('.'))+DECOY_DB_EXTENSION;
			File concatTargetDecoyDBFile = new File(databaseFile.getAbsoluteFile().getParent()+File.separator+concatDBFileName);
			if(!concatTargetDecoyDBFile.exists())
			{
				System.out.println("Creating " + concatTargetDecoyDBFile.getPath() + ".");
				if(ReverseDB.reverseDB(databaseFile.getPath(), concatTargetDecoyDBFile.getPath(), true) == false)
				{
					printUsageAndExit("Cannot create a decoy database file!");
				}
			}
			databaseFile = concatTargetDecoyDBFile;
		}
		
		if(!useUniformAAProb)
			DBScanner.setAminoAcidProbabilities(databaseFile.getPath(), aaSet);
		
		aaSet.registerEnzyme(enzyme);
		
		CompactFastaSequence fastaSequence = new CompactFastaSequence(databaseFile.getPath());
		if(useTDA)
		{
			float ratioUniqueProteins = fastaSequence.getRatioUniqueProteins();
			if(ratioUniqueProteins < 0.5f)
			{
				System.err.println("Error while indexing: " + databaseFile.getName() + " (too many redundant proteins)");
				System.err.println("If the database contains forward and reverse proteins, run MS-GFDB (or BuildSA) again with \"-tda 0\"");
				System.exit(-1);
			}
		}
		
		CompactSuffixArray sa = new CompactSuffixArray(fastaSequence, maxPeptideLength);
		System.out.print("Loading database finished ");
		System.out.format("(elapsed time: %.2f sec)\n", (float)(System.currentTimeMillis()-time)/1000);
		
		
		System.out.println("Reading spectra...");
    	Iterator<Spectrum> specItr = null;
		SpectrumAccessorBySpecIndex specMap = null;
		
		if(specFormat == SpecFileFormat.MZXML)
		{
			specItr = new MzXMLSpectraIterator(specFile.getPath());
			specMap = new MzXMLSpectraMap(specFile.getPath());
		}
		else if(specFormat == SpecFileFormat.DTA_TXT)
		{
			try {
				specItr = new PNNLSpectraIterator(specFile.getPath());
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			specMap = new PNNLSpectraMap(specFile.getPath());
		}
		else
		{
			SpectrumParser parser = null;
			if(specFormat == SpecFileFormat.MGF)
				parser = new MgfSpectrumParser();
			else if(specFormat == SpecFileFormat.MS2)
				parser = new MS2SpectrumParser();
			else if(specFormat == SpecFileFormat.PKL)
				parser = new PklSpectrumParser();
			
			try {
				specItr = new SpectraIterator(specFile.getPath(), parser);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			specMap = new SpectraMap(specFile.getPath(), parser);
		}
		
		if(specItr == null || specMap == null)
		{
			printUsageAndExit("Error while parsing spectrum file: " + specFile.getPath());
		}
		

		if(enzyme == null)
			numAllowedNonEnzymaticTermini = 2;
		
		// determine the number of spectra to be scanned together 
		long maxMemory = Runtime.getRuntime().maxMemory() - sa.getSize() - 1<<28;
		
		int avgPeptideMass = 2000;
		int numBytesPerMass = 12;
		int numSpecScannedTogether = (int)((float)maxMemory/avgPeptideMass/numBytesPerMass);
		ArrayList<SpecKey> specKeyList = SpecKey.getSpecKeyList(specItr, startSpecIndex, endSpecIndex, minCharge, maxCharge, activationMethod);
		int specSize = specKeyList.size();
		
		System.out.print("Reading spectra finished ");
		System.out.format("(elapsed time: %.2f sec)\n", (float)(System.currentTimeMillis()-time)/1000);
		
		numThreads = Math.min(numThreads, Math.round(Math.min(specSize, numSpecScannedTogether)/1000f));
		if(numThreads == 0)
			numThreads = 1;
		System.out.println("Using " + numThreads + (numThreads == 1 ? " thread." : " threads."));
		
		SpecDataType specDataType;
		if(!aaSet.containsPhosphorylation())
			specDataType = new SpecDataType(activationMethod, instType, enzyme);
		else
			specDataType = new SpecDataType(activationMethod, instType, enzyme, Modification.get("Phosphorylation"));
		int fromIndexGlobal = 0;
		
		String header = 
			"#SpecFile\tSpecIndex\tScan#\t"
			+"FragMethod\t"
//			+(showTitle ? "Title\t" : "")
			+"Precursor\tPMError("
			+(rightParentMassTolerance.isTolerancePPM() ? "ppm" : "Da")
			+")\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecProb\tP-value";
		MSGFDBResultGenerator gen = new MSGFDBResultGenerator(header);
		
		while(true)
		{
//			specSize = 728178;
			if(fromIndexGlobal >= specSize)
				break;
			int toIndexGlobal = Math.min(specSize, fromIndexGlobal+numSpecScannedTogether);
			System.out.println("Spectrum " + fromIndexGlobal + "-" + (toIndexGlobal-1) + " (total: " + specSize + ")");
			
			// Thread pool
			ExecutorService executor = Executors.newFixedThreadPool(numThreads);
			
			// Partition specKeyList
			int size = toIndexGlobal - fromIndexGlobal;
			int subListSize = size/numThreads;
			int residue = size % numThreads;
			
			int[] startIndex = new int[numThreads];
			int[] endIndex = new int[numThreads];
			
			for(int i=0; i<numThreads; i++)
			{
				startIndex[i] =  i > 0 ? endIndex[i-1] : fromIndexGlobal;
				endIndex[i] = startIndex[i] + subListSize + (i < residue ? 1 : 0);
			}
			
			for(int i=0; i<numThreads; i++)
			{
		    	ScoredSpectraMap specScanner = new ScoredSpectraMap(
		    			specMap,
		    			Collections.synchronizedList(specKeyList.subList(startIndex[i], endIndex[i])),
		    			leftParentMassTolerance,
		    			rightParentMassTolerance,
		    			numAllowedC13,
		    			specDataType
		    			);
		    	
				ConcurrentMSGFDB.RunMSGFDB msgfdbExecutor = new ConcurrentMSGFDB.RunMSGFDB(
							specScanner,
							sa,
							enzyme,
							aaSet,
							numMatchesPerSpec,
							minPeptideLength,
							maxPeptideLength,
							numAllowedNonEnzymaticTermini, 
							!useTDA,
							gen, 
							specFile.getName(),
							replicateMergedResults
							);
				executor.execute(msgfdbExecutor);
			}
			
			executor.shutdown();
			while(!executor.isTerminated()) {}	// wait until all threads terminate
			
			fromIndexGlobal += numSpecScannedTogether;
		}
		
    	time = System.currentTimeMillis();
		// Sort search results by spectral probabilities
		Collections.sort(gen);

		// Write results
		
    	if(showFDR && !useTDA && numMatchesPerSpec == 1)
    	{
    		PrintStream out = null;
    		if(outputFile == null)
    			out = System.out;
    		else
    		{
    			try {
    				out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
    			} catch (IOException e) {
    				e.printStackTrace();
    			}
    		}
    		System.out.println("Computing EFDRs...");
    		gen.computeEFDR();
    		System.out.print("Computing EFDRs finished");
    		System.out.format("(elapsed time: %.2f sec)\n", (float)(System.currentTimeMillis()-time)/1000);
    		gen.writeResults(out, true);
    		if(out != System.out)
    			out.close();
    	}
    	else if(!showFDR || !useTDA)
    	{
    		PrintStream out = null;
    		if(outputFile == null)
    			out = System.out;
    		else
    		{
    			try {
    				out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
    			} catch (IOException e) {
    				e.printStackTrace();
    			}
    		}
    		gen.writeResults(out, false);
    		if(out != System.out)
    			out.close();
    	}
    	else
    	{
    		System.out.println("Computing FDRs...");
    		try {
				File tempFile;
				if(outputFile != null)
					tempFile = new File(outputFile.getPath()+".temp.tsv");
				else
				{
					tempFile = File.createTempFile("MSGFDB", "tempResult");
					tempFile.deleteOnExit();
				}
				PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(tempFile)));
				gen.writeResults(out, false);
				out.flush();
				out.close();
				int specFileCol = 0;
				int specIndexCol = 1;
				int pepCol = 7;
				int dbCol = 8;
				int scoreCol = 11;
				fdr.ComputeFDR.computeFDR(tempFile, null, scoreCol, false, "\t", 
						specFileCol, specIndexCol, pepCol, null, true, true, 
						true, dbCol, "REV_",
						1, 1, outputFile);
				
			} catch (IOException e) {
				e.printStackTrace();
			}
    		System.out.print("Computing FDRs finished");
    		System.out.format("(elapsed time: %.2f sec)\n", (float)(System.currentTimeMillis()-time)/1000);
    	}
	}	    
}
