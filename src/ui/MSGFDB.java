package ui;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import fdr.Pair;

import parser.MgfSpectrumParser;
import parser.MzXMLSpectraMap;
import suffixarray.SuffixArraySequence;

import msdbsearch.DBScanner;
import msdbsearch.ReverseDB;
import msgf.MSGFDBResultGenerator;
import msgf.Tolerance;
import msscorer.NewRankScorer;
import msscorer.NewScorerFactory;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.ActivationMethod;
import msutil.InstrumentType;
import msutil.SpecFileFormat;
import msutil.SpectraMap;
import msutil.SpectrumAccessorByScanNum;

public class MSGFDB {
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();
		if(argv.length < 2 || argv.length % 2 != 0)
			printUsageAndExit("The number of parameters must be even.");

		File 	specFile = null;
		SpecFileFormat specFormat = null;
		File 	databaseFile 	= null;
		File	paramFile		= null;
		File	outputFile = null;
		Tolerance parentMassTolerance = null;
		Boolean isTolerancePPM = null;
		int numAllowedC13 = 1;
		int 	numMatchesPerSpec = 1;
		Enzyme	enzyme = Enzyme.TRYPSIN;
		ActivationMethod activationMethod = ActivationMethod.CID;
		InstrumentType instType = InstrumentType.LOW_RESOLUTION_LTQ;
		int numAllowedNonEnzymaticTermini = 1;
		boolean showTitle = false;
		boolean useTDA = false;
		int minPeptideLength = 6;
		int maxPeptideLength = 40;
		
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
					int posDot = specFile.getName().lastIndexOf('.');
					if(posDot >= 0)
					{
						String extension = specFile.getName().substring(posDot);
						if(extension.equalsIgnoreCase(".mzXML"))
							specFormat = SpecFileFormat.MZXML;
						else if(extension.equalsIgnoreCase(".mgf"))
							specFormat = SpecFileFormat.MGF;
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
				parentMassTolerance = Tolerance.parseToleranceStr(argv[i+1]);
				if(parentMassTolerance == null)
				{
					printUsageAndExit("Illegal tolerance value: " + argv[i+1]);
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
				if(argv[i+1].equalsIgnoreCase("1"))
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
			else if(argv[i].equalsIgnoreCase("-title"))
			{
				if(argv[i+1].equalsIgnoreCase("1"))
					showTitle = true;
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
//			else if(argv[i].equalsIgnoreCase("-err"))
//			{
//				if(argv[i+1].equalsIgnoreCase("1"))
//					useError = true;
//			}
			else
			{
				printUsageAndExit("Invalid option: " + argv[i]);
			}
		}
		
		if(specFile == null)
			printUsageAndExit("Spectrum is not specified.");
		if(databaseFile == null)
			printUsageAndExit("Database is not specified.");
		
		if(parentMassTolerance == null)
			printUsageAndExit("Parent mass tolerance is not specified.");

		if(isTolerancePPM != null)
			parentMassTolerance = new Tolerance(parentMassTolerance.getValue(), isTolerancePPM);
	
		if(aaSet == null)
			aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		if(activationMethod == ActivationMethod.HCD)
			instType = InstrumentType.HIGH_RESOLUTION_LTQ;
		
		if(useTDA)
		{
			String dbFileName = databaseFile.getName();
			String concatDBFileName = dbFileName.substring(0, dbFileName.lastIndexOf('.'))+".revConcat.fasta";
			File concatTargetDecoyDBFile = new File(databaseFile.getAbsoluteFile().getParent()+File.separator+concatDBFileName);
			if(!concatTargetDecoyDBFile.exists())
			{
				System.out.println("Creating " + concatDBFileName + ".");
				if(ReverseDB.reverseDB(databaseFile.getPath(), concatTargetDecoyDBFile.getPath(), true) == false)
				{
					printUsageAndExit("Cannot create decoy database file!");
				}
			}
			databaseFile = concatTargetDecoyDBFile;
		}
		
		DBScanner.setAminoAcidProbabilities(databaseFile.getPath(), aaSet);
		////////// Debug ////////////
//		aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		DBScanner.setAminoAcidProbabilities("/home/sangtaekim/Research/Data/CommonContaminants/IPI_human_3.79_withContam.fasta", aaSet);
//		DBScanner.setAminoAcidProbabilities("/home/sangtaekim/Research/Data/ABRF/StudyFiles/iPRG2011CCS.fasta", aaSet);
		//////////////////
		
		aaSet.registerEnzyme(enzyme);
		
		runMSGFDB(specFile, specFormat, databaseFile, paramFile, parentMassTolerance, numAllowedC13,
	    		outputFile, enzyme, numAllowedNonEnzymaticTermini,
	    		activationMethod, instType, aaSet, numMatchesPerSpec, showTitle, useTDA,
	    		minPeptideLength, maxPeptideLength);
		System.out.format("Time: %.3f sec\n", (System.currentTimeMillis()-time)/(float)1000);
	}
	
	public static void printUsageAndExit()
	{
		printUsageAndExit(null);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println(message);
		System.out.println("MSGFDB v2 (06/15/2011)");
		System.out.print("usage: java -Xmx3500M -jar MSGFDB.jar\n"
				+ "\t-s SpectrumFile (*.mzXML or *.mgf)\n" //, *.mgf, *.pkl, *.ms2)\n"
				+ "\t-d Database (*.fasta)\n"
				+ "\t-t ParentMassTolerance (e.g. 2.5Da or 30ppm, no space is allowed.)\n"
				+ "\t[-o outputFileName] (Default: stdout)\n"
				+ "\t[-tda 0/1 (0: don't search decoy database (default), 1: search decoy database to compute FDR]\n"
				+ "\t[-m FragmentationMethodID] (0: as written in the spectrum, 1: CID (Default), 2: ETD, 3: HCD)\n"//, 3: CID/ETD pair)\n"
				+ "\t[-inst InstrumentID] (0: Low-res LCQ/LTQ (Default for CID and ETD), 1: TOF , 2: High-res LTQ (Default for HCD))\n"
				+ "\t[-e EnzymeID] (0: No enzyme, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N)\n"
				+ "\t[-c13 0/1/2] (Number of allowed C13, Default: 1)\n"
				+ "\t[-nnet 0/1/2] (Number of allowed non-enzymatic termini, Default: 1)\n"
				+ "\t[-mod modificationFileName (Default: standard amino acids with fixed C+57)]\n"
				+ "\t[-minLength minPepLength] (Default: 6)\n"
				+ "\t[-maxLength maxPepLength] (Default: 40)\n"
				+ "\t[-n numMatchesPerSpec (Default: 1)]\n"
//				+ "\t[-param paramFile]\n"
//				+ "\t[-err 0/1 (0: don't use peak errors (default), 1: use peak errors for scoring]\n"
//				+ "\t[-title 0/1] (0: don't show title (default), 1: show title)\n"
				);
		System.exit(-1);
	}
	
    public static void runMSGFDB(
    		File specFile, 
    		SpecFileFormat specFormat, 
    		File databaseFile, 
    		File paramFile, 
    		Tolerance parentMassTolerance, 
    		int numAllowedC13,
    		File outputFile, 
    		Enzyme enzyme, 
    		int numAllowedNonEnzymaticTermini,
    		ActivationMethod activationMethod,  
    		InstrumentType instType,
    		AminoAcidSet aaSet, 
    		int numMatchesPerSpec,
    		boolean showTitle,
    		boolean useTDA,
    		int minPeptideLength,
    		int maxPeptideLength
    		)
	{
		SpectrumAccessorByScanNum specAccessor = null;
		if(specFormat == SpecFileFormat.MZXML)
			specAccessor = new MzXMLSpectraMap(specFile.getPath());
		else if(specFormat == SpecFileFormat.MGF)
			specAccessor = new SpectraMap(specFile.getPath(), new MgfSpectrumParser());
		
		NewRankScorer scorer = null;
		if(paramFile != null)
			scorer = new NewRankScorer(paramFile.getPath());
		else if(activationMethod != null)
			scorer = NewScorerFactory.get(activationMethod, instType, enzyme);

		if(enzyme == null)
			numAllowedNonEnzymaticTermini = 2;
		
//		NominalMassFactory factory = new NominalMassFactory(aaSet, enzyme, maxPeptideLength);
		
		// determine the number of spectra to be scanned together 
		long maxMemory = Runtime.getRuntime().maxMemory();
		int avgPeptideMass = 2000;
		int numBytesPerMass = 8;
		int numSpecScannedTogether = (int)((float)maxMemory/avgPeptideMass/numBytesPerMass);
		ArrayList<Integer> scanNumList = specAccessor.getScanNumList();
		
		////////// Debug ////////////
//		aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		DBScanner.setAminoAcidProbabilities("/home/sangtaekim/Research/Data/CommonContaminants/IPI_human_3.79_withContam.fasta", aaSet);
//		aaSet.printAASet();
//		scanNumList.clear();
//		scanNumList.add(21678);
//		scanNumList.add(31669);
//		scanNumList.add(3888);
//		scanNumList.add(3256);
//		scanNumList.add(6416);
//		scanNumList.add(3751);
//		scanNumList.add(338);	// decoytest
//		scanNumList.add(857);
//		scanNumList.add(685);	// Q-17
//		scanNumList.add(4378);	// Q-17
//		scanNumList.add(1162);	// M+16
//		scanNumList.add(3888);
//		for(int sn=1000; sn<1100; sn++)
//			scanNumList.add(sn);
		////////////////////////////////
		
		int fromIndex = 0;
		
		String header = 
			"#SpecFile\tScan#\t"
			+"FragMethod\t"
			+(showTitle ? "Title\t" : "")
			+"Precursor\tPMError("
			+(parentMassTolerance.isTolerancePPM() ? "ppm" : "Da")
			+")\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecProb\tP-value";
		MSGFDBResultGenerator gen = new MSGFDBResultGenerator(header);	
		while(true)
		{
			if(fromIndex >= scanNumList.size())
				break;
			int toIndex = Math.min(scanNumList.size(), fromIndex+numSpecScannedTogether);
			System.out.println("Spectrum " + fromIndex + "-" + (toIndex-1) + " (total: " + scanNumList.size() + ")");
			
			// spectrum preprocessing
	    	System.out.print("Preprocessing spectra...");
	    	long time = System.currentTimeMillis();
	    	DBScanner sa = new DBScanner(
	    			new SuffixArraySequence(databaseFile.getPath()), 
	    			aaSet,
	    			specAccessor,
	    			scanNumList.subList(fromIndex, toIndex),
	    			parentMassTolerance,
	    			numAllowedC13,
	    			scorer,
	    			activationMethod,
	    			enzyme,
	    			numMatchesPerSpec,
	    			minPeptideLength,
	    			maxPeptideLength
	    			);
	    	System.out.println(" " + (System.currentTimeMillis()-time)/(float)1000 + " sec");
	    	
	    	// db search
	    	time = System.currentTimeMillis(); 
			if(enzyme == null)
				sa.dbSearchNoEnzyme(true);	// currently not supported
			else if(enzyme.isCTerm())
			{
				if(!aaSet.containsModification())
					sa.dbSearchCTermEnzymeNoMod(numAllowedNonEnzymaticTermini, true);
				else
					sa.dbSearchCTermEnzyme(numAllowedNonEnzymaticTermini, true);
			}
			else
				sa.dbSearchNTermEnzyme(numAllowedNonEnzymaticTermini, true);
			
			System.out.println("Database search... " + (System.currentTimeMillis()-time)/(float)1000 + " sec");

	    	// running MS-GF
			System.out.print("Computing P-values...");
	    	time = System.currentTimeMillis(); 
	    	sa.computeSpecProb(!useTDA);
	    	System.out.println(" " + (System.currentTimeMillis()-time)/(float)1000 + " sec");
	    	
			System.out.print("Generating results...");
	    	time = System.currentTimeMillis(); 
	    	sa.addDBSearchResults(gen, specFile.getName(), showTitle);
	    	System.out.println(" " + (System.currentTimeMillis()-time)/(float)1000 + " sec");
			fromIndex += numSpecScannedTogether;
		}

		// Sort search results
		Collections.sort(gen);
		
		System.out.print("Computing EFDRs...");
    	long time = System.currentTimeMillis();
    	if(!useTDA && numMatchesPerSpec == 1)
    		gen.computeEFDR();
    	System.out.println(" " + (System.currentTimeMillis()-time)/(float)1000 + " sec");
    	
		System.out.print("Writing results...");
    	time = System.currentTimeMillis(); 
    	if(!useTDA)
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
    		gen.writeResults(out, numMatchesPerSpec == 1);
    		if(out != System.out)
    			out.close();
    	}
    	else
    	{
    		try {
//				File tempFile = File.createTempFile("MSGFDB", "tempResult");
				File tempFile = new File(outputFile.getPath()+".temp.tsv");
				PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(tempFile)));
//				tempFile.deleteOnExit();
				gen.writeResults(out, false);
				int scanNumCol = 1;
				int pepCol = 6;
				int dbCol = 7;
				int scoreCol = 11;
				if(showTitle)
				{
					++pepCol;
					++dbCol;
					++scoreCol;
				}
				fdr.ComputeFDR.computeFDR(tempFile, null, scoreCol, false, "\t", 
						scanNumCol, pepCol, null, true, true, 
						true, dbCol, "REV_",
						1, 1, outputFile);
//				tempFile.delete();
				
			} catch (IOException e) {
				e.printStackTrace();
			}
    	}
    	System.out.println(" " + (System.currentTimeMillis()-time)/(float)1000 + " sec");
		
	}	
}
