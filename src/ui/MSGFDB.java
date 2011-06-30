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

import parser.MS2SpectrumParser;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraIterator;
import parser.MzXMLSpectraMap;
import parser.PNNLSpectrumParser;
import parser.PklSpectrumParser;
import parser.SpectrumParser;
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
import msutil.SpecKey;
import msutil.SpectraIterator;
import msutil.SpectraMap;
import msutil.Spectrum;
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
		ActivationMethod activationMethod = null;
		InstrumentType instType = InstrumentType.LOW_RESOLUTION_LTQ;
		int numAllowedNonEnzymaticTermini = 1;
		boolean showTitle = false;
		boolean useTDA = false;
		boolean useUniformAAProb = false;
		int minPeptideLength = 6;
		int maxPeptideLength = 40;
		int minCharge = 2;
		int maxCharge = 3;
		int scanNum = -1;
		
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
			else if(argv[i].equalsIgnoreCase("-scan"))
			{
				try {
					scanNum = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal scanNum: " + argv[i+1]);
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
			else if(argv[i].equalsIgnoreCase("-title"))
			{
				if(argv[i+1].equalsIgnoreCase("1"))
					showTitle = true;
			}
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

		if(minPeptideLength > maxPeptideLength)
			printUsageAndExit("MinPepLength must not be larger than MaxPepLength!");
			
		if(minCharge > maxCharge)
			printUsageAndExit("MinPrecursorCharge must not be larger than MaxPrecursorCharge!");
		
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
		
		if(!useUniformAAProb)
			DBScanner.setAminoAcidProbabilities(databaseFile.getPath(), aaSet);
		
		aaSet.registerEnzyme(enzyme);
		
		runMSGFDB(specFile, specFormat, databaseFile, paramFile, parentMassTolerance, numAllowedC13,
	    		outputFile, enzyme, numAllowedNonEnzymaticTermini,
	    		activationMethod, instType, aaSet, numMatchesPerSpec, scanNum, showTitle, useTDA,
	    		minPeptideLength, maxPeptideLength, minCharge, maxCharge);
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
		System.out.println("MSGFDB v2 (06/29/2011)");
		System.out.print("Usage: java -Xmx2000M -jar MSGFDB.jar\n"
				+ "\t-s SpectrumFile (*.mzXML, *.mgf, *.ms2, *.pkl or *_dta.txt)\n" //, *.mgf, *.pkl, *.ms2)\n"
				+ "\t-d Database (*.fasta or *.fa)\n"
				+ "\t-t ParentMassTolerance (e.g. 2.5Da or 30ppm, no space is allowed.)\n"
				+ "\t[-o outputFileName] (Default: stdout)\n"
				+ "\t[-tda 0/1] (0: don't search decoy database (default), 1: search decoy database to compute FDR)\n"
				+ "\t[-m FragmentationMethodID] (0: as written in the spectrum or CID if no info (Default), 1: CID, 2: ETD, 3: HCD)\n"
				+ "\t[-inst InstrumentID] (0: Low-res LCQ/LTQ (Default for CID and ETD), 1: TOF , 2: High-res LTQ (Default for HCD))\n"
				+ "\t[-e EnzymeID] (0: No enzyme, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N)\n"
				+ "\t[-c13 0/1/2] (Number of allowed C13, Default: 1)\n"
				+ "\t[-nnet 0/1/2] (Number of allowed non-enzymatic termini, Default: 1)\n"
				+ "\t[-mod ModificationFileName] (Modification file, Default: standard amino acids with fixed C+57)\n"
				+ "\t[-minLength MinPepLength] (Minimum peptide length to consider, Default: 6)\n"
				+ "\t[-maxLength MaxPepLength] (Maximum peptide length to consider, Default: 40)\n"
				+ "\t[-minCharge MinPrecursorCharge] (Minimum precursor charge to consider if charges are not specified in the spectrum file, Default: 2)\n"
				+ "\t[-maxCharge MaxPrecursorCharge] (Maximum precursor charge to consider if charges are not specified in the spectrum file, Default: 3)\n"
				+ "\t[-n NumMatchesPerSpec] (Number of matches per spectrum to be reported, Default: 1)\n"
				+ "\t[-uniformAAProb 0/1] (0: use amino acid probabilities computed from the input database (default), 1: use probability 0.05 for all amino acids)\n"
//				+ "\t[-scan scanNum] (scan number to be searched)\n"
//				+ "\t[-param paramFile]\n"
//				+ "\t[-err 0/1 (0: don't use peak errors (default), 1: use peak errors for scoring]\n"
//				+ "\t[-title 0/1] (0: don't show title (default), 1: show title)\n"
				);
		System.out.println("Example: java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 30ppm -c13 1 -o testMSGFDB.tsv -tda 1");
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
    		int scanNum,
    		boolean showTitle,
    		boolean useTDA,
    		int minPeptideLength,
    		int maxPeptideLength,
    		int minCharge,
    		int maxCharge
    		)
	{
    	
    	Iterator<Spectrum> specItr = null;
		SpectrumAccessorByScanNum specAccessor = null;
		
		if(specFormat == SpecFileFormat.MZXML)
		{
			specItr = new MzXMLSpectraIterator(specFile.getPath());
			specAccessor = new MzXMLSpectraMap(specFile.getPath());
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
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			specAccessor = new SpectraMap(specFile.getPath(), parser);
		}
		
		if(specItr == null || specAccessor == null)
		{
			printUsageAndExit("Error while parsing spectrum file: " + specFile.getPath());
		}
		
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
		ArrayList<SpecKey> specKeyList = SpecKey.getSpecKeyList(specItr, minCharge, maxCharge);
		
		if(scanNum >= 0)
		{
			specKeyList.clear();
			Spectrum spec = specAccessor.getSpectrumByScanNum(scanNum);
			if(spec.getCharge() == 0)
			{
				for(int c=minCharge; c<=maxCharge; c++)
					specKeyList.add(new SpecKey(scanNum, c));
			}
			else
				specKeyList.add(new SpecKey(scanNum, spec.getCharge()));
		}
		
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
			if(fromIndex >= specKeyList.size())
				break;
			int toIndex = Math.min(specKeyList.size(), fromIndex+numSpecScannedTogether);
			System.out.println("Spectrum " + fromIndex + "-" + (toIndex-1) + " (total: " + specKeyList.size() + ")");
			
			// spectrum preprocessing
	    	System.out.print("Preprocessing spectra...");
	    	long time = System.currentTimeMillis();
	    	DBScanner sa = new DBScanner(
	    			new SuffixArraySequence(databaseFile.getPath()), 
	    			aaSet,
	    			specAccessor,
	    			specKeyList.subList(fromIndex, toIndex),
	    			parentMassTolerance,
	    			numAllowedC13,
	    			scorer,
	    			activationMethod,
	    			enzyme,
	    			numMatchesPerSpec,
	    			minPeptideLength,
	    			maxPeptideLength,
	    			minCharge,
	    			maxCharge
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

		
		System.out.print("Computing EFDRs...");
    	long time = System.currentTimeMillis();
		// Sort search results by spectral probabilities
		Collections.sort(gen);
    	if(!useTDA && numMatchesPerSpec == 1)
    	{
    		gen.computeEFDR();
    	}
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
				File tempFile;
				if(outputFile != null)
					tempFile = new File(outputFile.getPath()+".temp.tsv");
				else
				{
					tempFile = File.createTempFile("MSGFDB", "tempResult");
					tempFile.deleteOnExit();
				}
				PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(tempFile)));
//				tempFile.deleteOnExit();
				gen.writeResults(out, false);
				out.flush();
				out.close();
				int scanNumCol = 1;
				int pepCol = 6;
				int dbCol = 7;
				int scoreCol = 10;
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
