package ui;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import parser.MgfSpectrumParser;
import parser.MzXMLSpectraMap;
import suffixarray.SuffixArraySequence;

import msdbsearch.DBScanner;
import msgf.MSGFDBResultGenerator;
import msgf.NominalMassFactory;
import msgf.Tolerance;
import msscorer.NewRankScorer;
import msscorer.NewScorerFactory;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.ActivationMethod;
import msutil.SpecFileFormat;
import msutil.SpectraMap;
import msutil.SpectrumAccessorByScanNum;

public class MSGFDB {
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();
		if(argv.length < 2 || argv.length % 2 != 0)
			printUsageAndExit();

		File 	specFile = null;
		SpecFileFormat specFormat = null;
		File 	databaseFile 	= null;
		File	paramFile		= null;
		File	outputFile = null;
		Tolerance	parentMassTolerance = null;
		int numAllowedC13 = 0;
		int 	numMatchesPerSpec = 1;
		Enzyme	enzyme = Enzyme.TRYPSIN;
		ActivationMethod activationMethod = null;
		int numAllowedNonEnzymaticTermini = 0;
		boolean showTitle = false;
		boolean useError = false;
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
					printUsageAndExit("Illigal enzyme: " + argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-mod"))
			{
				File modFile = new File(argv[i+1]);
				if(!modFile.exists())
				{
					printUsageAndExit(modFile + " doesn't exist.");
				}
				aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
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
			else if(argv[i].equalsIgnoreCase("-err"))
			{
				if(argv[i+1].equalsIgnoreCase("1"))
					useError = true;
			}
			else
				printUsageAndExit();
		}
		
		if(specFile == null)
			printUsageAndExit("Spectrum is not specified.");
		if(databaseFile == null)
			printUsageAndExit("Database is not specified.");
		
		if(parentMassTolerance == null)
			printUsageAndExit("Parent mass tolerance is not specified.");
	
		if(aaSet == null)
			aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		DBScanner.setAminoAcidProbabilities(databaseFile.getPath(), aaSet);
		////////// Debug ////////////
//		aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		DBScanner.setAminoAcidProbabilities("/home/sangtaekim/Research/Data/CommonContaminants/IPI_human_3.79_withContam.fasta", aaSet);
//		DBScanner.setAminoAcidProbabilities("/home/sangtaekim/Research/Data/ABRF/StudyFiles/iPRG2011CCS.fasta", aaSet);
		//////////////////
		
		aaSet.registerEnzyme(enzyme);
		
		runMSGFDB(specFile, specFormat, databaseFile, paramFile, useError, parentMassTolerance, numAllowedC13,
	    		outputFile, enzyme, numAllowedNonEnzymaticTermini,
	    		activationMethod, aaSet, numMatchesPerSpec, showTitle,
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
			System.err.println(message);
		System.err.println("MSGFDB v20100323");
		System.out.print("usage: java -Xmx3500M -jar MSGFDB.jar\n"
				+ "\t-s spectrumFile (*.mzXML or *.mgf)\n" //, *.mgf, *.pkl, *.ms2)\n"
				+ "\t-d database (*.fasta)\n"
				+ "\t-t parentMassTolerance (e.g. 2.5Da or 50ppm, no space is allowed.)\n"
				+ "\t[-c13 0/1/2] (Number of allowed C13, default: 0)\n"
				+ "\t[-m FragmentationMethodID] (0: as written in the spectrum, 1: CID , 2: ETD, 3: HCD)\n"//, 3: CID/ETD pair)\n"
				+ "\t[-e EnzymeID] (0: No enzyme, 1: Trypsin (default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N)\n"
				+ "\t[-nnet 0/1/2] (Number of allowed non-enzymatic termini, default: 0)\n"
				+ "\t[-param paramFile]\n"
				+ "\t[-o outputFileName] (Default: stdout)\n"
//				+ "\t[-n numMatchesPerSpec (Default: 1)]\n"
				+ "\t[-err 0/1 (0: don't use peak errors (default), 1: use peak errors for scoring]\n"
				+ "\t[-mod modificationFileName (default: standard amino acids with fixed C+57)]\n"
				+ "\t[-title 0/1] (0: don't show title (default), 1: show title)\n"
				+ "\t[-minLength minPepLength] (default: 6)\n"
				+ "\t[-maxLength maxPepLength] (default: 40)\n"
				);
		System.exit(-1);
	}
	
    public static void runMSGFDB(
    		File specFile, 
    		SpecFileFormat specFormat, 
    		File databaseFile, 
    		File paramFile, 
    		boolean useError,
    		Tolerance parentMassTolerance, 
    		int numAllowedC13,
    		File outputFile, 
    		Enzyme enzyme, 
    		int numAllowedNonEnzymaticTermini,
    		ActivationMethod activationMethod,  
    		AminoAcidSet aaSet, 
    		int numMatchesPerSpec,
    		boolean showTitle,
    		int minPeptideLength,
    		int maxPeptideLength
    		)
	{
		SpectrumAccessorByScanNum specAccessor = null;
		if(specFormat == SpecFileFormat.MZXML)
			specAccessor = new MzXMLSpectraMap(specFile.getPath());
		else if(specFormat == SpecFileFormat.MGF)
			specAccessor = new SpectraMap(specFile.getPath(), new MgfSpectrumParser());
		
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
		
		NewRankScorer scorer = null;
		if(paramFile != null)
			scorer = new NewRankScorer(paramFile.getPath());
		else if(activationMethod != null)
			scorer = NewScorerFactory.get(activationMethod, enzyme);

		if(!useError && scorer != null)
			scorer.doNotUseError();
		
		NominalMassFactory factory = new NominalMassFactory(aaSet, enzyme, maxPeptideLength);
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
//		scanNumList.add(3888);
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
			+(scorer == null ? "FragMethod\t" : "")
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
	    			factory,
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
	    			maxPeptideLength,
	    			useError
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
	    	sa.computeSpecProb(true);
	    	System.out.println(" " + (System.currentTimeMillis()-time)/(float)1000 + " sec");
	    	
			System.out.print("Generating results...");
	    	time = System.currentTimeMillis(); 
	    	sa.addDBSearchResults(gen, specFile.getName(), showTitle);
	    	System.out.println(" " + (System.currentTimeMillis()-time)/(float)1000 + " sec");

			fromIndex += numSpecScannedTogether;
		}

		System.out.print("Computing EFDRs...");
    	long time = System.currentTimeMillis(); 
    	gen.computeEFDR();
    	System.out.println(" " + (System.currentTimeMillis()-time)/(float)1000 + " sec");
    	
		System.out.print("Writing results...");
    	time = System.currentTimeMillis(); 
    	gen.writeResults(out);
    	System.out.println(" " + (System.currentTimeMillis()-time)/(float)1000 + " sec");
		
		if(out != System.out)
			out.close();
	}	
}
