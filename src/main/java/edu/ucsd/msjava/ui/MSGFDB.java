package edu.ucsd.msjava.ui;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


import edu.ucsd.msjava.msdbsearch.CompactFastaSequence;
import edu.ucsd.msjava.msdbsearch.CompactSuffixArray;
import edu.ucsd.msjava.msdbsearch.ConcurrentMSGFDB;
import edu.ucsd.msjava.msdbsearch.DBScanner;
import edu.ucsd.msjava.msdbsearch.ReverseDB;
import edu.ucsd.msjava.msdbsearch.ScoredSpectraMap;
import edu.ucsd.msjava.msgf.MSGFDBResultGenerator;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.SpecFileFormat;
import edu.ucsd.msjava.msutil.SpecKey;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumAccessorBySpecIndex;
import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.IntRangeParameter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.params.ToleranceParameter;


public class MSGFDB {
	public static final String VERSION = "8091";
	public static final String RELEASE_DATE = "08/06/2012";
	
	public static final String DECOY_PROTEIN_PREFIX = "XXX";
	public static final String DECOY_DB_EXTENSION = ".revConcat.fasta";
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();

		ParamManager paramManager = new ParamManager("MSGFDB", MSGFDB.VERSION, MSGFDB.RELEASE_DATE, "java -Xmx2000M -jar MSGFDB.jar");
		paramManager.addMSGFDBParams();
		
		if(argv.length == 0)
		{
			paramManager.printUsageInfo();
			return;
		}
		
		// Parse parameters
		String errMessage = paramManager.parseParams(argv); 
		if(errMessage != null)
		{
			System.err.println("[Error] " + errMessage);
			System.out.println();
			paramManager.printUsageInfo();
			return;
		}
		
		// Running MS-GFDB
		paramManager.printToolInfo();
		String errorMessage = runMSGFDB(paramManager);
		if(errorMessage != null)
		{
			System.err.println("[Error] " + errorMessage);
			System.out.println();
		}
		else
			System.out.format("MS-GFDB complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis()-time)/(float)1000);
	}
	
    public static String runMSGFDB(ParamManager paramManager)
	{
		// Spectrum file
		FileParameter specParam = paramManager.getSpecFileParam();
		File specPath = specParam.getFile();
		
		if(!specPath.isDirectory())
		{
			// Spectrum format
			SpecFileFormat specFormat = (SpecFileFormat)specParam.getFileFormat();
			// Output file
			File outputFile = paramManager.getOutputFileParam().getFile();
			
			return runMSGFDB(specPath, specFormat, outputFile, paramManager);
		}
		else	// spectrum directory
		{
			for(File f : specPath.listFiles())
			{
				SpecFileFormat specFormat = SpecFileFormat.getSpecFileFormat(f.getName());
				if(specParam.isSupported(specFormat))
				{
					System.out.println("\nProcessing " + f.getAbsolutePath());
					String outputFileName = f.getName().substring(0, f.getName().lastIndexOf('.'))+".tsv";
					File outputFile = new File(outputFileName);
					if(outputFile.exists())
						return outputFile.getAbsolutePath() + " already exists!";
					System.out.println("Writing results to " + outputFile.getAbsolutePath());
					String errMsg = runMSGFDB(f, specFormat, outputFile, paramManager);
					if(errMsg != null)
						return errMsg;
				}
			}
			return null;
		}
	}
    
    private static String runMSGFDB(File specFile, SpecFileFormat specFormat, File outputFile, ParamManager paramManager)
    {
		long time = System.currentTimeMillis();
		
		// DB file
		File databaseFile = paramManager.getDBFileParam().getFile();
		
		// PM tolerance
		ToleranceParameter tol = ((ToleranceParameter)paramManager.getParameter("t"));
		Tolerance leftParentMassTolerance = tol.getLeftTolerance();
		Tolerance rightParentMassTolerance = tol.getRightTolerance();
		
		int toleranceUnit = paramManager.getIntValue("u");
		if(toleranceUnit != 2)
		{
			boolean isTolerancePPM;
			if(toleranceUnit == 0)
				isTolerancePPM = false;
			else 
				isTolerancePPM = true;
			leftParentMassTolerance = new Tolerance(leftParentMassTolerance.getValue(), isTolerancePPM);
			rightParentMassTolerance = new Tolerance(rightParentMassTolerance.getValue(), isTolerancePPM);
		}
		
		int numAllowedC13 = paramManager.getIntValue("c13");
		if(rightParentMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f)
			numAllowedC13 = 0;
		
		Enzyme enzyme = paramManager.getEnzyme();
		int numAllowedNonEnzymaticTermini = paramManager.getIntValue("nnet");
		ActivationMethod activationMethod = paramManager.getActivationMethod();
		InstrumentType instType = paramManager.getInstType();
		if(activationMethod == ActivationMethod.HCD)
			instType = InstrumentType.HIGH_RESOLUTION_LTQ;
		
		Protocol protocol = paramManager.getProtocol();
		
		AminoAcidSet aaSet = null;
		File modFile = paramManager.getModFileParam().getFile();
		if(modFile == null)
			aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		else
		{
			String modFileName = modFile.getName();
			String ext = modFileName.substring(modFileName.lastIndexOf('.')+1);
			if(ext.equalsIgnoreCase("xml"))
				aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(modFile.getAbsolutePath());
			else
				aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getAbsolutePath());
			if(aaSet.containsPhosphorylation())
			{
				protocol = Protocol.PHOSPHORYLATION;
			}
		}
		
		int numMatchesPerSpec = paramManager.getIntValue("n");
		
		int startSpecIndex = ((IntRangeParameter)paramManager.getParameter("index")).getMin();
		int endSpecIndex = ((IntRangeParameter)paramManager.getParameter("index")).getMax();
		
		boolean useTDA = paramManager.getIntValue("tda") == 1 ? true : false;
		boolean showFDR = paramManager.getIntValue("showFDR") == 1 ? true : false;
		boolean showDecoy = paramManager.getIntValue("showDecoy") == 1 ? true : false;
		
		int minPeptideLength = paramManager.getIntValue("minLength");
		int maxPeptideLength = paramManager.getIntValue("maxLength");
		if(minPeptideLength > maxPeptideLength)
		{
			return "MinPepLength must not be larger than MaxPepLength";
		}
		
		int minCharge = paramManager.getIntValue("minCharge");
		int maxCharge = paramManager.getIntValue("maxCharge");
		if(minCharge > maxCharge)
		{
			return "MinCharge must not be larger than MaxCharge";
		}
		
		int numThreads = paramManager.getIntValue("thread");
		boolean useUniformAAProb = paramManager.getIntValue("uniformAAProb") == 1 ? true : false;
		boolean replicateMergedResults = paramManager.getIntValue("replicate") == 1 ? true : false;
		boolean doNotDseEdgeScore = paramManager.getIntValue("edgeScore") == 1 ? true : false;
		
		System.out.println("Loading database files...");
		File dbIndexDir = paramManager.getFile("dd");
		if(dbIndexDir != null)
		{
			File newDBFile = new File(dbIndexDir.getAbsolutePath()+File.separator+databaseFile.getName());
			if(!useTDA)
			{
				if(!newDBFile.exists())
				{
					System.out.println("Creating " + newDBFile.getAbsolutePath() + ".");
					ReverseDB.copyDB(databaseFile.getAbsolutePath(), newDBFile.getAbsolutePath());
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
				System.out.println("Creating " + concatTargetDecoyDBFile.getAbsolutePath() + ".");
				if(ReverseDB.reverseDB(databaseFile.getAbsolutePath(), concatTargetDecoyDBFile.getAbsolutePath(), true, DECOY_PROTEIN_PREFIX) == false)
				{
					return "Cannot create a decoy database file!";
				}
			}
			databaseFile = concatTargetDecoyDBFile;
		}
		
		if(!useUniformAAProb)
			DBScanner.setAminoAcidProbabilities(databaseFile.getAbsolutePath(), aaSet);
		
		aaSet.registerEnzyme(enzyme);
		
		CompactFastaSequence fastaSequence = new CompactFastaSequence(databaseFile.getAbsolutePath()).truncateAnnotation();
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

		SpectraAccessor specAcc = new SpectraAccessor(specFile, specFormat);

		if(specAcc.getSpecMap() == null || specAcc.getSpecItr() == null)
			return "Error while parsing spectrum file: " + specFile.getPath();
		

		if(enzyme == null)
			numAllowedNonEnzymaticTermini = 2;
		
		// determine the number of spectra to be scanned together 
		long maxMemory = Runtime.getRuntime().maxMemory() - sa.getSize() - 1<<28;
		
		int avgPeptideMass = 2000;
		int numBytesPerMass = 12;
		int numSpecScannedTogether = (int)((float)maxMemory/avgPeptideMass/numBytesPerMass);
		ArrayList<SpecKey> specKeyList = SpecKey.getSpecKeyList(specAcc.getSpecItr(), startSpecIndex, endSpecIndex, minCharge, maxCharge, activationMethod);
		int specSize = specKeyList.size();
		
		System.out.print("Reading spectra finished ");
		System.out.format("(elapsed time: %.2f sec)\n", (float)(System.currentTimeMillis()-time)/1000);
		
		numThreads = Math.min(numThreads, Math.round(Math.min(specSize, numSpecScannedTogether)/1000f));
		if(numThreads == 0)
			numThreads = 1;
		System.out.println("Using " + numThreads + (numThreads == 1 ? " thread." : " threads."));
		
		SpecDataType specDataType = new SpecDataType(activationMethod, instType, enzyme, protocol);
		int fromIndexGlobal = 0;
		
		List<MSGFDBResultGenerator.DBMatch> resultList = Collections.synchronizedList(new ArrayList<MSGFDBResultGenerator.DBMatch>());
		
		while(true)
		{
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
		    			specAcc,
		    			Collections.synchronizedList(specKeyList.subList(startIndex[i], endIndex[i])),
		    			leftParentMassTolerance,
		    			rightParentMassTolerance,
		    			numAllowedC13,
		    			specDataType
		    			);
		    	if(doNotDseEdgeScore)
		    		specScanner.turnOffEdgeScoring();
		    	
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
							resultList, 
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
		Collections.sort(resultList);

		// Write results
		
		String header = 
			"#SpecFile\tSpecIndex\tScan#\t"
			+"FragMethod\t"
			+"Precursor\tPMError("
			+(rightParentMassTolerance.isTolerancePPM() ? "ppm" : "Da")
			+")\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecProb\tP-value";
		
		MSGFDBResultGenerator gen = new MSGFDBResultGenerator(header, resultList);
		
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
    		gen.writeResults(out, true, false);
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
    		gen.writeResults(out, false, false);
    		if(out != System.out)
    			out.close();
    	}
    	else
    	{
    		System.out.println("Computing FDRs...");
    		try {
				File tempFile = null;
				if(outputFile != null)
				{
					tempFile = new File(outputFile.getAbsolutePath()+".temp.tsv");
				}
				else
				{
					tempFile = File.createTempFile("MSGFDB", "tempResult");
					tempFile.deleteOnExit();
				}
				PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(tempFile)));
				gen.writeResults(out, false, false);
				out.flush();
				out.close();
				int specFileCol = 0;
				int specIndexCol = 1;
				int pepCol = 7;
				int dbCol = 8;
				int scoreCol = 11;
				edu.ucsd.msjava.fdr.ComputeFDR.computeFDR(tempFile, null, scoreCol, false, "\t", 
						specFileCol, specIndexCol, pepCol, null, true, showDecoy, 
						true, dbCol, DECOY_PROTEIN_PREFIX,
						1, 1, outputFile);
				
			} catch (IOException e) {
				e.printStackTrace();
			}
    		System.out.print("Computing FDRs finished");
    		System.out.format("(elapsed time: %.2f sec)\n", (float)(System.currentTimeMillis()-time)/1000);
    	}
    	return null;
	}	    
}
