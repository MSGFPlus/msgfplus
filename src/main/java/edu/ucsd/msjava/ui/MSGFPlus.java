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


import edu.ucsd.msjava.jmzparser.MzMLAdapter;
import edu.ucsd.msjava.msdbsearch.CompactFastaSequence;
import edu.ucsd.msjava.msdbsearch.CompactSuffixArray;
import edu.ucsd.msjava.msdbsearch.ConcurrentMSGFPlus;
import edu.ucsd.msjava.msdbsearch.DBScanner;
import edu.ucsd.msjava.msdbsearch.ReverseDB;
import edu.ucsd.msjava.msdbsearch.ScoredSpectraMap;
import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.DBSearchIOFiles;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.SpecFileFormat;
import edu.ucsd.msjava.msutil.SpecKey;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumAccessorBySpecIndex;
import edu.ucsd.msjava.mzid.MZIdentMLGen;
import edu.ucsd.msjava.params.ParamManager;


public class MSGFPlus {
	public static final String VERSION = "1.0 (v8226)";
	public static final String RELEASE_DATE = "08/20/2012";
	
	public static final String DECOY_DB_EXTENSION = ".revConcat.fasta";
	public static final String DECOY_PROTEIN_PREFIX = "XXX";
	
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		if(argv.length == 0)
		{
			paramManager.printUsageInfo();
			return;
		}
		
		MzMLAdapter.turnOffLogs();
		
		// Parse parameters
		String errMessage = paramManager.parseParams(argv); 
		if(errMessage != null)
		{
			System.err.println("[Error] " + errMessage);
			System.out.println();
			paramManager.printUsageInfo();
			return;
		}
		
		// Running MS-GF+
		paramManager.printToolInfo();
		String errorMessage = runMSGFPlus(paramManager);
		if(errorMessage != null)
		{
			System.err.println("[Error] " + errorMessage);
			System.out.println();
		}
		else
			System.out.format("MS-GF+ complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis()-time)/(float)1000);
	}
	
    public static String runMSGFPlus(ParamManager paramManager)
	{
    	SearchParams params = new SearchParams();
    	String errorMessage = params.parse(paramManager);
    	if(errorMessage != null)
    		return errorMessage;
    	else
    	{
    		List<DBSearchIOFiles> ioList = params.getDBSearchIOList();
    		boolean multFiles = false;
    		if(ioList.size() > 2)
    			multFiles = true;
    		
    		int ioIndex = -1;
    		for(DBSearchIOFiles ioFiles : ioList)
    		{
    			++ioIndex;
    			File specFile = ioFiles.getSpecFile();
    			SpecFileFormat specFormat = ioFiles.getSpecFileFormat();
    			File outputFile = ioFiles.getOutputFile();
    			
    			if(multFiles)
    			{
    				System.out.println("\nProcessing " + specFile.getPath());
    				System.out.println("Writing results to " + outputFile.getPath());
    			}
				String errMsg = runMSGFPlus(ioIndex, specFormat, outputFile, params);
				if(errMsg != null)
					return errMsg;
    		}
    	}
    	return null;
	}
    
    private static String runMSGFPlus(int ioIndex, SpecFileFormat specFormat, File outputFile, SearchParams params)
    {
		long time = System.currentTimeMillis();
		
		// DB file
		File databaseFile = params.getDatabaseFile();
		
		// PM tolerance
		Tolerance leftParentMassTolerance = params.getLeftParentMassTolerance();
		Tolerance rightParentMassTolerance = params.getRightParentMassTolerance();
		
		int min13C = params.getMin13C();	// inclusive
		int max13C = params.getMax13C();	// inclusive
		
		Enzyme enzyme = params.getEnzyme();
		
		ActivationMethod activationMethod = params.getActivationMethod();
		InstrumentType instType = params.getInstType();
		Protocol protocol = params.getProtocol();
		
		AminoAcidSet aaSet = params.getAASet();
		
		int startSpecIndex = params.getStartSpecIndex();
		int endSpecIndex = params.getEndSpecIndex();
		
		boolean useTDA = params.useTDA();
		
		int minCharge = params.getMinCharge();
		int maxCharge = params.getMaxCharge();
		
		int numThreads = params.getNumThreads();
		boolean doNotDseEdgeScore = params.doNotDseEdgeScore();
		
		System.out.println("Loading database files...");
		
		File dbIndexDir = params.getDBIndexDir();
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
				if(ReverseDB.reverseDB(databaseFile.getPath(), concatTargetDecoyDBFile.getPath(), true, DECOY_PROTEIN_PREFIX) == false)
				{
					return "Cannot create a decoy database file!";
				}
			}
			databaseFile = concatTargetDecoyDBFile;
		}
		
		DBScanner.setAminoAcidProbabilities(databaseFile.getPath(), aaSet);
		aaSet.registerEnzyme(enzyme);
		
		CompactFastaSequence fastaSequence = new CompactFastaSequence(databaseFile.getPath());
		if(useTDA)
		{
			float ratioUniqueProteins = fastaSequence.getRatioUniqueProteins();
			if(ratioUniqueProteins < 0.5f)
			{
				System.err.println("Error while indexing: " + databaseFile.getName() + " (too many redundant proteins)");
				System.err.println("If the database contains forward and reverse proteins, run MS-GF+ (or BuildSA) again with \"-tda 0\"");
				System.exit(-1);
			}
		}
		
		CompactSuffixArray sa = new CompactSuffixArray(fastaSequence, params.getMaxPeptideLength());
		System.out.print("Loading database finished ");
		System.out.format("(elapsed time: %.2f sec)\n", (float)(System.currentTimeMillis()-time)/1000);
		
		System.out.println("Reading spectra...");
		
		File specFile = params.getDBSearchIOList().get(ioIndex).getSpecFile();
		
		SpectraAccessor specAcc = new SpectraAccessor(specFile, specFormat);

    	Iterator<Spectrum> specItr = null;
    	try {
			specItr = specAcc.getSpecItr();
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		SpectrumAccessorBySpecIndex specMap = specAcc.getSpecMap();
		
		if(specItr == null || specMap == null)
		{
			return "Error while parsing spectrum file: " + specFile.getPath();
		}

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
		
		SpecDataType specDataType = new SpecDataType(activationMethod, instType, enzyme, protocol);
		int fromIndexGlobal = 0;
		
		MZIdentMLGen mzidGen = new MZIdentMLGen(params, aaSet, sa, specAcc, ioIndex);
		
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
		    			specMap,
		    			Collections.synchronizedList(specKeyList.subList(startIndex[i], endIndex[i])),
		    			leftParentMassTolerance,
		    			rightParentMassTolerance,
		    			min13C,
		    			max13C,
		    			specDataType
		    			);
		    	if(doNotDseEdgeScore)
		    		specScanner.turnOffEdgeScoring();
		    	
				ConcurrentMSGFPlus.RunMSGFPlus msgfdbExecutor = new ConcurrentMSGFPlus.RunMSGFPlus(
							specScanner,
							sa,
							params,
							mzidGen
							);
				executor.execute(msgfdbExecutor);
			}
			
			executor.shutdown();
			while(!executor.isTerminated()) {}	// wait until all threads terminate
			
			fromIndexGlobal += numSpecScannedTogether;
		}
		
    	time = System.currentTimeMillis();
    	
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

    	mzidGen.writeResults(out);

    	out.close();
    	
    	return null;
	}	    
}
