package edu.ucsd.msjava.ui;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;


import edu.ucsd.msjava.msdictionary.MSDicLauncher;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewAdditiveScorer;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.params.ParamParser;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;
import edu.ucsd.msjava.parser.PklSpectrumParser;
import edu.ucsd.msjava.parser.SpectrumParser;
import edu.ucsd.msjava.suffixarray.SuffixArray;
import edu.ucsd.msjava.suffixarray.SuffixArraySequence;

public class MSDictionary {
	public static void main(String argv[])
	{
		if(argv.length != 2 && argv.length != 4)
		{
			printUsageAndExit();
			System.exit(-1);
		}
		String paramFileName = null;
		String outputFileName = null;
		for(int i=0; i<argv.length; i+=2)
		{ 
			if(!argv[i].startsWith("-") || i+1 >= argv.length)
				printUsageAndExit();
			if(argv[i].equalsIgnoreCase("-i"))
			{
				paramFileName = argv[i+1];
			}
			else if(argv[i].equalsIgnoreCase("-o"))
			{
				outputFileName = argv[i+1];
			}
		}
		if(paramFileName == null)
		{
			System.out.println("Error: parameter file is missing.");
			printUsageAndExit();
		}
		
		runMSDictionary(paramFileName, outputFileName);
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("MS-Dictionary (v.20100201)\n" +
				"usage: java -jar MSDictionary.jar -i paramFile [-o outputFile]\n" +
				"(example: java -Xmx3000M -jar MSDictionary.jar -i sampleInput.txt -o test.txt)");
		System.exit(-1);
	}
	
	// default parameters
	public static void runMSDictionary(String paramFile, String outputFileName)
	{
		ParamParser.Parameters params = ParamParser.parseFromFile(paramFile);

		// Spectrum
		String specFileName = params.getParameter("Spectrum");
		if(!new File(specFileName).exists())
			printParsingErrorAndExit(specFileName + " doesn't exist.");
		Iterator<Spectrum> specIterator = null;
		String ext = specFileName.substring(specFileName.lastIndexOf('.')+1);
		if(ext.equalsIgnoreCase("mzxml"))
			specIterator = new MzXMLSpectraIterator(specFileName);
		else
		{
			SpectrumParser parser = null;
			if(ext.equalsIgnoreCase("mgf"))
				parser = new MgfSpectrumParser();
			else if(ext.equalsIgnoreCase("pkl"))
				parser = new PklSpectrumParser();
			else
				printParsingErrorAndExit(ext+" format is not supported");
			
			try {
				specIterator = new SpectraIterator(specFileName, parser);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		// Database
		String dbFileName = params.getParameter("Database");
		SuffixArray sa = null;
		if(dbFileName != null)	// if no database is specified, just generate reconstructions
		{
			if(!new File(dbFileName).exists())
				printParsingErrorAndExit(dbFileName + " doesn't exist.");
			sa = new SuffixArray(new SuffixArraySequence(dbFileName, edu.ucsd.msjava.sequences.Constants.AMINO_ACIDS_19));			
		}

		// Scoring parameters
		String scoringParamFile = params.getParameter("ScoringParams");
		NewAdditiveScorer scorer;
		if(scoringParamFile == null)
			scorer = NewScorerFactory.get(ActivationMethod.CID, Enzyme.TRYPSIN);
		else
			scorer = new NewRankScorer(scoringParamFile);
		
		MSDicLauncher msDicLauncher = new MSDicLauncher(specIterator, scorer, sa);
		
		// Parent mass tolerance
		Tolerance pmTolerance = null;
		String pmTolStr = params.getParameter("PMTolerance");
		if(pmTolStr != null)
			pmTolerance = Tolerance.parseToleranceStr(pmTolStr);
		if(pmTolerance == null)
			printParsingErrorAndExit("Input file parsing error: illegal parent mass tolerance.");
		else
			msDicLauncher.pmTolerance(pmTolerance);
		
		// Fragment mass tolerance
//		Tolerance fragTolerance = null;
//		String fragTolStr = params.getParameter("Tolerance");
//		if(fragTolStr != null)
//			fragTolerance = Tolerance.parseToleranceStr(fragTolStr);
//		if(fragTolerance == null)
//			printParsingErrorAndExit("Input file parsing error: illegal fragment mass tolerance.");
//		else
//			msDicLauncher.fragTolerance(fragTolerance);
			
		// Spectral probability
		Float specProb = null;
		String specProbStr = params.getParameter("SpecProb");
		if(specProbStr != null)
			specProb = Float.parseFloat(specProbStr);
		if(specProb == null)
			printParsingErrorAndExit("Input file parsing error: illegal spectral probability.");
		else
			msDicLauncher.specProb(specProb);
		
		// Number of reconstructions
		Float numRecs = null;
		String numRecsStr = params.getParameter("NumRecs");
		if(numRecsStr != null)
			numRecs = Float.parseFloat(numRecsStr);
		if(numRecs == null)
			printParsingErrorAndExit("Input file parsing error: illegal spectral probability.");
		else
			msDicLauncher.numRecs(numRecs);
		
		Integer isNumInclusive = params.getIntParameter("IsNumInclusive");
		if(isNumInclusive == null)
			printParsingErrorAndExit("Input file parsing error: illegal IsNumInclusive field.");
		else if(isNumInclusive == 1)
			msDicLauncher.setNumInclusive();
		
		Integer isTrypticOnly = params.getIntParameter("IsTrypticOnly");
		if(isTrypticOnly == null)
			printParsingErrorAndExit("Input file parsing error: illegal IsTrypticOnly field.");
		else if(isTrypticOnly == 0)
			msDicLauncher.allowNonTryptic();
		
		Integer msgfThreshold = params.getIntParameter("MSGFThreshold");
		if(msgfThreshold == null)
			printParsingErrorAndExit("Input file parsing error: illegal MSGFThreshold field.");
		else
			msDicLauncher.msgfScoreThreshold(msgfThreshold);
		
		Float minParentMass = params.getFloatParameter("MinParentMass");
		if(minParentMass == null)
			printParsingErrorAndExit("Input file parsing error: illegal minimum parent mass.");
		else
			msDicLauncher.minParentMass(minParentMass);
			
		Float maxParentMass = params.getFloatParameter("MaxParentMass");
		if(maxParentMass == null)
			printParsingErrorAndExit("Input file parsing error: illegal maximum parent mass.");
		else
			msDicLauncher.maxParentMass(maxParentMass);
		
		if(outputFileName != null)
			msDicLauncher.outputFileName(outputFileName);
		
		msDicLauncher.runMSDictionary();
	}
	
	
	public static void printParsingErrorAndExit(String message)
	{
		System.err.println(message);
		System.exit(-1);
	}
}
