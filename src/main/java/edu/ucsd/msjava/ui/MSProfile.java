package edu.ucsd.msjava.ui;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import edu.ucsd.msjava.msgf.AminoAcidGraph;
import edu.ucsd.msjava.msgf.GeneratingFunction;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msgf.NominalMassFactory;
import edu.ucsd.msjava.msgf.Profile;
import edu.ucsd.msjava.msgf.ProfileGF;
import edu.ucsd.msjava.msgf.ProfilePeak;
import edu.ucsd.msjava.msgf.ScoredSpectrum;
import edu.ucsd.msjava.msgf.ScoredSpectrumSum;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Constants;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Sequence;
import edu.ucsd.msjava.msutil.SpecFileFormat;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;

public class MSProfile {
	public static float MIN_PROF_PROB_REPORT_THRESHOLD = 0.01f;
	public static void main(String argv[])
	{
		if(argv.length < 2 || argv.length % 2 != 0)
			printUsageAndExit("Illegal parameters");

		// required
		File 	specFile = null;
		SpecFileFormat specFileFormat = null;
		File	paramFile = null;

		// optional
		boolean	isPaired = false;
		boolean pairSpecified = false;
		String	fragMethod = null;
		boolean methodSpecified = false;
		File	outputFile = null;
		File	profFile = null;
		File	prmFile = null;
		float	specProbThreshold = 1e-9f;
		float	delta = 0.03f;
		float	profProbThreshold = 0.3f;
		Enzyme	enzyme = Enzyme.TRYPSIN;
		boolean isFixedModSpecified = false;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length)
				printUsageAndExit("Illegal parameter: " + argv[i]);
			if(argv[i].equalsIgnoreCase("-i"))
			{
				specFile = new File(argv[i+1]);
				if(!specFile.exists())
				{
					printUsageAndExit(specFile + " doesn't exist.");
				}
				String specFileName = specFile.getName();
				String ext = specFileName.substring(specFileName.lastIndexOf('.')+1);
				if(ext.equalsIgnoreCase("mzxml"))
					specFileFormat = SpecFileFormat.MZXML;
				else if(ext.equalsIgnoreCase("mgf"))
					specFileFormat = SpecFileFormat.MGF;
//				else if(ext.equalsIgnoreCase("pkl"))
//					specFileFormat = SpecFileFormat.PKL;
//				else if(ext.equalsIgnoreCase("ms2"))
//					specFileFormat = SpecFileFormat.MS2;
				if(specFileFormat == null)
					printUsageAndExit("Illegal file format: " + specFileName);
			}
			else if(argv[i].equalsIgnoreCase("-param"))
			{
				paramFile = new File(argv[i+1]);
				if(!paramFile.exists())
					printUsageAndExit(paramFile + " doesn't exist.");
			}
			else if(argv[i].equalsIgnoreCase("-pair"))	// Fragmentation method
			{
				if(methodSpecified)
					printUsageAndExit("Illegal parameter: at most one of -m and -pair can be specified!");
				if(argv[i+1].equalsIgnoreCase("0"))
					isPaired = false;
				else if(argv[i+1].equalsIgnoreCase("1"))
					isPaired = true;
				else
					printUsageAndExit("Illegal parameter: -pair " + argv[i+1]);
				pairSpecified = true;
			}			
			else if(argv[i].equalsIgnoreCase("-m"))
			{
				if(pairSpecified)
					printUsageAndExit("Illegal parameter: at most one of -m and -pair can be specified!");
					
				fragMethod = argv[i+1];
				methodSpecified = true;
			}
			else if(argv[i].equalsIgnoreCase("-gp"))
			{
				outputFile = new File(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-prof"))
			{
				profFile = new File(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-prm"))
			{
				prmFile = new File(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-p"))
			{
				try {
					specProbThreshold = Float.parseFloat(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal specProbThreshold: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-delta"))
			{
				try {
					delta = Float.parseFloat(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal deltaScoreForGPTemplate: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-prob"))
			{
				try {
					profProbThreshold = Float.parseFloat(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal profProbThreshold: " + argv[i+1]);
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
			else if(argv[i].equalsIgnoreCase("-fixMod"))
			{
				// 0: No mod, 1: Carbamidomethyl C, 2: Carboxymethyl C
				isFixedModSpecified = true;
				if(argv[i+1].equalsIgnoreCase("0"))
					aaSet = AminoAcidSet.getStandardAminoAcidSet();
				else if(argv[i+1].equalsIgnoreCase("1"))
					aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
				else if(argv[i+1].equalsIgnoreCase("2"))
					aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarboxymethylatedCys();
				else
					printUsageAndExit("Illigal -fixMod parameter: " + argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-aaSet") && !isFixedModSpecified)
			{
				File aaSetFile = new File(argv[i+1]);
				if(!aaSetFile.exists())
				{
					printUsageAndExit(aaSetFile + " doesn't exist.");
				}
				aaSet = AminoAcidSet.getAminoAcidSet(aaSetFile.getPath());
			}
			else
				printUsageAndExit("Illegal parameter!");
		}		
		
		if(specFile == null)
			printUsageAndExit("specFileName is not specified!");
		
		runMSProfile(specFile, specFileFormat, outputFile, profFile, prmFile, enzyme, fragMethod, isPaired,
				aaSet, specProbThreshold, delta, profProbThreshold, paramFile);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.err.println(message);
		System.out.println("MSProfile 07/08/2011");
		System.out.println("usage: java -Xmx2000M -jar MSProfile.jar \n"
				+ "\t-i SpecFileName (*.mzXML)\n"
				+ "\t[-gp GappedPeptideOutputFileName] (Default: stdout)\n"
				+ "\t[-prof ProfileOutputFileName] (Default: no output)\n"
				+ "\t[-prm PrmScoreOutputFileName] (Default: no output)\n"
				+ "\t[-m FragMethod] (if specified, FragMethod will be exclusively considered, e.g. -m CID)\n"
				+ "\t[-pair 0/1 ] (0: not paired (default), 1: paired)\n"
				+ "\t[-e Enzyme 0/1/2/3/4/5/6/7] (0: No enzyme, 1: Trypsin (default), 2: Chymotrypsin, 3: LysC, 4: LysN, 5: GluC, 6: ArgC, 7: AspN)\n"
				+ "\t[-fixMod 0/1/2] (0: NoCysteineProtection, 1: CarbamidomethyC (default), 2: CarboxymethylC)\n"
				+ "\t[-aaSet AASetFileName (default: standard amino acids)]\n"
				+ "\t[-p SpecProbThreshold] (Default: 1e-9)\n"
				+ "\t[-delta DeltaScoreForGPTemplate] (Default: 0.03)\n"
				+ "\t[-prob ProfProbThreshold] (Default: 0.3)\n"
				+ "\t[-param ScoringParamFile]\n"
				);
		System.exit(-1);
	}

	public static void runMSProfile(File specFile, SpecFileFormat format, File gpOutFile, File profOutFile, File prmOutFile,
			Enzyme enzyme, String fragMethod, boolean isPaired, AminoAcidSet aaSet, 
			float specProbThreshold, float delta, float profProbThreshold, File paramFile)
	{
		// Currently, mzXML and mgf are supported.
		Iterator<Spectrum> iterator = null;
		try {
			if(format == SpecFileFormat.MZXML)
				iterator = new MzXMLSpectraIterator(specFile.getPath());
			else if(format == SpecFileFormat.MGF)
				iterator = new SpectraIterator(specFile.getPath(), new MgfSpectrumParser());
//			else if(format == SpecFileFormat.PKL)
//				iterator = new SpectraIterator(specFile.getPath(), new PklSpectrumParser());
//			else if(format == SpecFileFormat.MS2)
//				iterator = new SpectraIterator(specFile.getPath(), new MS2SpectrumParser());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		PrintStream out = null;
		if(gpOutFile != null)
			try {
				out = new PrintStream(new BufferedOutputStream(new FileOutputStream(gpOutFile)));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		else
			out = System.out;
		
		PrintStream profOut = null;
		if(profOutFile != null)
		{
			try {
				profOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(profOutFile)));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}

		PrintStream prmOut = null;
		if(prmOutFile != null)
		{
			try {
				prmOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(prmOutFile)));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		NewRankScorer customScorer = null;
		if(paramFile != null)
			customScorer = new NewRankScorer(paramFile.getPath());
		
		// for pair
		int prevScanNum = 0;
		float prevPrecursorMz = 0;
		int prevCharge = 0;
		ActivationMethod prevMethod = null;
		ScoredSpectrum<NominalMass> cachedScoredSpec = null;
		NominalMassFactory nominalMassFactory = new NominalMassFactory(aaSet, enzyme, 50);

		out.println("SpecFileName\tScanNum\tFragmentation\tPrecursorMz\tCharge\tMSGFScore\tGappedPeptide");
		while(iterator.hasNext())
		{
			Spectrum spec = iterator.next();
			
			// CID is the default activation method  
			if(spec.getActivationMethod() == null)
			{
				spec.setActivationMethod(ActivationMethod.CID);
			}
			// spectrum file must include "scan number" field
			if(spec.getScanNum() <= 0)
			{
				System.err.println("No scan number info in the spectrum!");
				System.exit(-1);
			}
			// if fragMethod is specified, spectrum must have the same fragmentation method
			if(fragMethod != null && !spec.getActivationMethod().getName().equalsIgnoreCase(fragMethod))
				continue;
			
			int nominalPepMass = Math.round((spec.getParentMass()-(float)Composition.H2O)*Constants.INTEGER_MASS_SCALER);

			int scanNum = spec.getScanNum();
			float precursorMz = spec.getPrecursorPeak().getMz();
			boolean paired;
			if(isPaired && scanNum == prevScanNum+1 && precursorMz == prevPrecursorMz && spec.getCharge() == prevCharge && spec.getActivationMethod() != prevMethod)
			{
				paired = true;
			}
			else
			{
				paired = false;
				prevScanNum = scanNum;
				prevMethod = spec.getActivationMethod();
				prevPrecursorMz = precursorMz;
				prevCharge = spec.getCharge();
			}
			
			NewRankScorer scorer;
			if(customScorer != null)
				scorer = customScorer;
			else
				scorer = NewScorerFactory.get(spec.getActivationMethod(), enzyme);
			ScoredSpectrum<NominalMass> curScoredSpec = scorer.getScoredSpectrum(spec);
			
			ScoredSpectrum<NominalMass> scoredSpec = null;		
			String scanNumStr = null;
			String methodStr = null;
			if(!isPaired)
			{
				scoredSpec = curScoredSpec;
				scanNumStr = String.valueOf(spec.getScanNum());
				methodStr = spec.getActivationMethod().getName();
			}
			else
			{
				if(paired)
				{
					scanNumStr = ""+(spec.getScanNum()-1)+"-"+spec.getScanNum();
					ArrayList<ScoredSpectrum<NominalMass>> scoredSpecList = new ArrayList<ScoredSpectrum<NominalMass>>();
					scoredSpecList.add(cachedScoredSpec);
					scoredSpecList.add(curScoredSpec);
					scoredSpec = new ScoredSpectrumSum<NominalMass>(scoredSpecList);
					methodStr = prevMethod.getName()+"/"+spec.getActivationMethod().getName();
				}
				else
				{
					cachedScoredSpec = curScoredSpec;
					continue;
				}
			}
			
			// Computation / Print
			AminoAcidGraph graph = new AminoAcidGraph(nominalMassFactory, spec.getParentMass(), scoredSpec);
			GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(graph).enzyme(enzyme);
			gf.computeGeneratingFunction();
			
			ProfileGF<NominalMass> profGf = new ProfileGF<NominalMass>(gf);
			profGf.computeProfile(specProbThreshold);
			
			Profile<NominalMass> profile = profGf.getSpectralProfile();
			Sequence<NominalMass> gappedPeptide = profGf.getGappedPeptideWithNominalMasses(1f-delta, profProbThreshold);
			
			// determine de novo sequence, find the peptide 1) following enzyme rule and 2) having minimum error
//			String deNovoSeq = null;
//			float minError = 100f;
//			ArrayList<String> optimalScoringPeptides;
//			if(gf.getNumEqualOrBetterPeptides(gf.getMaxScore()-1) < 100)
//				optimalScoringPeptides = gf.getReconstructions(gf.getMaxScore()-1);		
//			else
//			{
//				optimalScoringPeptides = new ArrayList<String>();
//				optimalScoringPeptides.add(gf.getOneReconstruction(gf.getMaxScore()-1));
//			}
//			for(int i=0; i<optimalScoringPeptides.size(); i++)
//			{
//				String seq = optimalScoringPeptides.get(i);
//				Peptide p = new Peptide(seq);
//				float err = Math.abs(p.getParentMass()-spec.getParentMass());
//				if(enzyme == null || !enzyme.isCleaved(p))
//					err += 10f;
//				if(err < minError)
//				{
//					deNovoSeq = seq;
//					minError = err;
//				}
//			}
			
			out.print(specFile.getName()+"\t"+scanNumStr+"\t"+methodStr+"\t"+spec.getPrecursorPeak().getMz()+"\t"+spec.getCharge()+"\t");
			out.print("\t"+(gf.getMaxScore()-1));
//			out.print("\t"+deNovoSeq);
			if(gappedPeptide.size() > 1)
			{
				out.print("\t");
				for(int i=gappedPeptide.size()-2; i>=0; i--)
					out.print(nominalPepMass-gappedPeptide.get(i).getNominalMass()+",");
				out.print(nominalPepMass);
			}
			out.println();

			if(profOut != null)
			{
				profOut.println("BEGIN IONS");
				profOut.println("TITLE=Profile_"+scanNumStr);
				profOut.println("PEPMASS=" + spec.getPrecursorPeak().getMz());
				profOut.println("SCANS=" + scanNumStr);
				profOut.println("CHARGE="+spec.getCharge()+"+");
				Profile<NominalMass> suffixProfile = profile.toNominalMasses();
				for(int i=suffixProfile.size()-2; i>=1; i--)
				{
					ProfilePeak<NominalMass> p = suffixProfile.get(i);
					if(p.getProbability() >= MIN_PROF_PROB_REPORT_THRESHOLD)
					{
						int mass = nominalPepMass - p.getNode().getNominalMass();
						profOut.println(mass+"\t"+p.getProbability());
					}
				}
				profOut.println(nominalPepMass+"\t"+1);
				profOut.println("END IONS");
				profOut.println();
			}
			
			if(prmOut != null)
			{
				prmOut.println("BEGIN IONS");
				prmOut.println("TITLE=PRM_"+scanNumStr);
				prmOut.println("PEPMASS=" + spec.getPrecursorPeak().getMz());
				prmOut.println("SCANS=" + scanNumStr);
				prmOut.println("CHARGE="+spec.getCharge()+"+");
				ArrayList<NominalMass> srmMassList = gf.getGraph().getIntermediateNodeList();
				Collections.sort(srmMassList, Collections.reverseOrder());
				for(NominalMass m : srmMassList)
				{
					int mass = nominalPepMass - m.getNominalMass();
					int score = gf.getGraph().getNodeScore(m);
					if(score > 100000)
						score = -10;
					prmOut.println(mass+"\t"+score);
				}
				prmOut.println("END IONS");
				prmOut.println();
			}
		}	
		if(out != null)
			out.close();
		if(profOut != null)
			profOut.close();
		if(prmOut != null)
			prmOut.close();
	}
}
