package ui;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Hashtable;

import msdbsearch.DBScanner;
import msgf.GeneratingFunction;
import msgf.GenericDeNovoGraph;
import msgf.IntMassFactory;
import msgf.ScoredSpectrum;
import msgf.Tolerance;
import msgf.IntMassFactory.IntMass;
import msscorer.NewRankScorer;
import msscorer.NewScoredSpectrum;
import msscorer.NewScorerFactory;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.CompositionFactory;
import msutil.Constants;
import msutil.Enzyme;
import msutil.SpectraMap;
import msutil.Spectrum;
import msutil.SpectrumAccessorByScanNum;

import parser.InsPecTPSM;
import parser.InsPecTParser;
import parser.MS2SpectrumParser;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraMap;
import parser.PSMList;
import parser.PklSpectrumParser;

public class MSGF {
	public static void main(String argv[])
	{
		if(argv.length < 2 || argv.length % 2 != 0)
		{
			System.out.println("Illegal parameter: #parameters must be even (" + argv.length + ")");
			for(int i=0; i<argv.length; i++)
				System.out.println((i+1)+": " + argv[i]);
			printUsageAndExit(null);
		}

		File 	specDir = null;
		File 	resultFile 	= null;
		File 	databaseFile 	= null;
		File	outputFile = null;
		File	paramFile = null;
		float	specProbThreshold = 1f;
		Enzyme	enzyme = Enzyme.TRYPSIN;
		ActivationMethod activationMethod = ActivationMethod.CID;
		
		int program = -1;	// 0: InsPecT, 1: AnnotatedMgf, 2: listFile
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		boolean onePerSpec = false;
		boolean isFixedModSpecified = false;
		boolean isAASetSpecified = false;
		boolean addMSGFColumn = false;
		
		float nTermFixedMod = 0;
		float cTermFixedMod = 0;
		
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length)
				printUsageAndExit("Illegal parameter: " + argv[i]);
			if(argv[i].equalsIgnoreCase("-i") || argv[i].equalsIgnoreCase("-inspect"))
			{
				program = 0;
				resultFile = new File(argv[i+1]);
				if(!resultFile.exists())
				{
					printUsageAndExit(resultFile + " doesn't exist.");
				}
			}
			else if(argv[i].equalsIgnoreCase("-d"))
			{
				specDir = new File(argv[i+1]);
				if(!specDir.isDirectory())
				{
					printUsageAndExit(argv[i+1]+" is not a directory.");
				}
			}
			else if(argv[i].equalsIgnoreCase("-param"))
			{
				paramFile = new File(argv[i+1]);
				if(!paramFile.exists())
				{
					printUsageAndExit(paramFile + " doesn't exist.");
				}
			}
			else if(argv[i].equalsIgnoreCase("-db"))
			{
				databaseFile = new File(argv[i+1]);
				if(!databaseFile.exists())
				{
					printUsageAndExit(argv[i+1]+" doesn't exist.");
				}
				String dbFileName = databaseFile.getName();
				String dbFileExt = dbFileName.substring(dbFileName.lastIndexOf('.')+1);
				if(!dbFileExt.equalsIgnoreCase("fasta"))
					printUsageAndExit(argv[i+1]+" should be a fasta file!");
			}
			else if(argv[i].equalsIgnoreCase("-o"))
			{
				outputFile = new File(argv[i+1]);
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
			else if(argv[i].equalsIgnoreCase("-m"))	// Fragmentation method
			{
				// 0: CID (default), 1: ETD
				if(argv[i+1].equalsIgnoreCase("0"))
					activationMethod = ActivationMethod.CID;
				else if(argv[i+1].equalsIgnoreCase("1"))
					activationMethod = ActivationMethod.ETD;
				else if(argv[i+1].equalsIgnoreCase("2"))
					activationMethod = ActivationMethod.HCD;
				else
					printUsageAndExit("Wrong fragmentation method: " + argv[i+1]);
			}			
			else if(argv[i].equalsIgnoreCase("-x"))
			{
				if(argv[i+1].equalsIgnoreCase("0"))
					onePerSpec = false;
				else if(argv[i+1].equalsIgnoreCase("1"))
					onePerSpec = true;
				else
					printUsageAndExit("Illigal parameter: " + argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-fixMod"))
			{
				if(isAASetSpecified)
					printUsageAndExit("-aaSet and -fixMod cannot be specified together!");
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
			else if(argv[i].equalsIgnoreCase("-nFixedMod"))
			{
				try {
					nTermFixedMod = Float.parseFloat(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal nFixedMod: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-cFixedMod"))
			{
				try {
					cTermFixedMod = Float.parseFloat(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal cFixedMod: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-aaSet"))
			{
				if(isFixedModSpecified)
					printUsageAndExit("-aaSet and -fixMod cannot be specified together!");
				File aaSetFile = new File(argv[i+1]);
				if(!aaSetFile.exists())
				{
					printUsageAndExit(aaSetFile + " doesn't exist.");
				}
				aaSet = AminoAcidSet.getAminoAcidSet(aaSetFile.getPath());
				isAASetSpecified = true;
			}
			else if(argv[i].equalsIgnoreCase("-addScore"))
			{
				if(argv[i+1].equalsIgnoreCase("0"))
					addMSGFColumn = false;
				else if(argv[i+1].equalsIgnoreCase("1"))
					addMSGFColumn = true;
				else
					printUsageAndExit("Illigal parameter: " + argv[i+1]);
			}
			else
			{
				printUsageAndExit("Illegal parameter!: " + argv[i]);
			}
		}
		
		if(databaseFile != null)
			DBScanner.setAminoAcidProbabilities(databaseFile.getPath(), aaSet);
		
		if(program < 0 || resultFile == null)
			printUsageAndExit("ResultFile is missing!");
		if(program == 0 && specDir == null)
			printUsageAndExit("SpecDir is missing!");
			
		runMSGF(resultFile, specDir, outputFile, aaSet, nTermFixedMod, cTermFixedMod, enzyme, activationMethod, onePerSpec, program, specProbThreshold, paramFile, addMSGFColumn);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.err.println(message);
		System.out.println("\nMSGFv2 v20101217");
		System.out.println("usage: java -Xmx2000M -jar MSGF.jar \n" +
				"\t-i ResultFile\n" // or -mascot mascotResult (*.dat) or -list idListFile\n"				 
				+ "\t-d SpecDir\n"// or -mgf annotatedMgfFile\n"
				+ "\t[-o OutputFileName] (Default: stdout)\n"
				+ "\t[-db DatabaseFileName] (for computing AA probabilities, if not specified, 1/20 is used for all AAs)\n"
				+ "\t[-m FragmentationMethod 0/1] (0: CID (default), 1: ETD, 2: HCD)\n"
				+ "\t[-e Enzyme 0/1/2/3/4/5/6/7] (0: No enzyme, 1: Trypsin (default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N)\n"
				+ "\t[-fixMod 0/1/2] (0: NoCysteineProtection, 1: CarbamidomethyC (default), 2: CarboxymethylC)\n"
//				+ "\t[-nFixedMod NTermFixedModMass] (Default: 0)\n"
//				+ "\t[-cFixedMod CTermFixedModMass] (Default: 0)\n"
				+ "\t[-aaSet AASetFileName (default: standard amino acids)]\n"
				+ "\t[-x 0/1] (0: all (default), 1: OnePerSpec)\n"
				+ "\t[-p SpecProbThreshold] (Default: 1)\n"
				+ "\t[-param ScoringParamFile]\n" 
				+ "\t[-addScore 0/1] (0: don't add MSGFScore (default), 1: add MSGFScore)\n"
//				+ "\t[-scoreFilter score]"
				);
		System.exit(-1);
	}
	
	public static void runMSGF(File resultFile, File specDir, File outputFileName, AminoAcidSet aaSet, float nTermFixedMod, float cTermFixedMod,
			Enzyme enzyme, ActivationMethod activationMethod, 
			boolean onePerSpec, int program, float specProbThreshold, File paramFile, boolean addMSGFColumn)
	{
		PrintStream out = null;
		if(outputFileName == null)
			out = System.out;
		else
		{
			try {
				out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		NewRankScorer scorer;
		if(paramFile != null)
		{
			scorer = new NewRankScorer(paramFile.getPath());
		}
		else
		{
			scorer = NewScorerFactory.get(activationMethod, enzyme);
		}

		float rescalingFactor = Constants.INTEGER_MASS_SCALER;;
		Tolerance pmTolerance = new Tolerance(0.1f);
		
		IntMassFactory factory = new IntMassFactory(aaSet, enzyme, Constants.MAX_PEPTIDE_LENGTH, rescalingFactor, true);
		int gfTableCapacity = factory.getMassIndex(factory.getAASet().getHeaviestAA().getMass());
		
		InsPecTParser parser = new InsPecTParser(aaSet);
		parser.parse(resultFile.getPath());
		
		String header = parser.getHeader();
		if(header != null)
		{
			out.print(header+"\tSpecProb");
			if(addMSGFColumn)
				out.print("\tMSGFScore\tDeNovoScore");
			out.println();
		}
			
		ArrayList<String> keyList = null;
		Hashtable<String, Double> minSpecProb = null;
		Hashtable<String, String> bestOut = null;
		if(onePerSpec)
		{
			minSpecProb = new Hashtable<String, Double>();
			bestOut = new Hashtable<String, String>();
			keyList = new ArrayList<String>();
		}
		
		PSMList<InsPecTPSM> psmList = parser.getPSMList(); 
		if(psmList == null)
		{
			out.println("The result file is empty!");
			out.close();
			return;
		}
//		Collections.sort(psmList, new PSM.PSMSpecFileAndScanNumComparator());
		SpectrumAccessorByScanNum specAccessor = null;
		String prevFileName = "";
		int prevScanNum = -1;
		Spectrum spec = null;
		
		for(InsPecTPSM psm : psmList)
		{
			if(psm.getPeptide() == null)
			{
				out.print(psm.getInsPecTString()+"\t"+"N/A: unrecognizable annotation");
				if(addMSGFColumn)
					out.print("\t\t");
				out.println();
				continue;
			}
			String fileName = psm.getSpecFileName();
			if(fileName.equalsIgnoreCase(prevFileName) && specAccessor != null)	// same spec file name
			{
				assert(specAccessor != null);
				if(psm.getScanNum() == prevScanNum)	// same spectrum
					assert(spec != null);
				else	// different spectrum
				{
					spec = specAccessor.getSpectrumByScanNum(prevScanNum = psm.getScanNum());
					prevScanNum = psm.getScanNum();
					if(onePerSpec)
						keyList.add(psm.getSpecFileName()+":"+psm.getScanNum());
				}
			}
			else	// different spec file name
			{
				prevFileName = fileName;
				String filePrefix = fileName.substring(0, fileName.lastIndexOf('.'));
				String ext = fileName.substring(fileName.lastIndexOf('.'));
				String specFilePath = specDir.getPath()+File.separatorChar+fileName;
				if(ext.equalsIgnoreCase(".mzxml"))	// mzXML
				{
					File specFile = new File(specFilePath);
					if(!specFile.exists())
					{
						for(File f : specDir.listFiles())
						{
							if(f.getName().startsWith(filePrefix))
							{
								if(f.getName().substring(f.getName().lastIndexOf('.')).equalsIgnoreCase(".mzxml"))
									specFile = f;
							}
						}
						if(!specFile.exists())
						{
							out.print(psm.getInsPecTString()+"\t"+"N/A: spectrum file is missing");
							if(addMSGFColumn)
								out.print("\t\t");
							out.println();
							continue;
						}
					}
					specAccessor = new MzXMLSpectraMap(specFile.getPath());
				}
				else if(ext.equalsIgnoreCase(".mgf"))
				{
					specAccessor = new SpectraMap(specFilePath, new MgfSpectrumParser());
				}
				else if(ext.equalsIgnoreCase(".pkl"))
				{
					specAccessor = new SpectraMap(specFilePath, new PklSpectrumParser());
				}
				else if(ext.equalsIgnoreCase(".ms2"))
				{			
					specAccessor = new SpectraMap(specFilePath, new MS2SpectrumParser());
				}
				else
				{
					out.print(psm.getInsPecTString()+"\t"+"N/A: unrecognizable spec format");
					if(addMSGFColumn)
						out.print("\t\t");
					out.println();
					continue;
				}
				spec = specAccessor.getSpectrumByScanNum(prevScanNum = psm.getScanNum());
				if(onePerSpec)
					keyList.add(psm.getSpecFileName()+":"+psm.getScanNum());
			}
			
			if(spec == null)
			{
				out.print(psm.getInsPecTString()+"\t"+"N/A: unrecognizable spec format");
				if(addMSGFColumn)
					out.print("\t\t");
				out.println();
				continue;
			}
			
			if(psm.getPeptide() == null || psm.getPeptide().contains(null))
			{			
				out.print(psm.getInsPecTString()+"\t"+"N/A: unrecognizable identification: " + psm.getPeptideStr());
				if(addMSGFColumn)
					out.print("\t\t");
				out.println();
				continue;
			}
			
			spec.getPrecursorPeak().setCharge(psm.getCharge());
			
			float expPM = spec.getParentMass();
			float calcPM = psm.getPeptide().getParentMass() + nTermFixedMod + cTermFixedMod;
			if(Math.abs(expPM - calcPM) > 10)
			{
				out.print(psm.getInsPecTString()+"\t"+"N/A: precursor mass != peptide mass (" + expPM + " vs " + calcPM + ")");
				if(addMSGFColumn)
					out.print("\t\t");
				out.println();
				continue;
			}
			
//			if(!psm.getAnnotation().toString().equalsIgnoreCase("R.IANLNKR.Y"))
//				continue;
			NewScoredSpectrum<IntMass> scoredSpec = scorer.getScoredSpectrum(spec);
			GenericDeNovoGraph<IntMass> graph = new GenericDeNovoGraph<IntMass>(factory, (psm.getPeptide().getNominalMass()+18)/Constants.INTEGER_MASS_SCALER, pmTolerance, enzyme, scoredSpec);
			
			GeneratingFunction<IntMass> gf = new GeneratingFunction<IntMass>(graph).enzyme(enzyme).doNotBacktrack().doNotCalcNumber().gfTableCapacity(gfTableCapacity);
			
			gf.computeGeneratingFunction();
			
			double specProb = gf.getSpectralProbability(psm.getAnnotation());
			int msgfScore = gf.getScore(psm.getAnnotation());
			assert(specProb > 0): psm.getInsPecTString()+"\t"+"SpecProb is zero!";
			String output = psm.getInsPecTString()+"\t"+specProb;
			if(addMSGFColumn)
				output += "\t"+msgfScore+"\t"+(gf.getMaxScore()-1);
			if(specProb <= specProbThreshold)
			{
				if(!onePerSpec)
					out.println(output);
				else
				{
					String specKey = psm.getSpecFileName()+":"+psm.getScanNum();
					Double prevBest = minSpecProb.get(specKey);
					if(prevBest == null || specProb < prevBest)
					{
						minSpecProb.put(specKey, specProb);
						bestOut.put(specKey, output);
					}
				}
			}
//			System.exit(0);
		}
		
		if(onePerSpec)
		{
			for(String key : keyList)
				out.println(bestOut.get(key));
		}
		
		out.flush();
		out.close();	
	}
	
	public static void runMSGFCompGraph(File resultFile, File specDir, File outputFileName, AminoAcidSet aaSet, float nTermFixedMod, float cTermFixedMod,
			Enzyme enzyme, ActivationMethod activationMethod, 
			boolean onePerSpec, int program, float specProbThreshold, File paramFile, boolean addMSGFColumn)
	{
		PrintStream out = null;
		if(outputFileName == null)
			out = System.out;
		else
		{
			try {
				out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		int maxLength = 20;
		
		Tolerance pmTolerance = new Tolerance(30, true);
		CompositionFactory factory = new CompositionFactory(aaSet, enzyme, maxLength+2);
		System.out.println("Composition Graph is built!");
		
		InsPecTParser parser = new InsPecTParser(aaSet);
		parser.parse(resultFile.getPath());
		
		String header = parser.getHeader();
		if(header != null)
		{
			out.print(header+"\tSpecProb");
			if(addMSGFColumn)
				out.print("\tMSGFScore");
			out.println();
		}
		
		NewRankScorer customScorer = null;
		if(paramFile != null)
			customScorer = new NewRankScorer(paramFile.getPath());
			
		ArrayList<String> keyList = null;
		Hashtable<String, Double> minSpecProb = null;
		Hashtable<String, String> bestOut = null;
		if(onePerSpec)
		{
			minSpecProb = new Hashtable<String, Double>();
			bestOut = new Hashtable<String, String>();
			keyList = new ArrayList<String>();
		}
		
		PSMList<InsPecTPSM> psmList = parser.getPSMList(); 
		if(psmList == null)
		{
			out.println("The result file is empty!");
			out.close();
			return;
		}
//		Collections.sort(psmList, new PSM.PSMSpecFileAndScanNumComparator());
		SpectrumAccessorByScanNum specAccessor = null;
		String prevFileName = "";
		int prevScanNum = -1;
		Spectrum spec = null;
		
		
		for(InsPecTPSM psm : psmList)
		{
			if(psm.getPeptide() == null)
			{
				out.print(psm.getInsPecTString()+"\t"+"N/A: unrecognizable annotation");
				if(addMSGFColumn)
					out.print("\t");
				out.println();
				continue;
			}
			String fileName = psm.getSpecFileName();
			if(fileName.equalsIgnoreCase(prevFileName) && specAccessor != null)	// same spec file name
			{
				assert(specAccessor != null);
				if(psm.getScanNum() == prevScanNum)	// same spectrum
					assert(spec != null);
				else	// different spectrum
				{
					spec = specAccessor.getSpectrumByScanNum(prevScanNum = psm.getScanNum());
					prevScanNum = psm.getScanNum();
					if(onePerSpec)
						keyList.add(psm.getSpecFileName()+":"+psm.getScanNum());
				}
			}
			else	// different spec file name
			{
				prevFileName = fileName;
				String filePrefix = fileName.substring(0, fileName.lastIndexOf('.'));
				String ext = fileName.substring(fileName.lastIndexOf('.'));
				String specFilePath = specDir.getPath()+File.separatorChar+fileName;
				if(ext.equalsIgnoreCase(".mzxml"))	// mzXML
				{
					File specFile = new File(specFilePath);
					if(!specFile.exists())
					{
						for(File f : specDir.listFiles())
						{
							if(f.getName().startsWith(filePrefix))
							{
								if(f.getName().substring(f.getName().lastIndexOf('.')).equalsIgnoreCase(".mzxml"))
									specFile = f;
							}
						}
						if(!specFile.exists())
						{
							out.print(psm.getInsPecTString()+"\t"+"N/A: spectrum file is missing");
							if(addMSGFColumn)
								out.print("\t");
							out.println();
							continue;
						}
					}
					specAccessor = new MzXMLSpectraMap(specFile.getPath());
				}
				else if(ext.equalsIgnoreCase(".mgf"))
				{
					specAccessor = new SpectraMap(specFilePath, new MgfSpectrumParser());
				}
				else if(ext.equalsIgnoreCase(".pkl"))
				{
					specAccessor = new SpectraMap(specFilePath, new PklSpectrumParser());
				}
				else if(ext.equalsIgnoreCase(".ms2"))
				{			
					specAccessor = new SpectraMap(specFilePath, new MS2SpectrumParser());
				}
				else
				{
					out.print(psm.getInsPecTString()+"\t"+"N/A: unrecognizable spec format");
					if(addMSGFColumn)
						out.print("\t");
					out.println();
					continue;
				}
				spec = specAccessor.getSpectrumByScanNum(prevScanNum = psm.getScanNum());
				if(onePerSpec)
					keyList.add(psm.getSpecFileName()+":"+psm.getScanNum());
			}
			
			if(spec == null)
			{
				out.print(psm.getInsPecTString()+"\t"+"N/A: unrecognizable spec format");
				if(addMSGFColumn)
					out.print("\t");
				out.println();
				continue;
			}
			
			if(psm.getPeptide() == null || psm.getPeptide().contains(null))
			{
				out.print(psm.getInsPecTString()+"\t"+"N/A: unrecognizable identification: " + psm.getPeptideStr());
				if(addMSGFColumn)
					out.print("\t");
				out.println();
				continue;
			}
			
			spec.getPrecursorPeak().setCharge(psm.getCharge());
			
			float expPM = spec.getParentMass();
			float calcPM = psm.getPeptide().getParentMass() + nTermFixedMod + cTermFixedMod;
			if(Math.abs(expPM - calcPM) > 10)
			{
				out.print(psm.getInsPecTString()+"\t"+"N/A: precursor mass != peptide mass (" + expPM + " vs " + calcPM + ")");
				if(addMSGFColumn)
					out.print("\t");
				out.println();
				continue;
			}
			
			if(psm.getPeptide().size() > maxLength)
			{
				out.print(psm.getInsPecTString()+"\t"+"N/A: peptide is too long for the composition graph");
				if(addMSGFColumn)
					out.print("\t");
				out.println();
				continue;
			}
			
			
			NewRankScorer scorer;
			if(customScorer != null)
				scorer = customScorer;
			else
				scorer = NewScorerFactory.get(activationMethod, enzyme);
			NewScoredSpectrum<Composition> scoredSpec = scorer.getScoredSpectrum(spec);
			GenericDeNovoGraph<Composition> graph = new GenericDeNovoGraph<Composition>(factory, spec.getParentMass(), pmTolerance, enzyme, scoredSpec);
			GeneratingFunction<Composition> gf = new GeneratingFunction<Composition>(graph).enzyme(enzyme).doNotBacktrack().doNotCalcNumber();
			
			gf.computeGeneratingFunction();
			
			double specProb = gf.getSpectralProbability(psm.getAnnotation());
			int msgfScore = gf.getScore(psm.getAnnotation());
			
			String output = psm.getInsPecTString()+"\t"+specProb;
			if(addMSGFColumn)
				output += "\t"+msgfScore;
			if(specProb <= specProbThreshold)
			{
				if(!onePerSpec)
					out.println(output);
				else
				{
					String specKey = psm.getSpecFileName()+":"+psm.getScanNum();
					Double prevBest = minSpecProb.get(specKey);
					if(prevBest == null || specProb < prevBest)
					{
						minSpecProb.put(specKey, specProb);
						bestOut.put(specKey, output);
					}
				}
			}
		}
		
		if(onePerSpec)
		{
			for(String key : keyList)
				out.println(bestOut.get(key));
		}
		
		out.flush();
		out.close();	
	}	
}
