package ui;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Hashtable;

import msdbsearch.DBScanner;
import msgf.DeNovoGraph;
import msgf.FlexAminoAcidGraph;
import msgf.GeneratingFunction;
import msgf.NominalMass;
import msscorer.NewRankScorer;
import msscorer.NewScoredSpectrum;
import msscorer.NewScorerFactory;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.SpectraMap;
import msutil.Spectrum;
import msutil.SpectrumAccessorBySpecIndex;

import params.ParamManager;
import parser.InsPecTPSM;
import parser.InsPecTParser;
import parser.MS2SpectrumParser;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraMap;
import parser.PSMList;
import parser.PklSpectrumParser;

public class MSGF {
	public static final String VERSION = "6432";
	public static final String RELEASE_DATE = "09/08/2011";
	
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();

		ParamManager paramManager = new ParamManager("MSGF", MSGFDB.VERSION, MSGFDB.RELEASE_DATE, "java -Xmx2000M -cp MSGFDB.jar ui.MSGF");
		paramManager.addMSGFParams();

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
		
		// Run MS-GF
		paramManager.printToolInfo();
		String errorMessage = runMSGF(paramManager);
		if(errorMessage != null)
		{
			System.err.println("[Error] " + errorMessage);
			System.out.println();
		}
		else
			System.out.format("MS-GF complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis()-time)/(float)1000);
		
		runMSGF(paramManager);
	}
	
	public static String runMSGF(ParamManager paramManager)
	{
		String outputFileName = paramManager.getFile("o").getPath();
		
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
		
		AminoAcidSet aaSet;
		int fixModID = paramManager.getIntValue("fixMod");
		
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
		SpectrumAccessorBySpecIndex specAccessor = null;
		String prevFileName = "";
		int prevScanNum = -1;
		Spectrum spec = null;
		
		for(InsPecTPSM psm : psmList)
		{
//			if(psm.getScanNum() != 4215 || !psm.getPeptideStr().equals("KHQLSQQENIQTY"))
//				continue;
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
					spec = specAccessor.getSpectrumBySpecIndex(prevScanNum = psm.getScanNum());
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
				spec = specAccessor.getSpectrumBySpecIndex(prevScanNum = psm.getScanNum());
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
			
			NewScoredSpectrum<NominalMass> scoredSpec = scorer.getScoredSpectrum(spec);
//			GenericDeNovoGraph<NominalMass> graph = new GenericDeNovoGraph<NominalMass>(factory, (psm.getPeptide().getNominalMass()+18)/Constants.INTEGER_MASS_SCALER, pmTolerance, enzyme, scoredSpec);
			AminoAcidSet modAASet = psm.getAASet(aaSet);
			modAASet.registerEnzyme(enzyme);
			DeNovoGraph<NominalMass> graph = new FlexAminoAcidGraph(
					modAASet, 
					psm.getPeptide().getNominalMass(),
					enzyme,
					scoredSpec,
					false,
					false
					);
			
			int msgfScore = graph.getScore(psm.getAnnotation());
			GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(graph).enzyme(enzyme).doNotBacktrack().doNotCalcNumber();
			gf.setUpScoreThreshold(msgfScore);
			gf.computeGeneratingFunction();
			double specProb = gf.getSpectralProbability(msgfScore);
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
		
		System.out.format("MS-GF complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis()-time)/(float)1000);
		return null;
	}
}
