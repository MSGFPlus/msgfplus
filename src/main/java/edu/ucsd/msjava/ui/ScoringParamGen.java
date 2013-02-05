package edu.ucsd.msjava.ui;

import java.io.File;

import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.ScoringParameterGeneratorWithErrors;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.AnnotatedSpectra;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.FileFormat;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.mzml.MzMLAdapter;
import edu.ucsd.msjava.params.FileListParameter;
import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.ObjectEnumParameter;
import edu.ucsd.msjava.params.ParamManager;

public class ScoringParamGen {

	public static final int VERSION = 8831;
	public static final String DATE = "02/04/2013";
	
	public static void main(String argv[])
	{
		ParamManager paramManager = new ParamManager("ScoringParamGen", String.valueOf(VERSION), DATE,
			"java -Xmx2000M -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen");
		
		MzMLAdapter.turnOffLogs();
		paramManager.addScoringParamGenParams();
		
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
		
		// Run program
		long time = System.currentTimeMillis();
		
		paramManager.printToolInfo();
		String errorMessage = runScoringParamGen(paramManager);
		if(errorMessage != null)
		{
			System.err.println("[Error] " + errorMessage);
			System.out.println();
		}
		else
			System.out.format("ScoringParamGen complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis()-time)/(float)1000);
		
	}

	public static String runScoringParamGen(ParamManager paramManager)
	{
		File[] resultFiles = paramManager.getFiles("i");
		File specDir = paramManager.getFile("d");
		
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
		
//		File modFile = paramManager.getFile("mod");
//		AminoAcidSet aaSet = null;
//		if(modFile == null)
//			aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		else
//		{
//			String modFileName = modFile.getName();
//			String ext = modFileName.substring(modFileName.lastIndexOf('.')+1);
//			if(ext.equalsIgnoreCase("xml"))
//				aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(modFile.getPath());
//			else
//				aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
//		}
		
		AnnotatedSpectra annotatedSpec = new AnnotatedSpectra(resultFiles, specDir, aaSet);
		System.out.print("Reading training PSMs...");
		String errMsg = annotatedSpec.parse();
		if(errMsg != null)
			return errMsg;
		System.out.println("Done.");
		
		ActivationMethod activationMethod = paramManager.getActivationMethod();
		InstrumentType instType = paramManager.getInstType();
		Enzyme enzyme = paramManager.getEnzyme();
		Protocol protocol = paramManager.getProtocol();
		SpecDataType dataType = new SpecDataType(activationMethod, instType, enzyme, protocol);
		
		ScoringParameterGeneratorWithErrors.generateParameters(
				annotatedSpec.getAnnotatedSpecContainer(), 
				dataType,
				aaSet,
				new File("."),
				false,
				true);
		return null;
	}
}
