package ui;

import java.io.File;

import msscorer.NewRankScorer;
import msscorer.ScoringParameterGeneratorWithErrors;
import msscorer.NewScorerFactory.SpecDataType;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.AnnotatedSpectra;
import msutil.Enzyme;
import msutil.FileFormat;
import msutil.InstrumentType;
import msutil.Protocol;
import params.FileListParameter;
import params.FileParameter;
import params.ObjectEnumParameter;
import params.ParamManager;

public class ScoringParamGen {

	public static final int VERSION = 7067;
	public static final String DATE = "01/19/2012";
	
	public static void main(String argv[])
	{
		ParamManager paramManager = new ParamManager("ScoringParamGen", String.valueOf(NewRankScorer.VERSION), NewRankScorer.DATE,
			"java -Xmx2000M -cp MSGFDB.jar ui.ScoringParamGen");
		
		FileListParameter resFileParam = new FileListParameter("i", "ResultPath", "MSGFDBResultFile (*.tsv) or MSGFDBResultDir");
		resFileParam.addFileFormat(new FileFormat(".tsv"));
		paramManager.addParameter(resFileParam);
		
		FileParameter specDirParam = new FileParameter("d", "SpecDir", "Path to directory containing spectrum files");
		specDirParam.mustBeADirectory();
		specDirParam.fileMustExist();
		paramManager.addParameter(specDirParam);

		// ActivationMethod
		ObjectEnumParameter<ActivationMethod> fragParam = new ObjectEnumParameter<ActivationMethod>("m", "FragmentMethodID");
		ActivationMethod[] methods = ActivationMethod.getAllRegisteredActivationMethods();
		for(int i=1; i<methods.length ;i++)
		{
			ActivationMethod m = methods[i];
			if(m != ActivationMethod.FUSION)
				fragParam.registerObject(m);
		}
		paramManager.addParameter(fragParam);
		
		// Instrument type
		paramManager.addInstTypeParam(null);
		
		// Enzyme
		ObjectEnumParameter<Enzyme> enzParam = new ObjectEnumParameter<Enzyme>("e", "EnzymeID");
		Enzyme[] allEnzymes = Enzyme.getAllRegisteredEnzymes();
		for(int i=1; i<allEnzymes.length; i++)
		{
			Enzyme e = allEnzymes[i];
			enzParam.registerObject(e);
		}
		paramManager.addParameter(enzParam);
		
		// Protocol
		paramManager.addProtocolParam();
		
		paramManager.addModFileParam();

//		StringParameter nlParam = new StringParameter("nl", "NeutralLosses", "Comma separated neutral losses to consider. Specify compositions or masses");
//		nlParam.setAdditionalDescription("E.g. '-nl H3PO4', '-nl 97.995,64.064'");
//		nlParam.defaultValue(null);
//		paramManager.addParameter(nlParam);
		
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
		File modFile = paramManager.getFile("mod");
		AminoAcidSet aaSet = null;
		if(modFile == null)
			aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		else
		{
			String modFileName = modFile.getName();
			String ext = modFileName.substring(modFileName.lastIndexOf('.')+1);
			if(ext.equalsIgnoreCase("xml"))
				aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(modFile.getPath());
			else
				aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
		}
		
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
