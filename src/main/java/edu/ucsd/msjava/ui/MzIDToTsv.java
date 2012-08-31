package edu.ucsd.msjava.ui;

import java.io.File;

import edu.ucsd.msjava.msutil.FileFormat;
import edu.ucsd.msjava.mzid.MzIDParser;
import edu.ucsd.msjava.mzml.MzMLAdapter;
import edu.ucsd.msjava.params.EnumParameter;
import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.ParamManager;

public class MzIDToTsv {
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();

		ParamManager paramManager = new ParamManager("MzIDToTsv", "8299", MSGFPlus.RELEASE_DATE, "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv");
		
		FileParameter inputFileParam = new FileParameter("i", "MzIDFile", "MS-GF+ output file (*.mzid)");  
		inputFileParam.addFileFormat(new FileFormat(".mzid"));
		inputFileParam.fileMustExist();
		inputFileParam.mustBeAFile();
		paramManager.addParameter(inputFileParam);
		
		FileParameter outputFileParam = new FileParameter("o", "TSVFile", "TSV output file (*.tsv)");  
		outputFileParam.addFileFormat(new FileFormat(".tsv"));
		paramManager.addParameter(outputFileParam);	
		
		EnumParameter showQValueParam = new EnumParameter("showQValue");
		showQValueParam.registerEntry("do not show Q-values");
		showQValueParam.registerEntry("show Q-values").setDefault();
		paramManager.addParameter(showQValueParam);

		EnumParameter showDecoyParam = new EnumParameter("showDecoy");
		showDecoyParam.registerEntry("do not show decoy PSMs").setDefault();
		showDecoyParam.registerEntry("show decoy PSMs");
		paramManager.addParameter(showDecoyParam);

		EnumParameter rank1OnlyParam = new EnumParameter("unroll");
		rank1OnlyParam.registerEntry("merge shared peptides").setDefault();
		rank1OnlyParam.registerEntry("unroll shared peptides");
		paramManager.addParameter(rank1OnlyParam);
		
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
		String errorMessage = convert(paramManager);
		if(errorMessage != null)
		{
			System.err.println("[Error] " + errorMessage);
			System.out.println();
		}
		else
			System.out.format("MzIDToTsv complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis()-time)/(float)1000);		
	}
	
	public static String convert(ParamManager paramManager)
	{
		// mzid File
		File mzIDFile = paramManager.getFile("i");
		
		// output tsv file
		File tsvFile = paramManager.getFile("o");
		
		boolean showQValue = paramManager.getIntValue("showQValue") == 1 ? true : false;
		boolean showDecoy = paramManager.getIntValue("showDecoy") == 1 ? true : false;
		boolean unroll = paramManager.getIntValue("unroll") == 1 ? true : false;
		
		MzIDParser parser = new MzIDParser(mzIDFile, showDecoy, !showQValue, unroll);
		parser.writeToTSVFile(tsvFile);
		
		return null;
		
	}
}
