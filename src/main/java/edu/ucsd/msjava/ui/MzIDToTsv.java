package edu.ucsd.msjava.ui;

import java.io.File;

import edu.ucsd.msjava.msutil.FileFormat;
import edu.ucsd.msjava.mzid.MzIDParser;
import edu.ucsd.msjava.mzml.MzMLAdapter;
import edu.ucsd.msjava.params.EnumParameter;
import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.ParamManager;

public class MzIDToTsv {
	public static final String VERSION = "v9108";
	
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();

		ParamManager paramManager = new ParamManager("MzIDToTsv", VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv");
		
		FileParameter inputFileParam = new FileParameter("i", "MzIDPath", "MS-GF+ output file (*.mzid) or directory containing mzid files");  
//		inputFileParam.addFileFormat(new FileFormat(".mzid"));
		inputFileParam.fileMustExist();
//		inputFileParam.mustBeAFile();
		paramManager.addParameter(inputFileParam);
		
		FileParameter outputFileParam = new FileParameter("o", "TSVFile", "TSV output file (*.tsv) (Default: MzIDFileName.tsv)");  
		outputFileParam.addFileFormat(new FileFormat(".tsv"));
		outputFileParam.setAsOptional();
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
		File mzIDPath = paramManager.getFile("i");
		
		boolean showQValue = paramManager.getIntValue("showQValue") == 1 ? true : false;
		boolean showDecoy = paramManager.getIntValue("showDecoy") == 1 ? true : false;
		boolean unroll = paramManager.getIntValue("unroll") == 1 ? true : false;
		
		if(mzIDPath.isDirectory())
		{
			for(File f : mzIDPath.listFiles())
			{
				if(f.getName().endsWith(".mzid"))
				{
					File tsvFile = new File(f.getPath().substring(0, f.getPath().lastIndexOf('.'))+".tsv");
					System.out.println("Converting " + f.getName() + " into " + tsvFile.getName());
					MzIDParser parser = new MzIDParser(f, showDecoy, !showQValue, unroll);
					parser.writeToTSVFile(tsvFile);
				}
			}
		}
		else
		{
			if(!mzIDPath.getName().endsWith(".mzid"))
			{
				return "Illegal file format: " + mzIDPath.getName();
			}
			// output tsv file
			File tsvFile = paramManager.getFile("o");
			if(tsvFile == null)
				tsvFile = new File(mzIDPath.getPath().substring(0, mzIDPath.getPath().lastIndexOf('.'))+".tsv");
			System.out.println("Converting " + mzIDPath.getName() + " into " + tsvFile.getName());
			MzIDParser parser = new MzIDParser(mzIDPath, showDecoy, !showQValue, unroll);
			parser.writeToTSVFile(tsvFile);
		}
		
		return null;
		
	}
}
