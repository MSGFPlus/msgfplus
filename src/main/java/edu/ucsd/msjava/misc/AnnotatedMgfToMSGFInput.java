package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.msutil.SpecFileFormat;
import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.parser.BufferedLineReader;

public class AnnotatedMgfToMSGFInput {
	public static final int VERSION = 7556;
	public static final String DATE = "04/05/2012";
	
	public static void main(String argv[]) throws Exception
	{
		ParamManager paramManager = new ParamManager("AnnotatedMgfToMSGFInput", String.valueOf(VERSION), DATE,
			"java -Xmx2000M -cp MSGFDB.jar misc.AnnotatedMgfToMSGFInput");

		FileParameter specParam = new FileParameter("i", "MSGFile", "Annotated mgf file (*.mgf)");
		specParam.addFileFormat(SpecFileFormat.MGF);
		specParam.fileMustExist();
		specParam.mustBeAFile();
		paramManager.addParameter(specParam);
		
		FileParameter outputFileParam = new FileParameter("o", "OutputFile", "MS-GF input file");
		paramManager.addParameter(outputFileParam);

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
		String errorMessage = convert(paramManager);
		if(errorMessage != null)
		{
			System.err.println("[Error] " + errorMessage);
			System.out.println();
		}
		else
			System.out.format("AnnotatedMgfToMSGFInput complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis()-time)/(float)1000);
		
	}

	public static String convert(ParamManager paramManager) throws Exception
	{
		File specFile = paramManager.getFile("i");
		File outputFile = paramManager.getFile("o");
		
//		SpectraIterator itr = new SpectraIterator(specFile.getPath(), new MgfSpectrumParser());
		BufferedLineReader in = new BufferedLineReader(specFile.getPath());
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		out.println("#SpectrumFile\tSpecIndex\tAnnotation\tCharge");
		String s;
		int specIndex = 0;
		int charge = 0;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("BEGIN"))
				specIndex++;
			else if(s.startsWith("END"))
			{
				charge = 0;
			}
			else if(s.startsWith("CHARGE="))
			{
				charge = Integer.parseInt(s.substring(s.indexOf('=')+1, s.lastIndexOf('+')));
			}
			else if(s.startsWith("SEQ="))
			{
				String pepSeq = s.substring(4);
				out.println(specFile.getName()+"\t"+specIndex+"\t"+"."+pepSeq+"."+"\t"+charge);
			}
		}

		in.close();
		out.close();
		
		return null;
	}	
}
