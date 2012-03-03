package misc;

import java.io.File;
import java.util.HashSet;

import msutil.DBFileFormat;
import msutil.FileFormat;
import params.*;
import parser.BufferedLineReader;

public class PepIdxToFasta {
	public static void main(String argv[]) throws Exception {
		ParamManager paramManager = new ParamManager("PepIdxToFasta", "1", "02/01/2012", "java -Xmx2000M -cp MSGFDB.jar misc.PepIdxToFasta");
		FileParameter sourceFileParam = new FileParameter("s", "*.pepidx", "pepidx file name");
		sourceFileParam.fileMustExist();
		sourceFileParam.addFileFormat(new FileFormat(".pepidx"));
		paramManager.addParameter(sourceFileParam);
		
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
		
		// Running MS-GFDB
//		paramManager.printToolInfo();
		String errorMessage = convertToFasta(paramManager);
		if(errorMessage != null)
		{
			System.err.println("[Error] " + errorMessage);
			System.out.println();
		}
		else
			System.out.format("PepIdxToFasta complete.");
		
	}
	
	public static String convertToFasta(ParamManager paramManager) throws Exception
	{
		File source = paramManager.getFile("s");
		String prevPep = "";
		BufferedLineReader in = new BufferedLineReader(source.getPath());
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#") || s.length() == 0)
				continue;
			String[] token = s.split("\t");
			if(token.length != 3)
				continue;
			String pep = token[0];
			if(pep.equals(prevPep))
				continue;
			
			prevPep = pep;
			System.out.println(">Fwd:"+pep.length());
			System.out.println(pep);
		}
		return null;
	}
	
}
