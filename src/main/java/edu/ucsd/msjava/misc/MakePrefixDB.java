package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.msutil.DBFileFormat;
import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.IntParameter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.parser.BufferedLineReader;

public class MakePrefixDB {
	public static final int VERSION = 7711;
	public static final String DATE = "05/03/2012";
	
	public static void main(String argv[]) throws Exception
	{
		ParamManager paramManager = new ParamManager("MakePrefixDB", String.valueOf(VERSION), DATE,
			"java -Xmx2000M -cp MSGFDB.jar misc.MakePrefixDB");

		FileParameter inputDBParam = new FileParameter("i", "DatabaseFile", "Input fasta file (*.fa, *.fasta)");
		inputDBParam.addFileFormat(DBFileFormat.FASTA);
		inputDBParam.fileMustExist();
		inputDBParam.mustBeAFile();
		paramManager.addParameter(inputDBParam);

		FileParameter outputDBParam = new FileParameter("o", "DatabaseFile", "Input fasta file (*.fa, *.fasta)");
		outputDBParam.addFileFormat(DBFileFormat.FASTA);
		outputDBParam.fileMustNotExist();
		paramManager.addParameter(outputDBParam);
		
		IntParameter prefixLengthParam = new IntParameter("l", "PrefixLength", "Prefix length, Default: 30");
		prefixLengthParam.minValue(1);
		prefixLengthParam.defaultValue(30);
		paramManager.addParameter(prefixLengthParam);
		
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
			System.out.format("MakePrefixDB complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis()-time)/(float)1000);
		
	}

	public static String convert(ParamManager paramManager) throws Exception
	{
		File inputDBFile = paramManager.getFile("i");
		File outputDBFile = paramManager.getFile("o");
		int prefixLength = paramManager.getIntValue("l");
		
		BufferedLineReader in = new BufferedLineReader(inputDBFile.getPath());
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputDBFile)));
		String s;
		int remaining = 0;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith(">"))
			{
				out.println(s);
				remaining = prefixLength;
			}
			else if(remaining > 0)
			{
				if(s.length() >= remaining)
				{
					out.println(s.substring(0, remaining));
					remaining = 0;
				}
				else
				{
					out.println(s);
					remaining -= s.length();
				}
			}
		}

		in.close();
		out.close();
		
		return null;
	}	
}
