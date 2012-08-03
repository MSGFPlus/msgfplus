package edu.ucsd.msjava.msgf.test;

import java.io.File;

import edu.ucsd.msjava.parser.BufferedLineReader;


public class MSGFValidation {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1)
		{
			printUsageAndExit(null);
		}
		
		File resultFile = new File(argv[1]);
		if(!resultFile.exists())
			printUsageAndExit("File is missing!");
		validateMSGFResult(resultFile);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.err.println(message);
		System.out.println("usage: java MSGFValidation MSGFOutputFile");
		System.exit(-1);
	}
	
	public static final String SPEC_FILE = "#SpectrumFile";
	public static final String SCAN_NUM = "Scan#";
	public static final String ANNOTATION = "Annotation";
	public static final String CHARGE = "Charge";
	public static final String SPEC_PROB = "SpecProb";
	
	private static int validateMSGFResult(File resultFile) throws Exception
	{
		BufferedLineReader in = new BufferedLineReader(resultFile.getPath());
		String header = in.readLine();
		if(header == null || header.length() == 0 || !header.startsWith("#"))
		{
			System.err.println("Header is missing!");
			return -1;
		}
		
		int specFileColumn = -1;
		int scanNumColumn = -1;
		int annotationColumn = -1;
		int chargeColumn = -1;
		int specProbColumn = -1;
		String[] headings = header.split("\t");
		for(int i=0; i<headings.length; i++)
		{
			if(headings[i].equalsIgnoreCase(SPEC_FILE))
				specFileColumn = i;
			else if(headings[i].equalsIgnoreCase(SCAN_NUM))
				scanNumColumn = i;
			else if(headings[i].equalsIgnoreCase(ANNOTATION))
				annotationColumn = i;
			else if(headings[i].equalsIgnoreCase(CHARGE))
				chargeColumn = i;
			else if(headings[i].equalsIgnoreCase(SPEC_PROB))
				specProbColumn = i;
		}
		if(specFileColumn == -1)
		{
			System.err.println(SPEC_FILE + " is missing");
			return -1;
		}
		if(scanNumColumn == -1)
		{
			System.err.println(SCAN_NUM + " is missing");
			return -1;
		}
		if(annotationColumn == -1)
		{
			System.err.println(ANNOTATION + " is missing");
			return -1;
		}
		if(chargeColumn == -1)
		{
			System.err.println(CHARGE + " is missing");
			return -1;
		}
		
		// Header is okay
		int lineNum = 0;
		String s;
		while((s=in.readLine()) != null)
		{
			lineNum++;
			if(s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length <= specFileColumn || token.length <= scanNumColumn || token.length <= annotationColumn || token.length <= chargeColumn)
			{
				System.err.println("Not enough tokens at line " + lineNum + ": ");
				return -1;
			}
			int charge = Integer.parseInt(token[chargeColumn]);
			if(charge <= 0)
			{
				System.err.println("charge is not positive at line " + lineNum + ": ");
				return -1;
			}
			float specProb = Float.parseFloat(token[specProbColumn]);
			if(specProb <= 0)
			{
				System.err.println("SpecProb is not positive at line " + lineNum + ": ");
				return -1;
			}
		}
		return 0;
	}
}
