package edu.ucsd.msjava.misc;

import java.io.File;

public class GetSearchParams {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1)
			printUsageAndExit(null);
		
		File mzidFile = new File(argv[0]);
		if(!mzidFile.exists())
		{
			printUsageAndExit("File does not exist.");
		}
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.err.println(message);
		System.out.println("Usage: java GetSearchParams MSGFPlusResultFile(*.mzid)");
		System.exit(0);
	}
			
}
