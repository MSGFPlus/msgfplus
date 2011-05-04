package misc;

import java.io.File;

public class GridMakeQsub {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1 && argv.length != 2)
			printUsageAndExit("Wrong arguments!");
		
		File dir = new File(argv[0]);
		if(!dir.exists() || !dir.isDirectory())
			printUsageAndExit("Wrong directory!");
		
		String option = null;
		if(argv.length == 2)
			option = argv[1];
		makeDriver(dir, option);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.err.println(message);
		System.out.println("usage: java GridMakeQsub dir [option (e.g. -l h_vmem=2G,virtual_free=2G)]");
		System.exit(-1);
	}
	
	public static void makeDriver(File dir, String option) throws Exception
	{
		System.out.println("#!/bin/bash");
		for(File f : dir.listFiles())
		{
			if(f.getName().endsWith(".sh"))
			{
				System.out.println("qsub " +
					(option != null ? option+" " : "") +
					f.getPath());
			}
		}
	}
}
