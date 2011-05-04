package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

public class GridIndexDB {
	public static void main(String argv[]) throws Exception
	{
		String javaPath = "/usr/bin/java";
		if(argv.length != 4 && argv.length != 5)
			printUsageAndExit("Wrong arguments!");
		
		File dir = new File(argv[0]);
		if(dir == null || !dir.exists() || !dir.isDirectory())
			printUsageAndExit("Wrong projectDir");
		
		File dbPath = new File(dir.getPath() + File.separator + "database");
		if(dbPath == null || !dbPath.exists() || !dbPath.isDirectory())
			printUsageAndExit(dir+File.separator+"database doesn't exist!");
		
		File scriptPath = new File(dir.getPath() + File.separator + "indexDBScripts");
		if(!scriptPath.exists())
			scriptPath.mkdir();
		else if(!scriptPath.isDirectory())
			printUsageAndExit(dir+File.separator+"indexDBScripts doesn't exist!");
		
		if(argv.length == 5)
			javaPath = argv[4];
		
		makeScripts(dir, dbPath, scriptPath, Integer.parseInt(argv[1]), Integer.parseInt(argv[2]), Integer.parseInt(argv[3]), javaPath);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.err.println(message);
		System.out.println("usage: java misc.GridIndexDB projectDir minMass maxMass numJobs [javaPath]");
		System.exit(-1);
	}
	
	public static void makeScripts(File dir, File dbPath, File scriptPath, int minMass, int maxMass, int numJobs, String javaPath) throws Exception
	{
		for(File dbFile : dbPath.listFiles())
		{
			if(!dbFile.getName().endsWith(".fasta"))
				continue;
			int numPartitions = maxMass - minMass + 1;
			int numMasses = numPartitions / (numJobs-1);
			int mass = minMass;
			for(int i=0; i<numJobs; i++)
			{
				int max = mass + numMasses - 1;
				if(i == numJobs-1)
					max = maxMass;
				String outputName = scriptPath.getPath()+File.separator+dbFile.getName().substring(0, dbFile.getName().lastIndexOf('.'))+"_"+i+".sh"; 
				PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputName)));
				out.println("#!/bin/bash");
				out.println("#");
				out.println("#$ -cwd");
				out.println("#$ -j y");
				out.println("#$ -S /bin/bash");
				out.println("#");
				out.println(javaPath + " -Xmx1500M -cp " + 
						dir.getPath()+File.separator+"MSGFDB.jar msdbsearch.MakeIndex "
						+ dbFile.getPath() + " " + mass + " " + max);
				mass += numMasses;
				out.close();
			}			
		}
	}
}
