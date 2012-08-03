package edu.ucsd.msjava.msdbsearch;

import java.io.File;

import edu.ucsd.msjava.ui.MSGFDB;


public class BuildSA {
	
	public static void main(String argv[])
	{
		if(argv.length < 2 || argv.length % 2 != 0)
			printUsageAndExit("The number of parameters must be even.");
		
		File dbPath = null;
		File outputDir = null;
		int mode = 2;
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length)
				printUsageAndExit("Illegal parameters");
			else if(argv[i].equalsIgnoreCase("-d"))
			{
				dbPath = new File(argv[i+1]);
				if(!dbPath.exists())
					printUsageAndExit(argv[i+1] + " doesn't exist.");
			}
			else if(argv[i].equalsIgnoreCase("-o"))
			{
				outputDir = new File(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-tda"))
			{
				if(argv[i+1].equals("0"))
					mode = 0;
				else if(argv[i+1].equals("1"))
					mode = 1;
				else if(argv[i+1].equals("2"))
					mode = 2;
				else
					printUsageAndExit("Illegal parameter: -tda " + argv[i+1]);
			}
		}
		if(dbPath == null)
			printUsageAndExit("Database must be specified!");
		
		buildSA(dbPath, outputDir, mode);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.out.println("Error: " + message);
		System.out.print("Usage: java -Xmx3500M BuildSA\n" +
				"\t-d DatabaseFile (*.fasta or *.fa)\n" +
//				"\t-o OutputDir\n" +
				"\t[-tda 0/1/2] (0: Target database only, 1: Concatenated target-decoy database only, 2: All (Default))\n");
		System.exit(-1);
	}
	
	public static void buildSA(File dbPath, File outputDir, int mode)
	{
		if(dbPath.isDirectory())
		{
			for(File f : dbPath.listFiles())
			{
				if(!f.getName().endsWith(".fasta") && !f.getName().endsWith(".fa"))
					continue;
				buildSAFiles(f, outputDir, mode);		
			}
		}
		else
		{
			if(dbPath.getName().endsWith(".fasta") || dbPath.getName().endsWith(".fa"))
			{
				buildSAFiles(dbPath, outputDir, mode);
			}
		}
		System.out.println("Done");
	}
	
	// mode => 0: target only, 1: target-decoy only, 2: both 
	public static void buildSAFiles(File databaseFile, File outputDir, int mode)
	{
		if(outputDir == null)
		{
			outputDir = databaseFile.getAbsoluteFile().getParentFile();
		}
		
		if(!outputDir.exists())
			outputDir.mkdir();
		
		String dbFileName = databaseFile.getName(); 
		
		// decoy
		if(mode == 1 || mode == 2)
		{
			String concatDBFileName = dbFileName.substring(0, dbFileName.lastIndexOf('.')) + MSGFDB.DECOY_DB_EXTENSION;
			File concatTargetDecoyDBFile = new File(outputDir.getPath()+File.separator+concatDBFileName);
			if(!concatTargetDecoyDBFile.exists())
			{
				System.out.println("Creating " + concatDBFileName + ".");
				if(ReverseDB.reverseDB(databaseFile.getPath(), concatTargetDecoyDBFile.getPath(), true, "XXX") == false)
				{
					System.err.println("Cannot create decoy database file!");
					System.exit(-1);
				}
			}
			System.out.println("Building suffix array: " + concatTargetDecoyDBFile.getPath());
			CompactFastaSequence tdaSequence = new CompactFastaSequence(concatTargetDecoyDBFile.getPath());
			float ratioUniqueProteins = tdaSequence.getRatioUniqueProteins();
			if(ratioUniqueProteins < 0.5f)
			{
				System.err.println("Error while indexing: " + concatTargetDecoyDBFile.getName() + " (too many redundant proteins)");
				System.err.println("Error while indexing: " + concatTargetDecoyDBFile.getName() + " (too many redundant proteins)");
				System.err.println("If the database contains forward and reverse proteins, run MS-GFDB (or BuildSA) again with \"-tda 0\"");
				System.exit(-1);
			}
			new CompactSuffixArray(tdaSequence);
		}
		
		if(mode == 0 || mode == 2)
		{
			File targetDBFile = new File(outputDir.getPath()+File.separator+dbFileName);
			if(!targetDBFile.exists())
			{
				System.out.println("Creating " + targetDBFile.getName() + ".");
				ReverseDB.copyDB(databaseFile.getPath(), targetDBFile.getPath());
			}
			System.out.println("Building suffix array: " + databaseFile.getPath());
			CompactFastaSequence sequence = new CompactFastaSequence(targetDBFile.getPath());
			new CompactSuffixArray(sequence);
		}
	}
}

