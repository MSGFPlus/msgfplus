package msdbsearch;

import java.io.File;

public class BuildSA {
	public static void main(String argv[])
	{
		if(argv.length < 2 || argv.length % 2 != 0)
			printUsageAndExit("The number of parameters must be even.");
		
		File dbPath = null;
		int mode = 2;
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length)
				printUsageAndExit("Illegal parameters");
			else if(argv[i].equalsIgnoreCase("-d"))
			{
				dbPath = new File(argv[i+1]);
				if(!dbPath.exists())
					printUsageAndExit(argv[0] + " doesn't exist.");
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
		
		buildSA(dbPath, mode);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.out.println("Error: " + message);
		System.out.print("Usage: java -Xmx3500M BuildSA\n" +
				"\t-d Database (*.fasta or *.fa)\n" +
				"\t[-tda 0/1/2] (0: Target database only, 1: Concatenated target-decoy database only, 2: All (Default))\n");
		System.exit(-1);
	}
	
	public static void buildSA(File path, int mode)
	{
		if(path.isDirectory())
		{
			for(File f : path.listFiles())
			{
				if(!f.getName().endsWith(".fasta") && !f.getName().endsWith(".fa"))
					continue;
				buildSAFiles(f, mode);		
			}
		}
		else
		{
			if(path.getName().endsWith(".fasta") || path.getName().endsWith(".fa"))
			{
				buildSAFiles(path, mode);
			}
		}
		System.out.println("Done");
	}
	
	public static void buildSAFiles(File databaseFile, int mode)
	{
		String dbFileName = databaseFile.getName(); 
		
		if(mode == 1 || mode == 2)
		{
			String concatDBFileName = dbFileName.substring(0, dbFileName.lastIndexOf('.'))+".revConcat.fasta";
			File concatTargetDecoyDBFile = new File(databaseFile.getAbsoluteFile().getParent()+File.separator+concatDBFileName);
			if(!concatTargetDecoyDBFile.exists())
			{
				System.out.println("Creating " + concatDBFileName + ".");
				if(ReverseDB.reverseDB(databaseFile.getPath(), concatTargetDecoyDBFile.getPath(), true) == false)
				{
					System.err.println("Cannot create decoy database file!");
					System.exit(-1);
				}
			}
			System.out.println("Building suffix array: " + databaseFile.getPath());
			CompactFastaSequence tdaSequence = new CompactFastaSequence(concatTargetDecoyDBFile.getPath());
			new CompactSuffixArray(tdaSequence);
		}
		
		if(mode == 0 || mode == 2)
		{
			System.out.println("Building suffix array: " + databaseFile.getPath());
			CompactFastaSequence sequence = new CompactFastaSequence(databaseFile.getPath());
			new CompactSuffixArray(sequence);
		}
	}
}

