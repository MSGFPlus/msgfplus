package msdbsearch;

import java.io.File;

import suffixarray.SuffixArraySequence;
import suffixarray.SuffixArray;

public class BuildSA {
	public static void main(String argv[])
	{
		if(argv.length != 1)
			System.out.println("usage: java BuildSA path(directory or fasta file)");
		else
		{
			File f = new File(argv[0]);
			if(!f.exists())
				System.out.println(argv[0] + " doesn't exist.");
			else
			{
				buildSA(f);
			}
		}
	}
	
	public static void buildSA(File path)
	{
		if(path.isDirectory())
		{
			for(File f : path.listFiles())
			{
				if(!f.getName().endsWith(".fasta") && !f.getName().endsWith(".fa"))
					continue;
				buildSAFiles(f);			}
		}
		else
		{
			if(path.getName().endsWith(".fasta") || path.getName().endsWith(".fa"))
			{
				buildSAFiles(path);
			}
		}
		System.out.println("Done");
	}
	
	public static void buildSAFiles(File databaseFile)
	{
		String dbFileName = databaseFile.getName(); 
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
		SuffixArraySequence sequence = new SuffixArraySequence(databaseFile.getPath());
		new SuffixArray(sequence);

		System.out.println("Building suffix array: " + databaseFile.getPath());
		SuffixArraySequence tdaSequence = new SuffixArraySequence(concatTargetDecoyDBFile.getPath());
		new SuffixArray(tdaSequence);
	}
}

