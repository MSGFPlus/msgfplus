package msdbsearch;

import java.io.File;

import suffixarray.SuffixArraySequence;
import suffixarray.SuffixArray;

public class BuildSA {
	public static void main(String argv[])
	{
		if(argv.length != 1)
			System.out.println("usage: java BuildSA [path (directory or fasta file)]");
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
				if(!f.getName().endsWith(".fasta"))
					continue;
				System.out.println(f.getName());
				SuffixArraySequence sequence = new SuffixArraySequence(f.getPath());
				new SuffixArray(sequence);
			}
		}
		else
		{
			if(path.getName().contains(".fa"))
			{
				System.out.println(path.getName());
				SuffixArraySequence sequence = new SuffixArraySequence(path.getPath());
				new SuffixArray(sequence);
			}
		}
		System.out.println("Done");
	}
}
