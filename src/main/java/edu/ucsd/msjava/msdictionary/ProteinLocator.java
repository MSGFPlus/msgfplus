package edu.ucsd.msjava.msdictionary;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.suffixarray.SuffixArray;
import edu.ucsd.msjava.suffixarray.SuffixArraySequence;




public class ProteinLocator {
	public static void main(String argv[])
	{
		processNitinMSDicResult();
	}
	
	public static void processNitinMSDicResult()
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/NitinSignalPep/msDic_11_6_1800_20.txt";
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		HashSet<String> pepSet = new HashSet<String>();
		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length == 9 && token[0].startsWith("null"))
			{
				pepSet.add(token[2].substring(token[2].indexOf('.')+1, token[2].lastIndexOf('.')));
			}
		}
		
		
		Hashtable<String, ArrayList<String>> table = new Hashtable<String, ArrayList<String>>();
		
		File dir = new File(System.getProperty("user.home")+"/Research/Data/HumanGenome/translated");
		assert(dir.isDirectory());
		for(File f : dir.listFiles())
		{
			if(f.getName().endsWith("fasta"))
			{
				SuffixArray sa = new SuffixArray(new SuffixArraySequence(f.getPath(), edu.ucsd.msjava.sequences.Constants.AMINO_ACIDS_19));
				for(String pep : pepSet)
				{
					ArrayList<String> matches = table.get(pep);
					if(matches == null)
						table.put(pep, sa.getAllMatchedStrings(pep, 20));
					else
						matches.addAll(sa.getAllMatchedStrings(pep, 20));
				}
				sa = null;
			}
		}
		
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length == 9 && token[0].startsWith("null"))
			{
				String pep = token[2].substring(token[2].indexOf('.')+1, token[2].lastIndexOf('.'));
				ArrayList<String> matches = table.get(pep);
				for(String m : matches)
				{
					System.out.print(token[0]);
					for(int i=1; i<token.length; i++)
					{
						if(i != 2)
							System.out.print("\t"+token[i]);
						else
							System.out.print("\t"+m);
					}
					System.out.println();
				}
			}
			else
				System.out.println(s);
		}
		
	}
}
