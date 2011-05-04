package msdbsearch;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

public class ShuffleDB {

	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit();
		
		String ext1 = argv[0].substring(argv[0].lastIndexOf('.')+1);
		String ext2 = argv[1].substring(argv[1].lastIndexOf('.')+1);
		if(!ext1.equalsIgnoreCase("fasta") || !ext2.equalsIgnoreCase("fasta"))
		{
			System.out.println(ext1 + "," + ext2);
			printUsageAndExit();
		}
		shuffleDB(argv[0], argv[1], false);
		
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("usage: java ShuffleDB input(fasta) output(fasta)");
		System.exit(0);
	}
	
	public static void shuffleDB(String inFileName, String outFileName, boolean concatenate) throws Exception
	{
		BufferedReader in = new BufferedReader(new FileReader(inFileName));
		String s;
		ArrayList<StringBuffer> proteinList = new ArrayList<StringBuffer>();
		ArrayList<String> annotation = new ArrayList<String>();
		StringBuffer protein = null;
		
		// read
		while((s = in.readLine()) != null)
		{
			if(s.startsWith(">"))	// start of a protein
			{
				if(protein != null)
					proteinList.add(protein);
				annotation.add(s.trim());
				protein = new StringBuffer();
			}
			else
			{
				protein.append(s.trim());
			}
		}
		if(protein != null)
			proteinList.add(protein);
		
		in.close();
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFileName)));
		if(concatenate)
		{
			for(int i=0; i<proteinList.size(); i++)
			{
				out.println(annotation.get(i));
				out.println(proteinList.get(i));
			}
		}
		
		// shuffle
		int numProteins = proteinList.size();
		Random rand = new Random();
		for(int i=0; i<numProteins; i++)
		{
			for(int j=0; j<proteinList.get(i).length(); j++)
			{
				int protIdx = rand.nextInt(numProteins);
				int aaIdx = rand.nextInt(proteinList.get(protIdx).length());
				// swap
				char curChar = proteinList.get(i).charAt(j);
				char newChar = proteinList.get(protIdx).charAt(aaIdx);
				proteinList.get(i).setCharAt(j, newChar);
				proteinList.get(protIdx).setCharAt(aaIdx, curChar);
			}
		}
		
		// write
		for(int i=0; i<proteinList.size(); i++)
		{
			out.println(">SHFL_"+annotation.get(i).substring(1));
			out.println(proteinList.get(i));
		}
		out.close();
	}
}
