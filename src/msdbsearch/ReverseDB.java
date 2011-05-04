package msdbsearch;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;

public class ReverseDB {

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
		reverseDB(argv[0], argv[1]);
		
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("usage: java ReverseDB input(fasta) output(fasta)");
		System.exit(0);
	}
	
	public static void reverseDB(String inFileName, String outFileName) throws Exception
	{
		BufferedReader in = new BufferedReader(new FileReader(inFileName));
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFileName)));
		String s;
		StringBuffer protein = null;
		String annotation = null;
		while((s = in.readLine()) != null)
		{
			if(s.startsWith(">"))	// start of a protein
			{
				if(annotation != null)
				{
					StringBuffer rev = new StringBuffer();
					for(int i=protein.length()-1; i>=0; i--)
						rev.append(protein.charAt(i));
					out.println(">REV_" + annotation);
					out.println(rev.toString().trim());
				}
				annotation = s.substring(1);
				protein = new StringBuffer();
			}
			else
				protein.append(s);
		}
		if(protein != null && annotation != null)
		{
			StringBuffer rev = new StringBuffer();
			for(int i=protein.length()-1; i>=0; i--)
				rev.append(protein.charAt(i));
			out.println(">REV_" + annotation);
			out.println(rev.toString().trim());
		}
		in.close();
		out.close();
	}
}
