package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.parser.BufferedLineReader;

public class PEMMRProcessor {
	public static void main(String argv[]) throws Exception
	{
		process();
		System.out.println("Done");
	}
	
	public static void process() throws Exception
	{
		String resultFileName = "/Users/kims336/Research/Data/PEMMR/iTRAQ_N33T34_10ug_100cm_300min_C2_061213_MX_PEMMR.tsv";
		String outputFileName = "/Users/kims336/Research/Data/PEMMR/iTRAQ_N33T34_10ug_100cm_300min_C2_061213_MX_PEMMR_UMCID.tsv";
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
		
		BufferedLineReader in = new BufferedLineReader(resultFileName);
		
		String[] header = in.readLine().split("\t");
		int titleIndex = -1;
		int index = -1;
		for(String h : header)
		{
			++index;
			if(index != 0)
				out.print("\t");
			if(h.equalsIgnoreCase("Title"))
			{
				titleIndex = index; 
				out.print("UMC\t");
			}
			out.print(h);
		}
		out.println();
		
		if(titleIndex == -1)
		{
			System.out.println("No Title!");
			System.exit(-1);
		}
		
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.length() == 0 || s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length < 15)
				continue;
			
			for(int i=0; i<token.length; i++)
			{
				if(i != 0)
					out.print('\t');
				if(i == titleIndex)
				{
					String title = token[i];
					int umcId = -1;
					try {
//						String umcIdStr = title.split("\\|")[0];
						umcId = Integer.parseInt(title.split("\\|")[0]);
					}
					catch (Exception e)
					{
						System.err.println("Error while parsing Title: " + title);
					}
					out.print(umcId+"\t");
				}
				out.print(token[i]);
			}
			out.println();
		}
		
		in.close();
		out.close();
	}
}
