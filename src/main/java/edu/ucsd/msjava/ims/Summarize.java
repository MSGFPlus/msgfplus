package edu.ucsd.msjava.ims;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.parser.BufferedLineReader;

public class Summarize {
	public static void main(String argv[]) throws Exception
	{
		summarize();
	}
	
	public static void summarize() throws Exception
	{
		File resultFile = new File("/Users/kims336/Research/Data/IMS/5milTargets/bestPeptides.tsv");
		String s;
		BufferedLineReader in = new BufferedLineReader(resultFile.getPath());
		
		String headerLine = in.readLine();
		String[] header = headerLine.split("\t");
		int specProbCol = -1;
		for(int i=0; i<header.length; i++)
		{
			if(header[i].equals("SpecProb"))
				specProbCol = i;
		}
		if(specProbCol < 0)
		{
			System.out.println("No SpecProb column");
			System.exit(-1);
		}
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream("/Users/kims336/Research/Data/IMS/5milTargets/goodPeptides.tsv")));
		out.println(headerLine);
		float threshold = 1e-10f;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#") || s.length() == 0)
				continue;
			String[] token = s.split("\t");
			if(token.length != header.length)
				continue;
			float specProb = Float.parseFloat(token[specProbCol]);
			if(specProb < threshold)
			{
				out.println(s);
			}
		}
		in.close();
		out.close();
	}
}
