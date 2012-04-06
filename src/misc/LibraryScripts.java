package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import msgf.Histogram;

import parser.BufferedLineReader;

public class LibraryScripts {
	public static void main(String argv[]) throws Exception
	{
//		convert();
		analyzeLibraryPSMs();
	}

	public static void analyzeLibraryPSMs() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/SpecLib/MSGFOut_Human.tsv";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		Histogram<Integer> hist = new Histogram<Integer>();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#") || s.length() == 0)
				continue;
			String[] token = s.split("\t");
			if(token[4].startsWith("N"))
				continue;
			String annotation = token[2];
			int pepLength = 0;
			for(int i=0; i<annotation.length(); i++)
				if(Character.isLetter(annotation.charAt(i)))
					pepLength++;
			
			int charge = Integer.parseInt(token[3]);
			float specProb = Float.parseFloat(token[4]);
			float specProbScore = -(float)Math.log10(specProb);
			hist.add(Math.round(specProbScore));
		}
		hist.printSortedRatio();
	}
	
	public static void convert() throws Exception
	{
		File inputFile = new File(System.getProperty("user.home")+"/Research/Data/SpecLib/human_target.mgf");
		File outputFile = new File(System.getProperty("user.home")+"/Research/Data/SpecLib/human_target_annotated.mgf");
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		String s;
		BufferedLineReader in = new BufferedLineReader(inputFile.getPath());
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SEQ="))
			{
				s = s.replaceAll("\\(C,57\\.02146\\)", "C");
				s = s.replaceAll("\\(C,39\\.99\\)", "\\(C,-17\\.026549\\)");
				StringBuffer buf = new StringBuffer();
				String[] token = s.split("[(,)]");
				for(int i=0; i<token.length; i++)
				{
					if(token[i].length() == 0)
						continue;
					if(Character.isDigit(token[i].charAt(0)))
						buf.append("+");
					buf.append(token[i]);
				}
				out.println(buf.toString());
			}
			else if(s.startsWith("PRECURSOR="))
			{
				out.println("PEPMASS="+s.substring("PRECURSOR=".length()));
			}
			else
				out.println(s);
		}
		in.close();
		out.close();
		System.out.println("Done");
	}
}
