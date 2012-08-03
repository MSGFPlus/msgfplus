package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.msgf.Histogram;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.MgfSpectrumParser;


public class LibraryScripts {
	public static void main(String argv[]) throws Exception
	{
//		convert();
//		analyzeLibraryPSMs();
		makeTh();
//		extractShortPeptides();
//		cleanMgf();
	}

	public static void cleanMgf() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/Heck_DDDT/CID_IT.mgf";
		String outputFileName = "/home/sangtaekim/Test/MSGFLib/CID_IT_Idx.mgf";
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
		int specIndex = 0;
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			++specIndex;
			String title = spec.getTitle();
			spec.setTitle(String.valueOf(specIndex));
			spec.outputMgf(out);
		}
		out.close();
		System.out.println("Done");
	}
	
	public static void extractShortPeptidesSpectraST() throws Exception
	{
		File dir = new File("/home/sangtaekim/Test/MSGFLib");
		String fileName = dir.getPath()+"/SpectraST_1.tsv";
		String outputFileName = dir.getPath() + "/SpectraST_1_Long.tsv";
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
		
		String s;
		BufferedLineReader in = new BufferedLineReader(fileName);
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			String peptide = token[2];
			int pepLength = 0;
			for(int i=0; i<peptide.length(); i++)
			{
				if(Character.isUpperCase(peptide.charAt(i)))
					pepLength++;
			}
			if(pepLength >= 10)
				out.println(s);
		}
		
		out.close();
		System.out.println("Done");
	}
	
	public static void extractShortPeptides() throws Exception
	{
		File dir = new File("/home/sangtaekim/Test/MSGFLib");
		String fileName = dir.getPath()+"/CID_IT_1Th.tsv";
		String outputFileName = dir.getPath() + "/CID_IT_1Th_Short.tsv";
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
		
		String s;
		BufferedLineReader in = new BufferedLineReader(fileName);
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				out.println(s);
			else
			{
				String[] token = s.split("\t");
				String peptide = token[7];
				int pepLength = 0;
				for(int i=0; i<peptide.length(); i++)
				{
					if(Character.isUpperCase(peptide.charAt(i)))
						pepLength++;
				}
				if(pepLength < 10)
					out.println(s);
			}
		}
		
		out.close();
		System.out.println("Done");
	}
	
	public static void makeTh() throws Exception
	{
		File dir = new File("/home/sangtaekim/Test/MSGFLib");
		dir = new File(System.getProperty("user.home")+"/Research/Data/MSGFLib");
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(dir.getPath()+"/CID_IT_1Th_2.tsv")));
		out.println("#SpecFile\tSpecIndex\tScan#\tFragMethod\tPrecursor\tPMError(Da)\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecProb\tP-value");
		for(File f : dir.listFiles())
		{
			if(!f.getName().endsWith("Da_2.tsv"))
				continue;
			
			String fileName = f.getName();
			System.out.println("Processing " + fileName);
			int tol = Integer.parseInt(fileName.substring("CID_IT_".length(), fileName.lastIndexOf("Da")));
			String s;
			BufferedLineReader in = new BufferedLineReader(f.getPath());
			while((s=in.readLine()) != null)
			{
				if(s.startsWith("#") || s.length() == 0)
					continue;
				String[] token = s.split("\t");
				if(token.length != 13)
					continue;
				int charge = Integer.parseInt(token[6]);
				if(charge == tol)
					out.println(s);
			}
			in.close();
		}
		
		out.close();
		System.out.println("Done");
	}
	
	public static void analyzeLibraryPSMs() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/SpecLib/MSGFOut_Human.tsv";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		Histogram<Integer> hist = new Histogram<Integer>();
		Histogram<Integer> histLength = new Histogram<Integer>();
		
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
			
//			if(specProbScore < 10)
				histLength.add(pepLength);
		}
//		hist.printSortedRatio();
		histLength.printSortedRatio();
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
