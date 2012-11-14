package edu.ucsd.msjava.ims;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.parser.BufferedLineReader;

public class DtaToMSGFInput {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit("Illegal parameter");
		
		File dtaFile = new File(argv[0]);
		if(!dtaFile.exists())
			printUsageAndExit("File does not exist.");
		
		File msgfInputFile = new File(argv[1]);
		if(msgfInputFile.exists())
			printUsageAndExit(msgfInputFile.getName()+ " already exists.");

		makeMSGFInput(dtaFile, msgfInputFile);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println(message);
		System.out.println("Usage: java DtaToMSGFInput DTAFile MSGFInputFile");
		System.exit(-1);
	}
	
	public static void makeMSGFInput(File dtaFile, File msgfInputFile) throws Exception
	{
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(msgfInputFile)));
		String header = "#SpectrumFile\tScan#\tAnnotation\tCharge\tFrameNum\tPrevSpecProb";
		out.println(header);
		
		String dtaFileName = dtaFile.getName();
		
		BufferedLineReader in = new BufferedLineReader(dtaFile.getPath());
		String s;
		
		int lineNum = 0;
		int specIndex = -1;
		while((s=in.readLine()) != null)
		{
			lineNum++;
			if(s.startsWith("==="))
			{
				++specIndex;
				String metaInfo = s.substring(s.indexOf('"')+1, s.lastIndexOf(".dta"));
				String[] token = metaInfo.split("\\.");
				if(token.length != 5)
				{
					System.out.println("Syntax Error in Line " + lineNum + ": " + s);
					System.exit(-1);
				}
				String pepSeq = token[0].replaceAll("!", "").replaceAll("@", "+15.995");
				int maxCharge = Integer.parseInt(token[1]);
				int frameNum = Integer.parseInt(token[2]);
				float prevSpecProb = Float.parseFloat(token[3]+token[4]);
				
				for(int charge=2; charge<=maxCharge; charge++)
					out.println(dtaFileName+"\t"+specIndex+"\t"+"."+pepSeq+".\t"+charge+"\t"+frameNum+"\t"+prevSpecProb);
			}
		}

		in.close();
		out.close();
	}
}
