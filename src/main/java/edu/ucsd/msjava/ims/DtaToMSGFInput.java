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
		PrintStream msgfOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(msgfInputFile)));
		String header = "#SpectrumFile\tScan#\tAnnotation\tCharge\tFrameNum\tFromScan\tToScan\tPrevSpecProb";
		msgfOut.println(header);
		
		String dtaFilePath = dtaFile.getAbsolutePath();
//		File mgfFile = new File(dtaFilePath.substring(0, dtaFilePath.lastIndexOf("_dta.txt")) + ".mgf");
//		String mgfFileName = mgfFile.getName();
		
		BufferedLineReader in = new BufferedLineReader(dtaFile.getPath());
		String s;
		
		int lineNum = 0;
		int origSpecIndex = 0;
		while((s=in.readLine()) != null)
		{
			lineNum++;
			if(s.startsWith("==="))
			{
				++origSpecIndex;
				String metaInfo = s.substring(s.indexOf('"')+1, s.lastIndexOf(".dta"));
				String[] token = metaInfo.split("\\.");
				if(token.length != 8 && token.length != 9)
				{
					System.out.println("Syntax Error in Line " + lineNum + ": " + s);
					System.exit(-1);
				}
				String annotation = token[0]+"."+token[1].replaceAll("!", "").replaceAll("@", "+15.995").replaceAll("\\*", "+15.995")+"."+token[2];
				int charge = Integer.parseInt(token[3]);
				int frameNum = Integer.parseInt(token[4]);
				int fromScan = Integer.parseInt(token[5]);
				int toScan = Integer.parseInt(token[6]);
				float prevSpecProb;
				if(token.length == 9)
					prevSpecProb = Float.parseFloat(token[7]+token[8]);
				else
					prevSpecProb = Float.parseFloat(token[7]);
				
				msgfOut.println(dtaFile.getName()+"\t"+origSpecIndex+"\t"+annotation+"\t"+charge+"\t"+frameNum+"\t"+fromScan+"\t"+toScan+"\t"+prevSpecProb);
			}
		}

		in.close();
		msgfOut.close();
	}
}
