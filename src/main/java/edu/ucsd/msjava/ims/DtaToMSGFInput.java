package edu.ucsd.msjava.ims;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Collections;

import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.PNNLSpectrumParser;

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
		String header = "#SpectrumFile\tScan#\tAnnotation\tCharge\tFrameNum\tDtaIndex\tPrevSpecProb";
		msgfOut.println(header);
		
		String dtaFilePath = dtaFile.getAbsolutePath();
		File mgfFile = new File(dtaFilePath.substring(0, dtaFilePath.lastIndexOf("_dta.txt")) + ".mgf");
		String mgfFileName = mgfFile.getName();
		
		BufferedLineReader in = new BufferedLineReader(dtaFile.getPath());
		String s;
		
		int lineNum = 0;
		int specIndex = 0;
		int origSpecIndex = 0;
		while((s=in.readLine()) != null)
		{
			lineNum++;
			if(s.startsWith("==="))
			{
				++origSpecIndex;
				String metaInfo = s.substring(s.indexOf('"')+1, s.lastIndexOf(".dta"));
				String[] token = metaInfo.split("\\.");
				if(token.length != 5 && token.length != 4)
				{
					System.out.println("Syntax Error in Line " + lineNum + ": " + s);
					System.exit(-1);
				}
				String pepSeq = token[0].replaceAll("!", "").replaceAll("@", "+15.995");
				int maxCharge = Integer.parseInt(token[1]);
				int frameNum = Integer.parseInt(token[2]);
				float prevSpecProb;
				if(token.length == 5)
					prevSpecProb = Float.parseFloat(token[3]+token[4]);
				else
					prevSpecProb = Float.parseFloat(token[3]);
				
				for(int charge=2; charge<=maxCharge; charge++)
				{		
					++specIndex;
					msgfOut.println(mgfFileName+"\t"+specIndex+"\t"+"."+pepSeq+".\t"+charge+"\t"+frameNum+"\t"+origSpecIndex+"\t"+prevSpecProb);
				}
			}
		}

		in.close();
		
		msgfOut.close();
		
		PrintStream mgfOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(mgfFile)));
		
		SpectraIterator itr = new SpectraIterator(dtaFile.getPath(), new PNNLSpectrumParser());
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int maxCharge = spec.getCharge();
			for(int charge=2; charge<=maxCharge; charge++)
			{		
				float precursorMz = spec.getParentMass()/charge+(float)Composition.PROTON;
				spec.setPrecursor(new Peak(precursorMz, 0, charge));
				spec.outputMgf(mgfOut);
			}			
		}
		mgfOut.close();
	}
}
