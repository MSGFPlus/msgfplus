package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import parser.BufferedLineReader;

public class MSGFDBToInspect {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit(null);
		else
		{
			File msgfOutput = new File(argv[0]);
			if(msgfOutput == null || !msgfOutput.exists())
				printUsageAndExit(argv[0] + " not found!");
			File inspectOutput = new File(argv[1]);
			convert(msgfOutput, inspectOutput);
		}
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println(message);
		System.out.println("usage: java MSGFDBToInspect MSGFDBOutput InsPecTOutput");
		System.exit(-1);
	}
	
	public static void convert(File msgfOutput, File inspectOutput) throws Exception
	{
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(inspectOutput)));
		
		String s;
		BufferedLineReader in = new BufferedLineReader(msgfOutput.getPath());
		String header = in.readLine();
		String inspectHeader = "#SpectrumFile\tScan#\tAnnotation\tProtein\tCharge\tMQScore\tLength\tTotalPRMScore\tMedianPRMScore\tFractionY\tFractionB\tIntensity\tNTT\tp-value\tF-Score\tDeltaScore\tDeltaScoreOther\tRecordNumber\tDBFilePos\tSpecFilePos\tPrecursorMZ\tPrecursorMZError";
		out.println(inspectHeader);
//		String[] headerToken = in.readLine().split("\t");
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length != 12)
				continue;
			String specFile = token[0];
			String scanNum = token[1];
			String title = token[2];
			String actMethod = token[3];
			String precursor = token[4];
			String pmError = token[5];
			String charge = token[6];
			String peptide = token[7];
			String pepStr = peptide.substring(peptide.indexOf('.')+1, peptide.lastIndexOf('.'));
			String protein = token[8];
			String deNovoScore = token[9];
			String msgfScore = token[10];
			String specProb = token[11];
			out.print(specFile+"\t"+scanNum+"\t"+peptide+"\t"+protein+"\t"+charge+"\t"+0+"\t"+pepStr.length()+"\t"+0+"\t"+0+"\t"+0+"\t"+0+"\t"+0+"\t"+0);
			out.print("\t"+specProb+"\t"+0+"\t"+0+"\t"+0+"\t"+0+"\t"+0+"\t"+0+"\t"+precursor+"\t"+pmError);
			out.println();
		}
		out.close();
	}
}
