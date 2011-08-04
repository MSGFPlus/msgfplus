package misc;

import java.io.File;

import fdr.TargetDecoyPSMSet;

import parser.BufferedLineReader;
import parser.InsPecTParser;

public class CountID {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length > 4 || argv.length == 0)
			printUsageAndExit();
			
		float threshold = 0.01f;
		if(argv.length >= 2)
			threshold = Float.parseFloat(argv[1]);
		String decoyPrefix = null;
		if(argv.length == 3)
			decoyPrefix = argv[2];
		countPeptide(argv[0], threshold, decoyPrefix);
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("usage: java CountPeptide MSGFDBORInsPecTResult [FDRThreshold] (0: PSMLevel, 1: PeptideLevel) [DecoyPrefix]");
		System.exit(-1);
	}
	
	public static void countPeptide(String fileName, float threshold, String decoyPrefix) throws Exception
	{
		File tempFile = File.createTempFile("MSGFDB", "tempResult");
		tempFile.deleteOnExit();
		
		int specFileColumn = -1;
		int specIndexColumn = -1;
		int scanNumColumn = -1;
		int annotationColumn = -1;
		int proteinColumn = -1;
		int specProbColumn = -1;
		int fScoreColumn = -1;
		int mqScoreColumn = -1;
		
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s = in.readLine();
		String[] label = s.split("\t");
		for(int i=0; i<label.length; i++)
		{
			if(label[i].equalsIgnoreCase(InsPecTParser.ANNOTATION) || label[i].equalsIgnoreCase("Peptide"))
				annotationColumn = i;
			else if(label[i].equalsIgnoreCase(InsPecTParser.PROTEIN))
				proteinColumn = i;
			else if(label[i].equalsIgnoreCase(InsPecTParser.SPEC_PROB))
				specProbColumn = i;
			else if(label[i].equalsIgnoreCase(InsPecTParser.F_SCORE))
				fScoreColumn = i;
			else if(label[i].equalsIgnoreCase(InsPecTParser.MQ_SCORE))
				mqScoreColumn = i;
			else if(label[i].equalsIgnoreCase(InsPecTParser.SPEC_FILE) || label[i].equalsIgnoreCase("#SpecFile"))
				specFileColumn = i;
			else if(label[i].equalsIgnoreCase(InsPecTParser.SPEC_INDEX))
				specIndexColumn = i;
			else if(label[i].equalsIgnoreCase(InsPecTParser.SCAN_NUM))
				scanNumColumn = i;
		}
		in.close();
		
		int scoreCol = -1;
		String scoreName = null;
		boolean isGreaterBetter = false;
		if(specProbColumn >= 0)
		{
			scoreCol = specProbColumn;
			isGreaterBetter = false;
			scoreName = "SpecProb";
		}
		else if(fScoreColumn >= 0)
		{
			scoreCol = fScoreColumn;
			isGreaterBetter = true;
			scoreName = "F-score";
		}
		else if(mqScoreColumn >= 0)
		{
			scoreCol = mqScoreColumn;
			isGreaterBetter = true;
			scoreName = "MQScore";
		}
		else
		{
			System.out.println("No score!");
			System.exit(-1);
		}
		
		int specFileCol = -1;
		if(specFileColumn >= 0)
		{
			specFileCol = specFileColumn;
		}
		else
		{
			System.out.println("No #SpectrumFile or #SpecFile!");
			System.exit(-1);
		}
		
		int specIndexCol = -1;
		if(specIndexColumn >= 0)
			specIndexCol = specIndexColumn;
		else if(scanNumColumn >= 0)
			specIndexCol = scanNumColumn;
		
		int pepCol = annotationColumn;
		int dbCol = proteinColumn;
		
		TargetDecoyPSMSet psmSet = new TargetDecoyPSMSet(
				new File(fileName), 
				"\t", 
				true,
				scoreCol, 
				isGreaterBetter, 
				specFileCol,
				specIndexCol, 
				pepCol,
				null,
				dbCol, 
				(decoyPrefix == null ? "REV" : decoyPrefix));

		System.out.println("Score: " + scoreName);
		System.out.println("Threshold: " + threshold);
		System.out.println("#PSMs: " + psmSet.getNumIdentifiedPSMs(threshold));
		System.out.println("#Peptides: " + psmSet.getNumIdentifiedPeptides(threshold));
	}
}
