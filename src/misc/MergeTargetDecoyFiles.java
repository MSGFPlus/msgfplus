package misc;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;

import parser.BufferedLineReader;

public class MergeTargetDecoyFiles {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit();
		mergeSearchResults(argv[0], argv[1]);
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("usage: java MergeTargetDecoyFiles targetFile decoyFile");
		System.exit(-1);
	}
	
	public static void mergeSearchResults(String targetResults, String decoyResults) throws Exception
	{
		HashMap<String,String> targetMap = new HashMap<String,String>();
		String s;
		BufferedLineReader in = new BufferedLineReader(targetResults);
		
		int specProbCol = -1;
		int specFileCol = -1;
		int scanNumCol = -1;
		String headerStr = in.readLine();	// header
		String[] header = headerStr.split("\t");
		for(int i=0; i<header.length; i++)
		{
			if(header[i].equalsIgnoreCase("SpecProb"))
				specProbCol = i;
			else if(header[i].contains("SpecFile") || header[i].contains("SpectrumFile"))
				specFileCol = i;
			else if(header[i].contains("Scan"))
				scanNumCol = i;
		}
		
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < specFileCol || token.length < specProbCol || token.length < scanNumCol)
				continue;
			String key = token[specFileCol]+":"+token[scanNumCol];
			if(targetMap.get(key) == null)
				targetMap.put(key, s);
		}
		in.close();
		
		in = new BufferedLineReader(decoyResults);
		System.out.println(in.readLine());
		int prevScanNum = -1;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < specFileCol || token.length < specProbCol || token.length < scanNumCol)
				continue;
			int scanNum = Integer.parseInt(token[scanNumCol]);
			if(scanNum == prevScanNum)
				continue;
			else
				prevScanNum = scanNum;
			String key = token[specFileCol]+":"+token[scanNumCol];
			
			String targetResult = targetMap.get(key);
			String decoyResult = s;
			if(targetResult == null)
			{
				System.out.println(decoyResult);
				continue;
			}
			else
			{
				float targetSpecProb = Float.parseFloat(targetResult.split("\t")[specProbCol]);
				float decoySpecProb = Float.parseFloat(token[specProbCol]);
				if(targetSpecProb <= decoySpecProb)
					System.out.println(targetResult);
				else if(targetSpecProb > decoySpecProb)
					System.out.println(decoyResult);
			}
		}		
		
		in.close();
	}
}
