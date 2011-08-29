package misc;

import java.util.HashMap;

import parser.BufferedLineReader;

public class MergeTargetDecoyFiles {
	public static void main(String argv[]) throws Exception
	{
		boolean isMascot = false;
		if(argv.length != 2 && argv.length != 3)
			printUsageAndExit();
		if(argv.length == 3 && argv[2].equals("1"))
			isMascot = true;
			
		if(!isMascot)
			mergeSearchResults(argv[0], argv[1]);
		else
			mergeMascotSearchResults(argv[0], argv[1]);
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("usage: java MergeTargetDecoyFiles targetFile decoyFile [0/1] (0: MSGFDB, 1: Mascot)");
		System.exit(-1);
	}
	
	public static void mergeSearchResults(String targetResults, String decoyResults) throws Exception
	{
		HashMap<String,String> targetMap = new HashMap<String,String>();
		String s;
		BufferedLineReader in = new BufferedLineReader(targetResults);
		
		int specProbCol = -1;
		int specFileCol = -1;
		int specIndexCol = -1;
		String headerStr = in.readLine();	// header
		String[] header = headerStr.split("\t");
		for(int i=0; i<header.length; i++)
		{
			if(header[i].equalsIgnoreCase("SpecProb"))
				specProbCol = i;
			else if(header[i].contains("SpecFile") || header[i].contains("SpectrumFile"))
				specFileCol = i;
			else if(header[i].contains("SpecIndex"))
				specIndexCol = i;
		}
		
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < specFileCol || token.length < specProbCol || token.length < specIndexCol)
				continue;
			String key = token[specFileCol]+":"+token[specIndexCol];
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
			if(token.length < specFileCol || token.length < specProbCol || token.length < specIndexCol)
				continue;
			int scanNum = Integer.parseInt(token[specIndexCol]);
			if(scanNum == prevScanNum)
				continue;
			else
				prevScanNum = scanNum;
			String key = token[specFileCol]+":"+token[specIndexCol];
			
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
	
	public static void mergeMascotSearchResults(String targetResults, String decoyResults) throws Exception
	{
		HashMap<String,String> resultMap = new HashMap<String,String>();
		String s;
		BufferedLineReader in = new BufferedLineReader(targetResults);
		
		int scoreCol = -1;
		int titleCol = -1;
		String headerStr = in.readLine();	// header
		String[] header = headerStr.split("\t");
		for(int i=0; i<header.length; i++)
		{
			if(header[i].equalsIgnoreCase("MascotScore"))
				scoreCol = i;
			else if(header[i].contains("Title"))
				titleCol = i;
		}
		
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < titleCol || token.length < scoreCol)
				continue;
			String title = token[titleCol];
			if(resultMap.get(title) == null)
				resultMap.put(title, s);
		}
		in.close();
		
		in = new BufferedLineReader(decoyResults);
		// header
		System.out.println(in.readLine());
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < titleCol || token.length < scoreCol)
				continue;
			String title = token[titleCol];
			
			String prevResult = resultMap.get(title);
			if(prevResult == null)
			{
				resultMap.put(title, s);
			}
			else
			{
				float prevScore = Float.parseFloat(prevResult.split("\t")[scoreCol]);
				float curScore = Float.parseFloat(token[scoreCol]);
				if(curScore > prevScore)
					resultMap.put(title, s);
			}
		}	
		
		// results
		for(String result : resultMap.values())
			System.out.println(result);
		
		in.close();
	}	
}
