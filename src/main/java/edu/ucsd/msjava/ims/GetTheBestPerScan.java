package edu.ucsd.msjava.ims;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.ucsd.msjava.parser.BufferedLineReader;

public class GetTheBestPerScan {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1)
			printUsageAndExit("Illegal parameter.");
		File tsvFile = new File(argv[0]);
		if(!tsvFile.exists())
			printUsageAndExit("File does not exist.");
		getTheBest(tsvFile);
	}
		
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println(message);
		System.out.println("Usage: java GetTheBestPerScan TSVFile");
		System.exit(-1);
	}
	
	public static void getTheBest(File tsvFile) throws Exception
	{
		BufferedLineReader in = new BufferedLineReader(tsvFile.getPath());
		String header = in.readLine();
		String[] headerToken = header.split("\t");
		
		int dtaIndexCol = -1;
		int specProbCol = -1;
		for(int i=0; i<headerToken.length; i++)
		{
			if(headerToken[i].equalsIgnoreCase("DtaIndex"))
				dtaIndexCol = i;
			if(headerToken[i].equalsIgnoreCase("SpecProb"))
				specProbCol = i;
		}
		
		if(dtaIndexCol == -1)
		{
			System.err.println("DtaIndex column does not exist.");
			System.exit(-1);
		}
		if(specProbCol == -1)
		{
			System.err.println("SpecProb column does not exist.");
			System.exit(-1);
		}
		
		Map<Integer,String> table = new HashMap<Integer,String>();

		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			int dtaIndex = Integer.parseInt(token[dtaIndexCol]);
			String prev = table.get(dtaIndex);
			if(prev == null)
				table.put(dtaIndex, s);
			else
			{
				String[] tokenPrev = prev.split("\t");
				float prevSpecProb = Float.parseFloat(tokenPrev[specProbCol]);
				float specProb = Float.parseFloat(token[specProbCol]);
				if(specProb < prevSpecProb)
					table.put(dtaIndex, s);
			}
		}
		
		in.close();
		
		List<Integer> indexList = new ArrayList<Integer>(table.keySet());
		Collections.sort(indexList);
		System.out.println(header);
		for(int dtaIndex : indexList)
			System.out.println(table.get(dtaIndex));
	}
}
