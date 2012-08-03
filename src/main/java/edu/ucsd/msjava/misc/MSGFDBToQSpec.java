package edu.ucsd.msjava.misc;

import java.util.HashMap;
import java.util.TreeSet;

import edu.ucsd.msjava.msgf.Histogram;
import edu.ucsd.msjava.parser.BufferedLineReader;

public class MSGFDBToQSpec {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2 && argv.length != 3)
			printUsageAndExit("Wrong parameters!");
		convert(argv[0], argv[1], null);
	}
	
	public static void printUsageAndExit(String message) throws Exception
	{
		System.out.println(message);
		System.out.println("Usage: java -jar MSGFDB.jar misc.MSGFDBToQSpec MSGFDBResultFileName DatabaseFileName [OutputFileName]");
		System.exit(-1);
	}
	
	public static void convert(String msgfdbFileName, String dbFileName, String outputFileName) throws Exception
	{
		BufferedLineReader in = new BufferedLineReader(msgfdbFileName);
		String s;
		String header = in.readLine();
		
		System.out.println("protid\tprotLen\t0\t1");
		Histogram<String> control = new Histogram<String>();
		Histogram<String> treatment = new Histogram<String>();
		
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			float fdr = Float.parseFloat(token[12]);
			if(fdr > 0.01)
				continue;
			String protein = token[7];
			protein = protein.split("\\s+")[0].trim();
			String annotation = token[6];
			String pepStr = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
			if(pepStr.endsWith("K+8.014") || pepStr.endsWith("R+10.008"))
				treatment.add(protein);
			else
				control.add(protein);
		}
		
		HashMap<String,Integer> protLengthMap = getAnnotationProtLengthMap(dbFileName);
		
		TreeSet<String> protSet = new TreeSet<String>();
		for(String prot : treatment.keySet())
			protSet.add(prot);
		for(String prot : control.keySet())
			protSet.add(prot);
		for(String prot : protSet)
			System.out.println(prot+"\t"+protLengthMap.get(prot)+"\t"+control.get(prot)+"\t"+treatment.get(prot));
	}
	
	public static HashMap<String,Integer> getAnnotationProtLengthMap(String dbFileName) throws Exception
	{
		HashMap<String,Integer> map = new HashMap<String,Integer>();
		
		BufferedLineReader in = new BufferedLineReader(dbFileName);
		String s;
		String annotation = null;
		int length = 0;
		while((s=in.readLine()) != null)
		{
			if(s.length() == 0)
				continue;
			if(s.startsWith(">"))
			{
				if(annotation != null)
					map.put(annotation, length);
				length = 0;
				annotation = s.substring(1).split("\\s+")[0].trim();
			}
			else
			{
				length += s.trim().length();
			}
		}
		if(annotation != null)
			map.put(annotation, length);
			
		return map;
	}
}
