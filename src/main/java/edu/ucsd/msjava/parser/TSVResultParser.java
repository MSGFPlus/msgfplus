package edu.ucsd.msjava.parser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class TSVResultParser {
	private File tsvFile;
	private Set<String> pepSet;
	private Set<String> scanSet;
	private Set<String> idSet;
	private Map<String, Float> idToSpecEValue;
	
	public TSVResultParser(File tsvFile)
	{
		this.tsvFile = tsvFile;
	}
	
	public Set<String> getPepSet()
	{
		return pepSet;
	}
	
	public Set<String> getScanSet()
	{
		return scanSet;
	}
	
	public Set<String> getIdSet()
	{
		return idSet;
	}
	
	public Float getSpecEValue(String id)
	{
		return idToSpecEValue.get(id);
	}
	
	public String parse(float fdrThreshold)
	{
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(tsvFile.getPath());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String header = in.readLine();
		if(!header.startsWith("#") && !header.startsWith("Result"))
			return "No header!";
		
		String[] headerToken = header.split("\t");
		int specQValueColNum = -1;
		int pepQValueColNum = -1;
		int pepColNum = -1;
		int scanNumCol = -1;
		int idCol = -1;
		int specEValueCol = -1;
		for(int i=0; i<headerToken.length; i++)
		{
			if(headerToken[i].equalsIgnoreCase("FDR") || headerToken[i].equalsIgnoreCase("QValue") || headerToken[i].equalsIgnoreCase("q-value"))
				specQValueColNum = i;
			if(headerToken[i].equalsIgnoreCase("PepFDR") || headerToken[i].equalsIgnoreCase("PepQValue"))
				pepQValueColNum = i;
			if(headerToken[i].equalsIgnoreCase("Peptide") || headerToken[i].equalsIgnoreCase("Annotation"))
				pepColNum = i;
			if(headerToken[i].equalsIgnoreCase("ScanNum") || headerToken[i].equalsIgnoreCase("Scan#") || headerToken[i].equalsIgnoreCase("Scan"))
				scanNumCol = i;
			if(headerToken[i].equalsIgnoreCase("SpecID"))
				idCol = i;
			if(headerToken[i].equalsIgnoreCase("SpecEValue") || headerToken[i].equalsIgnoreCase("SpecProb"))
				specEValueCol = i;
		}
		if(specQValueColNum < 0)
			return "QValue column is missing!";
		if(pepQValueColNum < 0)
			return "PepQValue column is missing!";
		if(pepColNum < 0)
			return "Annotation column is missing!";
		if(scanNumCol < 0)
			return "Scan column is missing!";
		if(idCol < 0)
			return "SpecID column is missing!";
		if(specEValueCol < 0)
			return "SpecEValue column is missing!";
		
		String s;
		pepSet = new HashSet<String>();
		scanSet = new HashSet<String>();
		idSet = new HashSet<String>();
		idToSpecEValue = new HashMap<String,Float>();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length <= specQValueColNum || token.length <= pepQValueColNum || token.length <= pepColNum 
					|| token.length <= idCol || token.length <= specEValueCol)
				continue;
			double specQValue = Double.parseDouble(token[specQValueColNum]);
			double pepQValue = Double.parseDouble(token[pepQValueColNum]);
			float specEValue = Float.parseFloat(token[specEValueCol]);
//			if(token[scanNumCol].equals("6804"))
//				System.out.println("Debug");
			idToSpecEValue.put(token[idCol], specEValue);

			if(specQValue <= fdrThreshold)
			{
				scanSet.add(token[scanNumCol]);
				idSet.add(token[idCol]);
			}
			if(pepQValue <= fdrThreshold)
			{
				String annotation = token[pepColNum];
				
				String pepStr;
				
				if(annotation.matches("[A-Z\\-_]?\\..+\\.[A-Z\\-_]?"))
					pepStr = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
				else
					pepStr = annotation;
				
				StringBuffer unmodStr = new StringBuffer();
				for(int i=0; i<pepStr.length(); i++)
					if(Character.isLetter(pepStr.charAt(i)))
						unmodStr.append(pepStr.charAt(i));
				
				pepSet.add(unmodStr.toString());
			}			
		}
		
		try {
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return null;
	}
}
