package edu.ucsd.msjava.misc;

import java.util.HashSet;

import edu.ucsd.msjava.msgf.Histogram;
import edu.ucsd.msjava.parser.BufferedLineReader;

public class CountPSMs {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1 && argv.length != 2)
			printUsageAndExit();
		double threshold = 0.01;
		if(argv.length == 2)
			threshold = Double.parseDouble(argv[1]);
		countID(argv[0], threshold);
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("usage: java CountID MSGFDBResult.txt FDRThreshold");
		System.exit(-1);
	}
	
	public static void countID(String fileName, double threshold) throws Exception
	{
		BufferedLineReader in = new BufferedLineReader(fileName);
		String header = in.readLine();
		if(header == null) // || (!header.startsWith("#") && !header.startsWith("PSMId")))
		{
			System.out.println("Not a valid MSGFDB result file!");
			System.exit(0);
		}
		String[] headerToken = header.split("\t");
		int specQValueColNum = -1;
		int pepQValueColNum = -1;
		int pepColNum = -1;
		int specIDCol = -1;
		for(int i=0; i<headerToken.length; i++)
		{
			if(headerToken[i].equalsIgnoreCase("FDR") || headerToken[i].equalsIgnoreCase("QValue") || headerToken[i].equalsIgnoreCase("q-value"))
				specQValueColNum = i;
			if(headerToken[i].equalsIgnoreCase("PepFDR") || headerToken[i].equalsIgnoreCase("PepQValue"))
				pepQValueColNum = i;
			if(headerToken[i].equalsIgnoreCase("Peptide") || headerToken[i].equalsIgnoreCase("Annotation"))
				pepColNum = i;
			if(headerToken[i].equalsIgnoreCase("SpecID") || headerToken[i].equalsIgnoreCase("ScanNum") || headerToken[i].equalsIgnoreCase("Scan#") || headerToken[i].equalsIgnoreCase("Scan"))
				specIDCol = i;
		}
		if(specQValueColNum < 0)
		{
			System.out.println("QValue column is missing!");
			System.exit(0);
		}
		if(pepQValueColNum < 0)
		{
			System.out.println("PepQValue column is missing!");
			System.exit(0);
		}
		if(pepColNum < 0)
		{
			System.out.println("Annotation column is missing!");
			System.exit(0);
		}
		if(specIDCol < 0)
		{
			System.out.println("SpecID or Scan# column is missing!");
			System.exit(0);
		}
		
		String s;
		HashSet<String> pepSet = new HashSet<String>();
		HashSet<String> idScanSet = new HashSet<String>();
		HashSet<String> scanSet = new HashSet<String>();
		Histogram<Integer> nttHist = new Histogram<Integer>();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length <= specQValueColNum || token.length <= pepQValueColNum || token.length <= pepColNum)
				continue;
			double specQValue = Double.parseDouble(token[specQValueColNum]);
			double pepQValue = Double.parseDouble(token[pepQValueColNum]);
			scanSet.add(token[specIDCol]);
			if(specQValue <= threshold)
			{
				idScanSet.add(token[specIDCol]);
			}
			if(pepQValue <= threshold)
			{
				int ntt=0;
				String annotation = token[pepColNum];
				
				String pepStr;
				
				if(annotation.matches("[A-Z\\-_]?\\..+\\.[A-Z\\-_]?"))
				{
					char pre = annotation.charAt(0);
					if(pre == 'K' || pre == 'R' || pre == '_')
						ntt+=2;
					pepStr = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
					char last = annotation.charAt(annotation.lastIndexOf('.')-1);
					if(last == 'K' || last == 'R')
						ntt+=1;
					nttHist.add(ntt);
				}
				else
				{
					pepStr = annotation;
				}
				
				StringBuffer unmodStr = new StringBuffer();
				for(int i=0; i<pepStr.length(); i++)
					if(Character.isLetter(pepStr.charAt(i)))
						unmodStr.append(pepStr.charAt(i));
				
				pepSet.add(unmodStr.toString());
			}			
		}
		
		System.out.println("TotalPSM\t" + scanSet.size());
		System.out.println("NumID\t" + idScanSet.size()+"\t"+idScanSet.size()/(float)scanSet.size());
		System.out.println("NumPeptides\t" + pepSet.size());
		System.out.println("Cleavage hist");
		nttHist.printSortedRatio();
		
		in.close();
	}
}
