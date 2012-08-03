package edu.ucsd.msjava.misc;

import java.util.HashSet;

import edu.ucsd.msjava.msgf.Histogram;
import edu.ucsd.msjava.parser.BufferedLineReader;

public class CountPSMs {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit();
		countID(argv[0], Double.parseDouble(argv[1]));
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
		if(header == null || (!header.startsWith("#") && !header.startsWith("PSMId")))
		{
			System.out.println("Not a valid MSGFDB result file!");
			System.exit(0);
		}
		String[] headerToken = header.split("\t");
		int eFDRColNum = -1;
		int pepColNum = -1;
		for(int i=0; i<headerToken.length; i++)
		{
			if(headerToken[i].equalsIgnoreCase("EFDR") || headerToken[i].equalsIgnoreCase("FDR") || headerToken[i].equalsIgnoreCase("q-value"))
				eFDRColNum = i;
			if(headerToken[i].equalsIgnoreCase("Peptide") || headerToken[i].equalsIgnoreCase("Annotation"))
				pepColNum = i;
		}
		if(eFDRColNum < 0)
		{
			System.out.println("FDR column is missing!");
			System.exit(0);
		}
		if(pepColNum < 0)
		{
			System.out.println("Annotation column is missing!");
			System.exit(0);
		}
		
		int totalID = 0;
		int numID = 0;
		String s;
		HashSet<String> pepSet = new HashSet<String>();
		Histogram<Integer> nttHist = new Histogram<Integer>();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length <= eFDRColNum || token.length <= pepColNum)
				continue;
			double eFDR = Double.parseDouble(token[eFDRColNum]);
			totalID++;
			if(eFDR <= threshold)
			{
				numID++;
				int ntt=0;
				String annotation = token[pepColNum];
				char pre = annotation.charAt(0);
				if(pre == 'K' || pre == 'R' || pre == '_')
					ntt+=2;
				String pepStr = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
				StringBuffer unmodStr = new StringBuffer();
				for(int i=0; i<pepStr.length(); i++)
					if(Character.isLetter(pepStr.charAt(i)))
						unmodStr.append(pepStr.charAt(i));
				char last = unmodStr.charAt(unmodStr.length()-1);
				if(last == 'K' || last == 'R')
					ntt+=1;
				nttHist.add(ntt);
				pepSet.add(unmodStr.toString());
			}
			
//			if(ntt == 0) {
//				System.out.println(s);
//				System.exit(0);
//			}
		}
		
		System.out.println("TotalPSM\t" + totalID);
		System.out.println("NumID\t" + numID+"\t"+numID/(float)totalID);
		System.out.println("NumPeptides\t" + pepSet.size());
		System.out.println("Cleavage hist");
		nttHist.printSortedRatio();
	}
}
