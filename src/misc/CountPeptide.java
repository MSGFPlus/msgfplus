package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import msgf.Histogram;
import parser.BufferedLineReader;

public class CountPeptide {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit();
		countPeptide(argv[0], Double.parseDouble(argv[1]));
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("usage: java CountPeptide MSGFDBResult.txt FDRThreshold");
		System.exit(-1);
	}
	
	public static void countPeptide(String fileName, double threshold) throws Exception
	{
		BufferedLineReader in = new BufferedLineReader(fileName);
		String header = in.readLine();
		if(header == null || !header.startsWith("#"))
		{
			System.out.println("Not a valid MSGFDB result file!");
			System.exit(0);
		}
//		String[] headerToken = header.split("\t");
//		int eFDRColNum = -1;
//		int pepColNum = -1;
//		for(int i=0; i<headerToken.length; i++)
//		{
//			if(headerToken[i].equalsIgnoreCase("EFDR") || headerToken[i].equalsIgnoreCase("FDR"))
//				eFDRColNum = i;
//			if(headerToken[i].equalsIgnoreCase("Peptide") || headerToken[i].equalsIgnoreCase("Annotation"))
//				pepColNum = i;
//		}
//		if(eFDRColNum < 0)
//		{
//			System.out.println("FDR column is missing!");
//			System.exit(0);
//		}
//		if(pepColNum < 0)
//		{
//			System.out.println("Annotation column is missing!");
//			System.exit(0);
//		}

		File tempFile = File.createTempFile("MSGFDB", "tempResult");
		tempFile.deleteOnExit();
		
		int specIndexCol = 1;
		int pepCol = 7;
		int dbCol = 8;
		int scoreCol = 11;
		fdr.ComputeFDR.computeFDR(tempFile, null, scoreCol, false, "\t", 
				specIndexCol, pepCol, null, true, true, 
				true, dbCol, "REV_",
				1, 0.011, fileName);

	} catch (IOException e) {
		e.printStackTrace();
		
		int totalID = 0;
		int numID = 0;
		String s;
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
			}
			
//			if(ntt == 0) {
//				System.out.println(s);
//				System.exit(0);
//			}
		}
		
		System.out.println("TotalPSM\t" + totalID);
		System.out.println("NumID\t" + numID+"\t"+numID/(float)totalID);
		System.out.println("Cleavage hist");
		nttHist.printSortedRatio();
	}
}
