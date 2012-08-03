package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashSet;

import edu.ucsd.msjava.parser.MascotParser;
import edu.ucsd.msjava.parser.PSM;
import edu.ucsd.msjava.parser.PSMList;


public class DatToTxt {
	public static void main(String argv[]) throws Exception
	{
//		String fileName = System.getProperty("user.home")+"/Research/Data/HeckWhole/ResultsShabaz/LysN_ETD_F245063.dat";
		if(argv.length != 2 && argv.length != 3)
		{
			printUsageAndExit("Illegal parameters");
		}
		File datFile = new File(argv[0]);
		if(!datFile.isFile())
			printUsageAndExit(argv[0] + " is not a file.");
		File txtFile = new File(argv[1]);
		boolean isDecoy = false;
		if(argv.length == 3 && argv[2].equalsIgnoreCase("1"))
			isDecoy = true;
			
		datToTxt(datFile, txtFile, isDecoy);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.err.println(message);
		System.out.println("usage: java -Xmx3000M DatToTxt datFileName outputFileName [0/1] (0: target, 1:decoy)");
		System.exit(-1);
	}
	public static void datToTxt(File datFile, File txtFile, boolean isDecoy) throws Exception
	{
		PrintStream out;
		if(txtFile == null)
			out = System.out;
		else
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(txtFile)));
		out.println("#Title\tCharge\tAnnotation\tMascotScore\tProtein");
		HashSet<String> titleSet = new HashSet<String>();
		PSMList<PSM> psmList = MascotParser.parseFromDat(datFile.getPath(), isDecoy);
		for(PSM psm : psmList)
		{
			/*
			System.out.print(psm.getTitle());
			System.out.print("\t"+psm.getCharge());
			System.out.print(psm.getPrecedingResidue());
			System.out.print("."+psm.getPeptide()+".");
			System.out.print(psm.getSucceedingResidue());
			System.out.print("\t"+psm.getRawScore()+"\t"+psm.getProtein());
			System.out.println();
			*/
			if(titleSet.contains(psm.getTitle()))
				continue;
			else
				titleSet.add(psm.getTitle());
			out.print(psm.getTitle()+"\t"+psm.getCharge()+"\t"+psm.getPrecedingResidue()+"."+psm.getPeptideStr()+"."+psm.getSucceedingResidue());
			out.print("\t"+psm.getRawScore()+"\t"+psm.getProtein());
			out.println();
		}
		out.close();
	}
}
