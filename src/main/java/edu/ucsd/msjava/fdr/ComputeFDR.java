package edu.ucsd.msjava.fdr;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

public class ComputeFDR {
	public static final float FDR_REPORT_THRESHOLD = 0.1f;
	public static void main(String argv[]) throws Exception
	{
		// required
		File targetFile = null;
		int scoreCol = -1;
		int specFileCol = -1;
		
		// optional
		File outputFile = null;
		boolean isGreaterBetter = false;
		boolean hasHeader = true;
		File decoyFile = null;
		String delimeter = "\t";
		int pepCol = -1;
		int specIndexCol = -1;
		boolean isConcatenated = false;
		boolean includeDecoy = false;
		
		int dbCol = -1;
		String decoyPrefix = null;
		float fdrThreshold = 1;
		float pepFDRThreshold = 1;
		
		ArrayList<Pair<Integer,ArrayList<String>>> reqStrList = new ArrayList<Pair<Integer,ArrayList<String>>>();
		
		int i=0;
		while(i<argv.length)
		{
     		// 	-f resuleFileName dbCol decoyPrefix or -f targetFileName decoyFileName 
     		if(argv[i].equalsIgnoreCase("-f"))
			{
				if(i+2 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				targetFile = new File(argv[i+1]);
				if(!targetFile.exists())
					printUsageAndExit(argv[i+1] + " doesn't exist.");
				else if(!targetFile.isFile())
					printUsageAndExit(argv[i+1] + " is not a file.");
				if(i+3 < argv.length && !argv[i+3].startsWith("-"))	// concatenated; -f resultFileName dbCol decoyPrefix
				{
					dbCol = Integer.parseInt(argv[i+2]);
					decoyPrefix = argv[i+3];
					isConcatenated = true;
					i+=4;
				}
				else	// separate; -f targetFileName decoyFileName
				{
					decoyFile = new File(argv[i+2]);
					if(!decoyFile.exists())
						printUsageAndExit(argv[i+2] + " doesn't exist.");
					else if(!decoyFile.isFile())
						printUsageAndExit(argv[i+2] + " is not a file.");
					isConcatenated = false;
					i+=3;
				}
			}
			else if(argv[i].equalsIgnoreCase("-s"))
			{
				if(i+2 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				try {
					scoreCol = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal scoreCol: " + argv[i+1]);
				}
				if(argv[i+2].equalsIgnoreCase("1"))
					isGreaterBetter = true;
				else
					isGreaterBetter = false;
				i+=3;
			}
			else if(argv[i].equalsIgnoreCase("-o"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				outputFile = new File(argv[i+1]);
				i += 2;
			}
			else if(argv[i].equalsIgnoreCase("-h"))
			{
				if(argv[i+1].equalsIgnoreCase("0"))
					hasHeader = false;
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-decoy"))
			{
				if(argv[i+1].equalsIgnoreCase("1"))
					includeDecoy = true;
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-delim"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				delimeter = argv[i+1];
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-p"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				try {
					pepCol = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal pepCol: " + argv[i+1]);
				}
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-n"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				try {
					specIndexCol = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal pepCol: " + argv[i+1]);
				}
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-i"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				try {
					specFileCol = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal pepCol: " + argv[i+1]);
				}
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-m"))
			{
				int matchCol = -1;
				if(i+2 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				try {
					matchCol = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal matchCol: " + argv[i+1]);
				}
				String[] token = argv[i+2].split(",");
				ArrayList<String> reqStrOrList = new ArrayList<String>();
				for(String s : token)
					reqStrOrList.add(s);
				reqStrList.add(new Pair<Integer,ArrayList<String>>(matchCol, reqStrOrList));
				i+=3;
			}
			else if(argv[i].equalsIgnoreCase("-fdr"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				try {
					fdrThreshold = Float.parseFloat(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal pepCol: " + argv[i+1]);
				}
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-pepfdr"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				try {
					pepFDRThreshold = Float.parseFloat(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal pepCol: " + argv[i+1]);
				}
				i+=2;
			}
			else
			{
				printUsageAndExit("Illegal parameter");
			}
		}
		
		if(targetFile == null)
			printUsageAndExit("Target is missing!");
//		if(specFileCol < 0)
//			printUsageAndExit("specFileCol is missing or illegal!");
		if(scoreCol < 0)
			printUsageAndExit("scoreCol is missing or illegal!");
		if(pepCol < 0)
			printUsageAndExit("pepCol is missing or illegal!");
		if(specIndexCol < 0)
			printUsageAndExit("specIndexCol is missing or illegal!");
		
		computeFDR(targetFile, decoyFile, 
				scoreCol, isGreaterBetter, 
				delimeter, specFileCol, specIndexCol, pepCol, reqStrList, 
				isConcatenated, includeDecoy, hasHeader, dbCol, decoyPrefix, fdrThreshold, pepFDRThreshold, outputFile);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.err.println(message);
		System.out.print("Usage: java -cp MSGFDB.jar fdr.ComputeFDR\n" +
				"\t -f resuleFileName protCol decoyPrefix or -f targetFileName decoyFileName\n" +
				"\t -i specFileCol (SpecFile column number)\n" +
				"\t -n specIndexCol (specIndex column number)\n" +
				"\t -p pepCol (peptide column number)\n" +
				"\t -s scoreCol 0/1 (0: smaller better, 1: greater better)\n" +
				"\t [-o outputFileName (default: stdout)]\n" +
				"\t [-delim delimeter] (default: \\t)\n" +
				"\t [-m colNum keyword (the column 'colNum' must contain 'keyword'. If 'keyword' is delimeted by ',' (e.g. A,B,C), then at least one must be matched.)]\n" +
				"\t [-h 0/1] (0: no header, 1: header (default))\n" +
				"\t [-fdr fdrThreshold]\n" +
				"\t [-pepfdr pepFDRThreshod]\n" +
				"\t [-decoy 0/1 (0: include decoy, 1: don't include decoy (default))\n"
				);
		System.exit(-1);
	}

	public static void computeFDR(File targetFile, File decoyFile, int scoreCol, boolean isGreaterBetter, String delimeter, 
			int specFileCol, int specIndexCol, int pepCol, ArrayList<Pair<Integer,ArrayList<String>>> reqStrList, boolean isConcatenated, boolean includeDecoy, 
			boolean hasHeader, int dbCol, String decoyPrefix,
			float fdrThreshold, float pepFDRThreshold, File outputFile)
	{
		TargetDecoyAnalysis tda;
		TSVPSMSet target, decoy;
		if(dbCol >= 0)	// both target and decoy are in the same file
		{
			target = new TSVPSMSet(targetFile, delimeter, hasHeader, scoreCol, isGreaterBetter, specFileCol, specIndexCol, pepCol, reqStrList).decoy(dbCol, decoyPrefix, true);
			target.read();
			decoy = new TSVPSMSet(targetFile, delimeter, hasHeader, scoreCol, isGreaterBetter, specFileCol, specIndexCol, pepCol, reqStrList).decoy(dbCol, decoyPrefix, false);
			decoy.read();
		}
		else
		{
			target = new TSVPSMSet(targetFile, delimeter, hasHeader, scoreCol, isGreaterBetter, specFileCol, specIndexCol, pepCol, reqStrList);
			target.read();
			decoy = new TSVPSMSet(decoyFile, delimeter, hasHeader, scoreCol, isGreaterBetter, specFileCol, specIndexCol, pepCol, reqStrList);
			decoy.read();
		}
		tda = new TargetDecoyAnalysis(target, decoy);

		
		PrintStream out = null;
		if(outputFile != null)
			try {
				out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		else
			out = System.out;
		
		target.writeResults(tda, out, fdrThreshold, pepFDRThreshold);
		if(includeDecoy)
			decoy.writeResults(tda, out, fdrThreshold, pepFDRThreshold);
		
		if(out != System.out)
			out.close();
	}
}
