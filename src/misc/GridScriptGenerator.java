package misc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

public class GridScriptGenerator {
	/*
	public static final String BASE_DIR = "/scratchbig/sak008/";
	public static final String JAVA_PATH = BASE_DIR+"jre1.6.0_15/bin/java";
	public static final String SPEC_DIR = BASE_DIR+"spectra/";
	public static final String DB_DIR = BASE_DIR+"databases/";
	public static final String PARAM_DIR = BASE_DIR+"msgf/";
	public static final String PARAM_NAME = "ionstat_CID_LysN.txt";
	public static final String BIN_PATH = BASE_DIR+"msgf/MSGFDBSearch.jar";
	public static final String IDENTIFIER = "CIDLysN";
	public static final String OUTPUT_DIR = System.getProperty("user.home")+"/Research/MSGF2D/scripts/";
	*/
	
	public static final String BASE_DIR = "/data/sangtae/MSGF2D/";
	public static final String JAVA_PATH = "java";
	public static final String SPEC_DIR = BASE_DIR+"spectra/";
	public static final String DB_DIR = BASE_DIR+"databases/";
	public static final String PARAM_DIR = BASE_DIR+"msgf/";
//	public static final String PARAM_NAME = "ionstat_ETD_tryp.txt";
	public static final String BIN_PATH = BASE_DIR+"msgf/MSGFDBSearchSum.jar";
	public static final String IDENTIFIER = "LysNSum";
	public static final String OUTPUT_DIR = System.getProperty("user.home")+"/Research/MSGF2D/scripts/";
	
	public static void main(String argv[])
	{
		makePairedSearch();
		makeDriver();
		System.out.println("Done");
	}
	
	public static void makeDriver()
	{
		boolean memorySpecified = true;
		String memory = "3G";
		PrintStream out = null; 
		try {
			out = new PrintStream(new File(OUTPUT_DIR+
					"run"+IDENTIFIER+"Search"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		out.println("#!/bin/bash");
		for(int specNum=0; specNum<2; specNum++)
			for(int dbNum=0; dbNum<2; dbNum++)
				for(int seqNum=0; seqNum<25; seqNum++)
					out.println("qsub" + (memorySpecified ? " -l h_vmem="+memory : "") + " scripts/"+IDENTIFIER+seqNum+"_"+specNum+"_"+dbNum+".sh");
		out.close();
	}
	
	public static void makePairedSearch()
	{
		for(int specNum=0; specNum<2; specNum++)
			for(int dbNum=0; dbNum<2; dbNum++)
				for(int seqNum=0; seqNum<25; seqNum++)
					makeSH(seqNum, specNum, dbNum);
	}
	
	public static void makeSH(int seqNum, int specNum, int dbNum)
	{
		int startSpecNum = seqNum*313+1;
		int endSpecNum = startSpecNum+312;
//		if(seqNum == 24)
//			endSpecNum += 2;
		
		String specFileNameCID = null;
		String specFileNameETD = null;
		if(specNum == 0)
		{
			specFileNameCID = SPEC_DIR+"CID_paired_01_LysN.mgf";
			specFileNameETD = SPEC_DIR+"ETD_paired_01_LysN.mgf";
		}
		else
		{
			specFileNameCID = SPEC_DIR+"CID_paired_03_LysN.mgf";
			specFileNameETD = SPEC_DIR+"ETD_paired_03_LysN.mgf";
		}
		
		String dbFileName = null;
		if(dbNum == 0)
			dbFileName = DB_DIR+"uniprot_sprot.fasta";
		else if(dbNum == 1)
			dbFileName = DB_DIR+"reverse_uniprot_sprot.fasta";
		
		
		PrintStream out = null; 
		try {
			out = new PrintStream(OUTPUT_DIR+
					IDENTIFIER+seqNum+"_"+specNum+"_"+dbNum+".sh");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		out.println("#!/bin/bash");
		out.println("#");
		out.println("#$ -cwd");
		out.println("#$ -j y");
		out.println("#$ -S /bin/bash");
		out.println("#");
		
		
		// paired
		out.println(JAVA_PATH +
				" -Xmx1800M " +
				"-jar " +
				BIN_PATH + " " +
				specFileNameCID + " " +
				specFileNameETD + " " +
				dbFileName + " " +
				PARAM_DIR + " " +
				startSpecNum + " " +
				endSpecNum + " " +
				"> res"+IDENTIFIER+seqNum+"_"+specNum+"_"+dbNum+".txt");
		
		/**
		out.println(JAVA_PATH +
				" -Xmx1800M " +
				"-jar " +
				BIN_PATH + " " +
				PARAM_DIR + PARAM_NAME + " " +
				dbFileName + " " +
//				specFileNameCID + " " +
				specFileNameETD + " " +
				startSpecNum + " " +
				endSpecNum + " " +
				"> res"+IDENTIFIER+seqNum+"_"+specNum+"_"+dbNum+".txt");
		*/
		out.close();
	}
	
	public static void makeHumanSearch()
	{
		String dbDir = "/Users/sangtaekim/Research/Data/Human";
		File dir = new File(dbDir);
		if(!dir.isDirectory())
			System.exit(0);
		File[] dbFiles = dir.listFiles();
		int index = 0;
		for(File file : dbFiles)
		{
			index++;
			String fileName = file.getName();
			if(!fileName.startsWith("HS_"))
				continue;
			int fileNum = Integer.parseInt(fileName.substring(fileName.indexOf("frame")+5, fileName.lastIndexOf('.')));
			if(fileNum < 40 || fileNum >= 100)

				continue;
			PrintStream out = null; 
			try {
				out = new PrintStream(new File("scripts/hSearch"+fileNum+".sh"));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			
			out.println("#!/bin/bash");
			out.println("#");
			out.println("#$ -cwd");
			out.println("#$ -j y");
			out.println("#$ -S /bin/bash");
			out.println("#");
			out.println("/usr/java/jdk1.5.0_07/bin/java -Xmx1800M compatiblepeptides/HybridSearch -d dictionaries -db human/" +
					fileName + " > searchResult/human" + fileNum + ".txt");
			out.close();
		}
		
	}
}
