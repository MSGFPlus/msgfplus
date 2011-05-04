package misc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

public class RunOMSSAOnCCMS {
	public static final String BASE_DIR = "/data/sangtae/omssa-2.1.7.linux/";
	public static final String GRID_SPEC_DIR = "/data/sangtae/HeckWhole/spectra/";
	public static final String SPEC_DIR = "/home/sangtaekim/Research/Data/HeckWhole/Spectra/";
	public static final String BIN_PATH = BASE_DIR+"omssacl";
	public static final String SCRIPT_DIR = "/home/sangtaekim/Research/Data/HeckWhole/scripts/";
	public static final String RESULT_DIR = BASE_DIR+"results/";
	public static final String TARGET_DB = "/data/sangtae/Research/Data/SProt/uniprot_sprot.fasta";
	public static final String DECOY_DB = "/data/sangtae/Research/Data/SProt/reverse_uniprot_sprot.fasta";
	public static final String TOLERANCE = "50ppm";

	public static void main(String argv[])
	{
		makeScripts();
		makeDriver();
	}
	public static void makeDriver()
	{
		boolean memorySpecified = true;
		String memory = "3G";
		PrintStream out = null; 
		try {
			out = new PrintStream(new File(SCRIPT_DIR+
					"runOMSSA"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		out.println("#!/bin/bash");
		
		File scriptDir = new File(SCRIPT_DIR);
		for(File f : scriptDir.listFiles())
		{
			if(f.getName().endsWith(".sh"))
			{
				out.println("qsub" + (memorySpecified ? " -l h_vmem="+memory : "") + " " + BASE_DIR+"scripts/"+f.getName());
			}
		}
		out.close();
	}
	
	public static void makeScripts()
	{
		File specDir = new File(SPEC_DIR);
		for(File f : specDir.listFiles())
		{
			String fileName = f.getName();
			if(!fileName.endsWith(".mgf"))
				continue;
			
			boolean isLysN = false;
			if(fileName.contains("LysN"))
				isLysN = true;
			
			boolean isETD = false;
			if(fileName.contains("ETD"))
				isETD = true;
			
			for(int i=0; i<2; i++)
			{
				File dbFile;
				if(i==0)
					dbFile = new File(TARGET_DB);
				else
					dbFile = new File(DECOY_DB);
				PrintStream out = null; 
				try {
					out = new PrintStream(new File(SCRIPT_DIR+
							"omssa"+fileName.substring(fileName.indexOf('_')+4, fileName.lastIndexOf('.'))+"_"+i+".sh"));
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
				out.println("#!/bin/bash");
				out.println("#");
				out.println("#$ -cwd");
				out.println("#$ -j y");
				out.println("#$ -S /bin/bash");
				out.println("#");
				out.println(
						BIN_PATH
						+ " -fm " + GRID_SPEC_DIR+fileName
						+ " -oc " + "OMSSA_" + fileName.substring(0, fileName.lastIndexOf('.')) + "_" + i + ".csv"
						+ " -to 0.5"
						+ " -te 0.05"
						+ " -tez 0"
						+ " -hl 10"
						+ " -he 100"
						+ " -mf 3"
						+ " -d " + dbFile.getPath()
						+ (!isETD ? "" : " -cp 1 -i 2,5")
						+ (!isLysN ? "" : " -e 21")
						+ " -ht 10"
						+ " -zcc 1"
						+ " -zl 2"
						+ " -zh 8"
						+ " -v 3"
						+ " -tem 0"
						+ " -tom 0"
						+ " -sb1 0"
				);		
			}
			
		}		
	}
}

