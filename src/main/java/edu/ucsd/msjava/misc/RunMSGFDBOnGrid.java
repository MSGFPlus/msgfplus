package edu.ucsd.msjava.misc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

import edu.ucsd.msjava.msgf.Tolerance;

public class RunMSGFDBOnGrid {
	/*
	public static void main(String argv[])
	{
		File dir = null;
		int heapSize = 1300;
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length)
				printUsageAndExit("Illegal parameters");
			else if(argv[i].equalsIgnoreCase("-d"))
			{
				dir = new File(argv[i+1]);
				if(!dir.isDirectory())
					printUsageAndExit(argv[i+1] + " is not a directory!");
			}
			else if(argv[i].equalsIgnoreCase("-h"))
			{
				heapSize = Integer.parseInt(argv[i+1]);
				if(heapSize < 100 || heapSize > 4096)
					printUsageAndExit("Wrong parameters (heap size must between 1024 and 8192): " + argv[i] + " " + argv[i+1]);
			}
			else
				printUsageAndExit("Wrong parameters!");
		}
		if(dir == null)
			printUsageAndExit("-d parameter is missing!");
		makeScripts(dir, heapSize);
		makeDriver(dir, heapSize);
		System.out.println("Done");
	}
	
	public static void printUsageAndExit(String message)
	{
		System.err.println(message);
		System.out.println("usage: java RunMSGFOnCCMS\n" +
				"\t-d directory\n" +
				"\t-[h heapSize in MB (default: 1300M)]\n");
		System.exit(-1);
	}
	
	public static void makeDriver(File dir, int heapSize)
	{
		File scriptDir = new File(dir.getPath()+File.separator+"scripts");
		PrintStream out = null; 
		try {
			out = new PrintStream(new File(dir.getPath()+File.separator+
					"runMSGFDB"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		out.println("#!/bin/bash");
		
		for(File f : scriptDir.listFiles())
		{
			if(f.getName().endsWith(".sh"))
			{
				out.println("qsub" + " -l h_vmem="+heapSize+"G,virtual_free="+heapSize+"G " + f.getPath());
			}
		}
		out.flush();
		out.close();
	}
	
	public static void makeScripts(File dir, int heapSize)
	{
		File specDir = null;
		File databaseDir = null;
		File scriptDir = null;
		File excutablePath = null;
		
		for(File f : dir.listFiles())
		{
			
		}
		
		// sanity check
		if(!specDir.isDirectory())
		{
			System.err.println("spectra directory is missing!");
			System.exit(-1);
		}
		if(!databaseDir.isDirectory())
		{
			System.err.println("database directory is missing!");
			System.exit(-1);
		}
		if(!excutablePath.exists())
		{
			System.err.println("MSGFDB executable doesn't exist!");
			System.exit(-1);
		}
		if(!scriptDir.exists())
			scriptDir.mkdir();
		
		for(int method : methods)
		{
			String methodStr = "";
			if(method==1)
				methodStr = "CID";
			else if(method==2)
				methodStr = "ETD";
			else if(method==3)
				methodStr = "Sum";
			for(File f : specDir.listFiles())
			{
				String fileName = f.getName();
				if(!fileName.endsWith(".mzXML") && !fileName.endsWith(".mgf"))
					continue;

				String enzymeName = "";
				if(enzymeNum == 1)
					enzymeName = "Tryp";
				else if(enzymeNum == 3)
					enzymeName = "LysC";
				else if(enzymeNum == 4)
					enzymeName = "LysN";

				int heap = heapSize*1000 - 500;
				int i=-1;
				for(File dbFile : databaseDir.listFiles())
				{
					if(!dbFile.getName().endsWith(".fasta"))
						continue;
					i++;
					String name = methodStr+enzymeName+"_"+fileName.substring(fileName.indexOf('_')+1, fileName.lastIndexOf('.'))+"_"+i;
					PrintStream out = null; 
					try {
						out = new PrintStream(new File(scriptDir.getPath()+File.separator+
								name+".sh"));
					} catch (FileNotFoundException e) {
						e.printStackTrace();
					}
					out.println("#!/bin/bash");
					out.println("#");
					out.println("#$ -cwd");
					out.println("#$ -j y");
					out.println("#$ -S /bin/bash");
					out.println("#");
					out.print(javaExe.getPath() +
							" -Xmx"+heap+"M " +
							"-jar " +
							excutablePath.getPath() + 
							" -s " + f.getPath() +
							" -d " + dbFile.getPath() + 
							" -t " + tolerance.toString() +
							" -e " + enzymeNum +
							" -filter " + msgfThreshold +
							" -m " + method +
							" -n " + numMatches +
							" -o " + "MSGFDB_" + name + ".txt"
							);	
					if(aaSetFile != null)
						out.print(" -aaSet " + aaSetFile.getPath());
					out.println();
					out.flush();
					out.close();
				}
			}			
		}
	}	
	*/
}
