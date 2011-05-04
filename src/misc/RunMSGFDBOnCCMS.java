package misc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

import msgf.Tolerance;
import msutil.Enzyme;

public class RunMSGFDBOnCCMS {
	public static final String JAVA_PATH = "/usr/bin/java";

	public static void main(String argv[])
	{
		File dir = null;
		int[] methods = {0};
		int enzymeNum = 1;
		int msgfThreshold = 0;
		int numMatches = 10;
		int heapSize = 5;
		Tolerance tolerance = null;
		File javaExe = new File(JAVA_PATH);
		File aaSetFile = null;
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
			else if(argv[i].equalsIgnoreCase("-t"))
			{
				tolerance = Tolerance.parseToleranceStr(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-m"))
			{
				String[] token = argv[i+1].split(",");
				methods = new int[token.length];
				for(int tokenNum=0; tokenNum<token.length; tokenNum++)
				{
					methods[tokenNum] = Integer.parseInt(token[tokenNum]);
					if(methods[tokenNum] < 0 || methods[tokenNum] > 3)
						printUsageAndExit("Wrong argument: " + argv[i] + " " + argv[i+1]);
				}
			}
			else if(argv[i].equalsIgnoreCase("-e"))
			{
				enzymeNum = Integer.parseInt(argv[i+1]);
				if(enzymeNum < 0 || enzymeNum > 7)
					printUsageAndExit("Wrong argument: " + argv[i] + " " + argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-j"))
			{
				javaExe = new File(argv[i+1]);
				if(!javaExe.exists())
					printUsageAndExit("Wrong java path: " + argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-h"))
			{
				heapSize = Integer.parseInt(argv[i+1]);
				if(heapSize < 1 || heapSize > 30)
					printUsageAndExit("Wrong argument: " + argv[i] + " " + argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-f"))
			{
				msgfThreshold = Integer.parseInt(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-n"))
			{
				numMatches = Integer.parseInt(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-aaSet"))
			{
				aaSetFile = new File(argv[i+1]);
				if(!aaSetFile.exists())
					printUsageAndExit("Wrong aaSet file: " + argv[i+1]);
			}
		}
		if(dir == null)
			printUsageAndExit("-d option is missing!");
		if(tolerance == null)
			printUsageAndExit("-t option is missing!");
		makeScripts(javaExe, dir, methods, enzymeNum, tolerance, msgfThreshold, numMatches, heapSize, aaSetFile);
		makeDriver(dir, heapSize);
		System.out.println("Done");
	}
	
	public static void printUsageAndExit(String message)
	{
		System.err.println(message);
		System.out.println("usage: java RunMSGFOnCCMS \n" +
				"\t-d dir\n" +
				"\t-t tolerance (e.g. 50.0ppm, 0.5Da)\n" +
				"\t[-m methods (0:all,1:cid,2:etd,3:cid/etd e.g. 1,2,3 or 2)\n" +
				"\t[-e Enzyme 1/3/4] (1: Trypsin (default), 3: Lys-C, 4: Lys-N)\n" +
				"\t[-f msgfScoreThreshold (default: 0)]\n" +
				"\t[-h heapSizeAsGB (Default: 5)]\n" +
				"\t[-n numMatchesPerSpec (Default: 10)]\n" +
				"\t[-aaSet aaSetFilePath (Default: Standard with fixed C+57)]\n" +
				"\t[-j javaPath]");
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
	
	public static void makeScripts(File javaExe, File dir, int[] methods, int enzymeNum, 
			Tolerance tolerance, int msgfThreshold, int numMatches, int heapSize, File aaSetFile)
	{
		File specDir = new File(dir.getPath()+File.separator+"spectra");
		File databaseDir = new File(dir.getPath()+File.separator+"database");
		File scriptDir = new File(dir.getPath()+File.separator+"scripts");
		File excutablePath = new File(dir.getPath()+File.separator+"MSGFDB.jar");
		
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
}
