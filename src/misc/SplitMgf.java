package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import parser.BufferedLineReader;

public class SplitMgf {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit("Illegal parameters!");
		
		File mgfFile = new File(argv[0]);
		if(!mgfFile.exists())
			printUsageAndExit(argv[0] + " does not exist!");
		String fileName = mgfFile.getName();
		String ext = fileName.substring(fileName.lastIndexOf('.')+1);
		if(!ext.equalsIgnoreCase("mgf"))
			printUsageAndExit(argv[0] + ": Illegal file extension!");
		
		int numParts = Integer.parseInt(argv[1]);
		split(mgfFile, numParts);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.err.println(message);
		System.out.println("usage: java SplitMgf mgfFile(*.mgf) numParts");
		System.exit(-1);
	}
	
	public static void split(File mgfFile, int numParts) throws Exception
	{
		BufferedLineReader in = new BufferedLineReader(mgfFile.getPath());
		String s;
		
		// count the number of spectra
		int numSpecs = 0;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("BEGIN"))
				numSpecs++;
		}
		
		int numSpecsPerPart = numSpecs/numParts;
		
		PrintStream[] out = new PrintStream[numParts];
		for(int i=0; i<out.length; i++)
		{
			String name = mgfFile.getPath().substring(0, mgfFile.getPath().lastIndexOf('.')) + "_Part" + i + ".mgf";
			out[i] = new PrintStream(new BufferedOutputStream(new FileOutputStream(name)));
		}

		int curSize = 0;
		int fileNum = 0;
		in = new BufferedLineReader(mgfFile.getPath());
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("BEGIN"))
				curSize++;
			if(curSize > numSpecsPerPart)
			{
				if(fileNum < numParts-1)
					fileNum++;
				curSize = 0;
			}
			out[fileNum].println(s);
		}
		
		for(int i=0; i<out.length; i++)
			out[i].close();
	}
}
