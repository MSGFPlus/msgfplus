package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.parser.BufferedLineReader;


public class SplitFasta {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
		{
			System.out.println("usage: java SplitFasta *.fasta splitNum");
			System.exit(0);
		}
		if(!argv[0].contains("fa"))
		{
			System.out.println(argv[0] + " is not a fasta format");
			System.exit(0);
		}
		split(argv[0], Integer.parseInt(argv[1]));
		System.out.println("Done");
	}
	
	public static void split(String fileName, int splitNum) throws Exception
	{
		File file = new File(fileName);
		long fileSize = file.length();
		long splitSize = fileSize / splitNum;
		
		PrintStream[] out = new PrintStream[splitNum];
		for(int i=0; i<out.length; i++)
		{
			String name = file.getPath().substring(0, file.getPath().lastIndexOf('.')) + "_" + i + ".fasta";
			out[i] = new PrintStream(new BufferedOutputStream(new FileOutputStream(name)));
		}
		
		int fileNum = -1;
		long curSize = Long.MAX_VALUE;
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith(">"))
			{
				if(curSize >= splitSize)
				{
					fileNum++;
					curSize = 0;
				}
			}
			out[fileNum].println(s);
			curSize += s.length()+1;
		}
		
		for(int i=0; i<out.length; i++)
			out[i].close();
	}
}
