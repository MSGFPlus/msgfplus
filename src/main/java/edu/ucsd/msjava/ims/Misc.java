package edu.ucsd.msjava.ims;

import java.io.File;
import java.util.HashSet;

import edu.ucsd.msjava.parser.BufferedLineReader;

public class Misc {
	public static void main(String argv[]) throws Exception
	{
		countPeptides();
	}
	
	public static void countPeptides() throws Exception
	{
		File file = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\IMS\\Sarcopenia\\peptides.txt");
		HashSet<String> pepSet = new HashSet<String>();
		
		String s;
		BufferedLineReader in = new BufferedLineReader(file.getPath());
		while((s=in.readLine()) != null)
		{
			s = s.trim();
			if(s.isEmpty())
				continue;
			pepSet.add(s.substring(s.indexOf('.')+1, s.lastIndexOf('.')));
		}
		System.out.println("NumPeptides: " + pepSet.size());
	}
}
