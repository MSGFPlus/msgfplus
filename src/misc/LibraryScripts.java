package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import parser.BufferedLineReader;

public class LibraryScripts {
	public static void main(String argv[]) throws Exception
	{
		convert();
	}
	
	public static void convert() throws Exception
	{
		File inputFile = new File("/Users/sangtaekim/Research/Data/SpecLib/human_target.mgf");
		File outputFile = new File("/Users/sangtaekim/Research/Data/SpecLib/human_target_annotated.mgf");
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		String s;
		BufferedLineReader in = new BufferedLineReader(inputFile.getPath());
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SEQ="))
			{
				StringBuffer buf = new StringBuffer();
				String[] token = s.split("[(,)]");
				for(int i=0; i<token.length; i++)
				{
					if(token[i].length() == 0)
						continue;
					if(Character.isDigit(token[i].charAt(0)))
						buf.append("+");
					buf.append(token[i]);
				}
				out.println(buf.toString());
			}
			else if(s.startsWith("PRECURSOR="))
			{
				out.println("PEPMASS="+s.substring("PRECURSOR=".length()));
			}
			else
				out.println(s);
		}
		in.close();
		out.close();
		System.out.println("Done");
	}
}
