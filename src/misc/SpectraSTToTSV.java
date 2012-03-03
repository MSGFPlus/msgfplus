package misc;

import java.io.File;

import parser.BufferedLineReader;

public class SpectraSTToTSV {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1)
			printUsageAndExit(null);
		File resultFile = new File(argv[0]);
		if(!resultFile.exists())
			printUsageAndExit("File does not exist!");
		convert(resultFile);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.err.println(message);
		System.out.println("usage: java SpectraSTToTSV SpectraSTResult");
		System.exit(-1);
	}
	
	public static void convert(File resultFile) throws Exception
	{
		String s;
		BufferedLineReader in = new BufferedLineReader(resultFile.getPath());
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\\s+");
			if(token.length != 17)
			{
				if(token.length == 16)
				{
					String[] newToken = new String[17];
					String peptideCol = token[2];
					String pep = peptideCol.substring(0, peptideCol.lastIndexOf('0'));
					String score = peptideCol.substring(peptideCol.lastIndexOf('0'));
					
					for(int i=0; i<token.length; i++)
					{
						if(i < 2)
							newToken[i] = token[i];
						if(i == 2)
						{
							newToken[i] = pep;
							newToken[i+1] = score;
						}
						else
						{
							newToken[i+1] = token[i];
						}
					}
					token = newToken;
				}
				else
				{
					System.err.println("Illegal output: " + s);
					System.err.println("Token length: " + token.length);
					System.exit(-1);
				}
			}
			
			System.out.print(token[0]);
			for(int i=1; i<token.length; i++)
				System.out.print("\t"+token[i]);
			System.out.println();
		}
	}
}
