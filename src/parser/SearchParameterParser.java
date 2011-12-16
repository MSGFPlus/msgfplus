package parser;

import java.io.File;

public class SearchParameterParser {
	private boolean spectrumFile = true;
	private boolean database = true;
	private boolean parentMassTolerance = true;

	
	public void readParamDir(File paramDir)
	{
		if(paramDir.exists() && paramDir.isDirectory())
		{
			for(File f : paramDir.listFiles())
			{
				if(f.getName().endsWith("txt"))
					System.out.println(f.getName());
			}
		}
	}
	

}
