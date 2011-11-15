package misc;

import java.io.File;

public class VennDiagram {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit();
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("usage: java VennDiagram result1 result2 [threshold]");
		System.exit(-1);
	}
	
	public static void vennDiagram(File result1, File result2) throws Exception
	{
		
	}
}
