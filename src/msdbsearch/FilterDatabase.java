package msdbsearch;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;

import parser.BufferedLineReader;

import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class FilterDatabase {
	public static void main(String argv[]) throws Exception
	{
		File dbFile = null;
		File resultFile = null;
		int pepColumn = -1;
		String delimeter = "\t";
		
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length)
				printUsageAndExit("Illegal parameters");
			else if(argv[i].equalsIgnoreCase("-d"))
			{
				dbFile = new File(argv[i+1]);
				if(!dbFile.exists())
					printUsageAndExit(argv[i+1] + " doesn't exist");
			}
			else if(argv[i].equalsIgnoreCase("-r"))
			{
				resultFile = new File(argv[i+1]);
				if(!resultFile.exists())
					printUsageAndExit(argv[i+1] + " doesn't exist");
			}
			else if(argv[i].equalsIgnoreCase("-p"))
			{
				pepColumn = Integer.parseInt(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-delim"))
			{
				delimeter = argv[i+1];
			}
		}
		
		if(dbFile == null || !dbFile.exists())
			printUsageAndExit("Illegal dbFile!");
		if(resultFile == null || !resultFile.exists())
			printUsageAndExit("Illegal resultFile!");
		if(pepColumn < 0)
			printUsageAndExit("Illegal pepColumn!");
		filterDatabase(dbFile, resultFile, pepColumn, delimeter);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.out.println(message);
		System.out.println("usage: java FilterDatabase\n" +
				"\t-d database(*.fasta)\n" +
				"\t-r searchResult\n" +
				"\t-p pepColumn\n" +
				"\t[-delim delimeter] (default: \\t)");
		System.exit(-1);
	}
	
	public static void filterDatabase(File dbFile, File resultFile, int pepColumn, String delimeter) throws Exception
	{
	      SuffixArraySequence sequence = new SuffixArraySequence(dbFile.getPath());
	      SuffixArray sa = new SuffixArray(sequence);
	      HashSet<String> matchedEntrySet = new HashSet<String>();
	      
	      String s;
	      BufferedLineReader in = new BufferedLineReader(resultFile.getPath());
	      in.readLine();	// header
	      while((s=in.readLine()) != null)
	      {
	    	  String[] token = s.split(delimeter);
	    	  if(token.length < pepColumn)
	    		  continue;
	    	  String annotation = token[pepColumn];
	    	  String pepStr = annotation;
	    	  if(annotation.matches("[A-Z]\\.[A-Z]+\\.[A-Z]"))
	    		  pepStr = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
	    	  ArrayList<String> matchedEntries = sa.getAllMatchingAnnotations(pepStr);
	    	  ArrayList<String> matchedProtSeq = sa.getAllMatchingEntries(pepStr);
	    	  for(int i=0; i<matchedEntries.size(); i++)
	    	  {
	    		  String entry = matchedEntries.get(i);
	    		  if(!matchedEntrySet.contains(entry))
				  {
	    			  matchedEntrySet.add(entry);
	    			  System.out.println(">"+entry);
	    			  System.out.println(matchedProtSeq.get(i));
				  }
	    	  }
	      }
	}
	
}
