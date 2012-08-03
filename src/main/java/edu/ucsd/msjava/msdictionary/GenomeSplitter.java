package edu.ucsd.msjava.msdictionary;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;

import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.LineReader;



public class GenomeSplitter {
	private final String fileName;
	private Hashtable<String, Long> lengthTable;
	private long sum;
	private long[] subsetSize;
	private StringHashSet[] subsetAnnotation;
	
	private static class StringHashSet extends HashSet<String>{
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;}
	
	public GenomeSplitter(String fileName)
	{
		this.fileName = fileName;
	}
	
	public long getSum()	{ return sum; }
	public long[] getSubsetSizeArr()	{ return subsetSize; }
	
	private GenomeSplitter scan()
	{
		LineReader reader = null;
		try {
			reader = new BufferedLineReader(fileName);
		} catch (Exception e) {
			e.printStackTrace();
		}

		lengthTable = new Hashtable<String, Long>();
		long sum = 0;
		String annotation = null;
		String s;
		long chrSize = 0;
		
		while((s = reader.readLine()) != null)
		{
			if(s.length() == 0)
				continue;
			if(s.startsWith(">"))
			{
				if(annotation != null)
				{
					lengthTable.put(annotation, chrSize);
					sum += chrSize;
				}
				annotation = s.substring(1);
				chrSize = 0;
			}
			else
			{
				for(int i=0; i<s.length(); i++)
				{
					char c = s.charAt(i); 
					if(c == 'A' || c == 'T' || c == 'C' || c == 'G')
						chrSize++;
					else 
						continue;
				}
			}
		}		
		if(annotation != null)
		{
			lengthTable.put(annotation, chrSize);
			sum += chrSize;
		}
		this.sum = sum;
		return this;
	}
	
	public GenomeSplitter split(int numSubsets)
	{
		if(lengthTable == null)
			return this;
		subsetAnnotation = new StringHashSet[numSubsets];
		
		long lengthOfSubset = sum/numSubsets;
		float allowedError = 0.1f;
		ArrayList<String> annotationList = new ArrayList<String>(lengthTable.keySet());
		
		subsetSize = new long[numSubsets];
		int[] subsetBoundaries = new int[numSubsets];
		
		
		boolean success = false;
		while(!success)
		{
			Collections.shuffle(annotationList);
			int subsetIndex = 0;
			long size = 0;
			long sizeSoFar = 0;
			for(int i=0; i<annotationList.size(); i++)
			{
				if(subsetIndex == numSubsets-1)
				{
					long remainingSize = sum - sizeSoFar;
					if(remainingSize >= (long)(lengthOfSubset*(1-allowedError)) && 
							remainingSize <= (long)(lengthOfSubset*(1+allowedError)))
					{
						subsetSize[numSubsets-1] = remainingSize;
						subsetBoundaries[numSubsets-1] = annotationList.size()-1;
						success = true;
					}
					break;
				}
				else
				{
					String annotation = annotationList.get(i);
					long len = lengthTable.get(annotation);
					if(size+len < (int)(lengthOfSubset*(1-allowedError)))
					{
						size += len;
					}
					else if(size+len >= (long)(lengthOfSubset*(1-allowedError)) && 
							size+len <= (long)(lengthOfSubset*(1+allowedError)))
					{
						subsetSize[subsetIndex] = size+len;
						subsetBoundaries[subsetIndex] = i;
						sizeSoFar += size+len;
						subsetIndex++;
						size = 0;
					}
					else
						break;
				}
			}
		}
		int startIndex = 0;
		for(int i=0; i<subsetSize.length; i++)
		{
			subsetAnnotation[i] = new StringHashSet();
			for(int j=startIndex; j<=subsetBoundaries[i]; j++)
				subsetAnnotation[i].add(annotationList.get(j));
			startIndex = subsetBoundaries[i]+1;
		}
		return this;
	}
	
	// saved as prefix0.fasta, prefix1.fasta,...
	public void save(String prefix)
	{
		if(subsetAnnotation == null)
			return;
		
		LineReader reader = null;
		try {
			reader = new BufferedLineReader(fileName);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		PrintStream[] out = new PrintStream[subsetAnnotation.length];
		for(int i=0; i<out.length; i++)
		{
			try {
				out[i] = new PrintStream(new BufferedOutputStream(new FileOutputStream(prefix+i+".fasta")));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		// scan and save
		String s;
		PrintStream curOut = null;
		
		while((s = reader.readLine()) != null)
		{
			if(s.length() == 0)
				continue;
			if(s.startsWith(">"))
			{
				curOut = null;
				String annotation = s.substring(1);
				int index = -1;
				for(StringHashSet hash : subsetAnnotation)
				{
					index++;
					if(hash.contains(annotation))
						curOut = out[index];
				}
			}
			curOut.println(s);
		}
		for(PrintStream o : out)
			o.close();
	}
	
	public void printLengthDistribution()
	{
		if(lengthTable == null)
			return;
		ArrayList<String> annotationList = new ArrayList<String>(lengthTable.keySet());
		Collections.sort(annotationList);
		
		for(String annotation : annotationList)
		{
			System.out.println(annotation + "\t" +lengthTable.get(annotation));
		}
	}
	
	public static void main(String argv[])
	{
		if(argv.length != 2)
		{
			printUsage();
			System.exit(0);
		}
		else
		{
			String source = argv[0];
			if(!source.contains(".fa"))
			{
				System.out.println("Source is not the fasta format.");
				System.exit(0);
			}
			
			String targetPrefix = source.substring(0, source.lastIndexOf('.')+1);
			int numSubsets = Integer.parseInt(argv[1]);
			
			GenomeSplitter splitter = new GenomeSplitter(argv[0]);
			splitter.scan().split(numSubsets).save(targetPrefix);
			System.out.println("Total unmasked length: " + splitter.getSum());
			int subsetIndex = -1;
			for(long s : splitter.getSubsetSizeArr())
				System.out.println("Subset " + (++subsetIndex) + ": " + s);
			System.out.println("Done");
		}
	}
	
	public static void printUsage()
	{
		System.out.println("usage: java -Xmx(HeapSize) GenomeSplitter source(fasta) #files\n" +
				"(Example: java GenomeSplitter test.fasta 4)");
	}
}
