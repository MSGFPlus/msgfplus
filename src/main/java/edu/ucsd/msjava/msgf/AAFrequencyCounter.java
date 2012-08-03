package edu.ucsd.msjava.msgf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

public class AAFrequencyCounter {
	Histogram<String> frequencyTable;
	int nMer;
	int sizeNMer;

	public AAFrequencyCounter()
	{
		frequencyTable = new Histogram<String>();
		sizeNMer = 0;
	}
	
	public void setNMer(int nMer)
	{
		this.nMer = nMer;
	}
	
	public void readFromFreqFile(String fileName)
	{
		BufferedReader in = null;
		try {
			in = new BufferedReader(new FileReader(fileName));
			String s;
			
			s = in.readLine();
			String[] token = s.split("\t");
			assert(token[0].equalsIgnoreCase("n"));
			this.nMer = Integer.parseInt(token[1]);
			
			s = in.readLine();
			token = s.split("\t");
			assert(token[0].equalsIgnoreCase("size"));
			this.sizeNMer = Integer.parseInt(token[1]);
			
			while((s = in.readLine()) != null)
			{
				token = s.split("\t");
				assert(token.length == 2);
				frequencyTable.put(token[0], Integer.parseInt(token[1]));
			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	
	public void readFromFasta(String fileName) 
	{
		BufferedReader in = null;
		try {
			in = new BufferedReader(new FileReader(fileName));
			String s;
			while((s = in.readLine()) != null)
			{
				if(s.startsWith(">"))
					continue;
				StringBuffer buf = new StringBuffer();
				for(int i=0; i<s.length(); i++)
				{
					if(i >= nMer)
					{
						frequencyTable.add(buf.toString());
						sizeNMer++;
						buf.deleteCharAt(0);
					}
					buf.append(s.charAt(i));
				}
			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	
	public static float getRandomFrequency(String str)
	{
		float uniFreq = 0.05f;
		int numLI = 0;
		for(int i=0; i<str.length(); i++)
			if(str.charAt(i) == 'L' || str.charAt(i) == 'I')
				numLI++;
		return (float)(Math.pow(2, numLI)*Math.pow(uniFreq, str.length()));
	}
	
	public float getFrequency(String str)
	{
		ArrayList<String> strSet = new ArrayList<String>();
		strSet.add(str);
		for(int i=0; i<str.length(); i++)
		{
			char c = str.charAt(i);
			if(c == 'L')
			{
				int size = strSet.size();
				for(int j=0; j<size; j++)
				{
					String s = strSet.get(j);
					strSet.add(s.substring(0,i)+"I"+s.substring(i+1));
				}
			}
			else if(c == 'I')
			{
				int size = strSet.size();
				for(int j=0; j<size; j++)
				{
					String s = strSet.get(j);
					strSet.add(s.substring(0,i)+"L"+s.substring(i+1));
				}
			}
		}
		int occ = 0;
		for(String s : strSet)
			occ += getOccurrence(s);
		return occ / (float)sizeNMer;
	}
	
	public int getOccurrence(String str)
	{
		Integer occ = frequencyTable.get(str);
		if(occ == null)
			return 0;
		else
			return occ;
	}
	
	public static void main(String argv[])
	{
		System.out.println(getRandomFrequency("AAA"));
//		generate(3);
		/*
		AAFrequencyCounter counter = new AAFrequencyCounter();
		counter.readFromFreqFile("/home/sangtaekim/Research/Data/AAFrequency/SProt_2mer.txt");
		counter.frequencyTable.printSorted();
		*/
	}
	public static void generate(int nMer)
	{
		AAFrequencyCounter counter = new AAFrequencyCounter();
		counter.setNMer(nMer);
//		counter.readFromFasta("/home/sangtaekim/Research/Data/SProt/uniprot_sprot.fasta");
		counter.readFromFasta("/home/sangtaekim/Research/Data/EColiDB/Ecol_protein_formatted.fasta");

		System.out.println("n\t" + nMer);
		System.out.println("size\t" + counter.sizeNMer);
		String allAA = "GASPVTCLINDQKEMHFRYW";
		
		if(nMer == 1)
		{
			for(int i=0; i<allAA.length(); i++)
			{
				char c = allAA.charAt(i);
				System.out.println(c + "\t" + counter.getOccurrence(String.valueOf(c)));
			}
			
		}
		else if(nMer == 2)
		{
			for(int i=0; i<allAA.length(); i++)
			{
				for(int j=0; j<allAA.length(); j++)
				{
					String s = ""+ allAA.charAt(i) + allAA.charAt(j);
					System.out.println(s + "\t" + counter.getOccurrence(s));
				}
			}
		}
		else if(nMer == 3)
		{
			for(int i=0; i<allAA.length(); i++)
			{
				for(int j=0; j<allAA.length(); j++)
				{
					for(int k=0; k<allAA.length(); k++)
					{
						String s = ""+ allAA.charAt(i) + allAA.charAt(j) + allAA.charAt(k);
						System.out.println(s + "\t" + counter.getOccurrence(s));
					}
				}
			}
		}
	}
}
