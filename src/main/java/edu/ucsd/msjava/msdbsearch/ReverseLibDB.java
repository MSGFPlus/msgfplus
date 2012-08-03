package edu.ucsd.msjava.msdbsearch;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

public class ReverseLibDB {

	public static void main(String argv[])
	{
		if(argv.length != 2)
			printUsageAndExit();
		
		String ext1 = argv[0].substring(argv[0].lastIndexOf('.')+1);
		String ext2 = argv[1].substring(argv[1].lastIndexOf('.')+1);
		if(!ext1.equalsIgnoreCase("fasta") || !ext2.equalsIgnoreCase("fasta"))
		{
			System.out.println(ext1 + "," + ext2);
			printUsageAndExit();
		}
		reverseDB(argv[0], argv[1], true);
		
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("usage: java ReverseLibDB input(fasta) output(fasta)");
		System.exit(0);
	}
	
	public static boolean reverseDB(String inFileName, String outFileName, boolean concat)
	{
		BufferedReader in = null;
		PrintStream out = null;
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFileName)));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		
		String s;
		if(concat)
		{
			try {
				in = new BufferedReader(new FileReader(inFileName));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			try {
				while((s = in.readLine()) != null)
				{
					out.println(s);
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		try {
			in.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		try {
			in = new BufferedReader(new FileReader(inFileName));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		StringBuffer peptide = null;
		String annotation = null;
		try {
			while((s = in.readLine()) != null)
			{
				if(s.startsWith(">"))	// start of a protein
				{
					if(annotation != null)
					{
						StringBuffer rev = new StringBuffer();
						rev.append(peptide.charAt(0));
						for(int i=peptide.length()-2; i>0; i--)
							rev.append(peptide.charAt(i));
						rev.append(peptide.charAt(peptide.length()-1));
						out.println(">REV_" + annotation);
						out.println(rev.toString().trim());
					}
					annotation = s.substring(1);
					peptide = new StringBuffer();
				}
				else
					peptide.append(s);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		if(peptide != null && annotation != null)
		{
			StringBuffer rev = new StringBuffer();
			rev.append(peptide.charAt(0));
			for(int i=peptide.length()-2; i>0; i--)
				rev.append(peptide.charAt(i));
			rev.append(peptide.charAt(peptide.length()-1));
			out.println(">REV_" + annotation);
			out.println(rev.toString().trim());
		}
		try {
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		out.close();
		
		return true;
	}
	
	public static boolean copyDB(String inFileName, String outFileName)
	{
		BufferedReader in = null;
		PrintStream out = null;
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFileName)));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		
		String s;
		try {
			in = new BufferedReader(new FileReader(inFileName));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			while((s = in.readLine()) != null)
			{
				out.println(s);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		out.flush();
		out.close();
		
		return true;
	}	
}
