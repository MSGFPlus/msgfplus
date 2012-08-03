package edu.ucsd.msjava.msdictionary;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;

public class GenomeLocator {
	
	File[] genomeFiles;
	RandomAccessFile[] rafs;
	
	// genome can be either a fasta file or a directory containing multiple fasta files.
	public GenomeLocator(String genomePath)
	{
		File path = new File(genomePath);
		
		class FastaFilter implements FileFilter {
			public boolean accept(File pathname) {
				String name = pathname.getName();
				if(name.substring(name.lastIndexOf('.')+1).equalsIgnoreCase("fasta"))
					return true;
				return false;
			}
		}
		if(path.isDirectory())
			genomeFiles = path.listFiles(new FastaFilter());
		else
		{
			genomeFiles = new File[1];
			genomeFiles[0] = path;
		}
		
		assert(genomeFiles.length > 0);
		rafs = new RandomAccessFile[genomeFiles.length];
		int i = 0;
		for(File f : genomeFiles)
		{
			try {
				rafs[i++] = new RandomAccessFile(f, "r");
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
	}
	
	public String getGenome(String annotation, long filePos, int shift)
	{
		StringBuffer segment = new StringBuffer();
		for(RandomAccessFile raf : rafs)
		{
			StringBuffer s = new StringBuffer();
			try {
				raf.seek(filePos-shift);
				char residue;
				while(true)
				{
					residue = (char)raf.readByte();
					if(residue == 'A' || residue == 'T' || residue == 'C' || residue == 'G')
						s.append(residue);
					else if(residue == '\r' || residue == '\n')
						continue;
					else
						break;
				}
				System.out.println(s);
				System.out.println(GenomeTranslator.translate(s.toString(), shift));
				System.out.println(GenomeTranslator.translateReverseComplement(s.toString(), shift));
				System.out.println();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return segment.toString();
	}
	
	public static void main(String argv[])
	{
		GenomeLocator loc = new GenomeLocator(System.getProperty("user.home")+"/Research/Data/HumanGenome/splitted");
		loc.getGenome("14", 146908676L, 1);
	}
}
