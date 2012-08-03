package edu.ucsd.msjava.msdictionary;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.parser.FullyBufferedLineReader;



public class GenomeTranslator {
	private static final int MIN_PROTEIN_SIZE = 6;
	
	private final String fileName;
	public GenomeTranslator(String genomeFileName)
	{
		this.fileName = genomeFileName;
	}

	public void translateAndWriteTo(String outputFileName)
	{
		FullyBufferedLineReader reader = new FullyBufferedLineReader(fileName);
		
		PrintStream out = null;
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		String s;
		String annotation = null;
		StringBuffer segment = null;
		int segmentPos = 0;
		int curPosition = 0;
		//int curLinePosition = 0;
		while((s = reader.readLine()) != null)
		{
			if(s.startsWith(">"))
			{
			  System.out.println("Processing " + s);
				if(segment != null) {
				  //System.out.println("Writing out...");
					writeSixFrameTranslation(annotation, segmentPos, segment.toString(), out);
				}
				String[] token = s.split("\\s+");
				assert(token[0].length() > 1);
				annotation = token[0].substring(1);
				segment = null;
				curPosition = 0;
			}
			else
			{
				for(int i=0; i<s.length(); i++)
				{
					char c = Character.toUpperCase(s.charAt(i));
					if(c == 'A' || c == 'T' || c == 'C' || c == 'G')
					{
						if(segment==null)
						{
							segment = new StringBuffer();
							segmentPos = curPosition;
						}
						segment.append(c);
					}
					else
					{
						if(segment != null) {
						  //System.out.println("Writing out...");
							writeSixFrameTranslation(annotation, segmentPos, segment.toString(), out);
							segment = null;
						}
					}
					curPosition++;
				}
			}
			//curLinePosition = reader.getPosition();
		}
		if(segment != null) {
			writeSixFrameTranslation(annotation, segmentPos, segment.toString(), out);
		}
	}
	
	private void writeSixFrameTranslation(String chromosome, int segmentPos, String segment, PrintStream out)
	{
		if(segment.length() < MIN_PROTEIN_SIZE*3)
			return;
	
//		System.out.println(segment);
//		out = System.out;
		
		String annotation = null;
		// forward
		for(int shift=0; shift<3; shift++)
		{
			annotation = ">"+chromosome+" "+ (segmentPos+shift) + " " + 0 + " " + shift;
			StringBuffer buf = new StringBuffer();

			for(int i=shift; i<segment.length()-2; i+=3)
			{
				char aa = Codon.translate(segment.substring(i, i+3));
				buf.append(aa);
			}
			if(buf.length() > MIN_PROTEIN_SIZE)
			{
				out.println(annotation);
				out.println(buf);
			}
		}
		
		// reverse
		for(int shift=0; shift<3; shift++)
		{
			
			StringBuffer buf = new StringBuffer();
			
			int j=0;
			for(int i=segment.length()-1-shift; i>1; i-=3)
			{
				char aa = Codon.translateRevComplement(segment.substring(i-2, i+1));
				buf.append(aa);
				j=i-2;
			}
			annotation = ">"+chromosome+" "+ (segmentPos+j) + " " + 1 + " " + shift;
			
			if(buf.length() > MIN_PROTEIN_SIZE)
			{
				out.println(annotation);
				out.println(buf);
			}
		}
	}	
	
	private void writeSixFrameTranslationConsideringTermCodon(String chromosome, int segmentPos, String segment, PrintStream out)
	{
		if(segment.length() < MIN_PROTEIN_SIZE*3)
			return;
		
//		out = System.out;
		
		String annotation = null;
		// forward
		for(int shift=0; shift<3; shift++)
		{
			annotation = null;
			StringBuffer buf = null;
			for(int i=shift; i<segment.length()-2; i+=3)
			{
				char aa = Codon.translate(segment.substring(i, i+3));
				if(aa != '*')
				{
					if(buf == null)
					{
						buf = new StringBuffer();
						annotation = ">"+chromosome+" "+ (segmentPos+i) + " " + 0 + " " + shift;
					}
					buf.append(aa);
				}
				else 	// stop codon
				{
					if(buf != null && buf.length() > MIN_PROTEIN_SIZE)
					{
						out.println(annotation);
						out.println(buf);
					}
					buf = null;
				}
			}
			if(buf != null && buf.length() > MIN_PROTEIN_SIZE)
			{
				out.println(annotation);
				out.println(buf);
			}
		}
		
		// reverse
		for(int shift=0; shift<3; shift++)
		{
			annotation = null;
			StringBuffer buf = null;
			for(int i=segment.length()-1-shift; i>1; i-=3)
			{
				char aa = Codon.translateRevComplement(segment.substring(i-2, i+1));
				if(aa != '*')
				{
					if(buf == null)
					{
						buf = new StringBuffer();
						annotation = ">"+chromosome+" "+ (segmentPos+i) + " " + 1 + " " + shift;
					}
					buf.append(aa);
				}
				else 	// stop codon
				{
					if(buf != null && buf.length() > MIN_PROTEIN_SIZE)
					{
						out.println(annotation);
						out.println(buf);
					}
					buf = null;
				}
			}
			if(buf != null && buf.length() > MIN_PROTEIN_SIZE)
			{
				out.println(annotation);
				out.println(buf);
			}
		}
	}
	
	public static String translate(String genome, int shift)
	{
		StringBuffer prot = new StringBuffer();
		for(int i=shift; i<genome.length()-2; i+=3)
		{
			char aa = Codon.translate(genome.substring(i, i+3));
			prot.append(aa);
		}
		return prot.toString();
	}

	public static String translateReverseComplement(String genome, int shift)
	{
		StringBuffer prot = new StringBuffer();
		for(int i=genome.length()-1-shift; i>1; i-=3)
		{
			char aa = Codon.translateRevComplement(genome.substring(i-2, i+1));
			prot.append(aa);
		}
		return prot.toString();
	}
	
	public static void main(String argv[]) 
	{

    //new GenomeTranslator("/home/jung/Data/Databases/Asp/gen/Asp.fasta").translateAndWriteTo("/home/jung/Desktop/test.fasta");
    
		if(argv.length != 2 || !argv[0].contains(".fa") || !argv[1].contains(".fa"))
		{
			System.out.println("usage: java -Xmx(HeapSize) GenomeTranslator genome(*.fasta) translation(*.fasta)");
			System.exit(-1);
		}
		new GenomeTranslator(argv[0]).translateAndWriteTo(argv[1]);
		System.out.println("Done.");
		
		
		/*
		for(int splitNum=0; splitNum<4; splitNum++)
		{
			GenomeTranslator translator = new GenomeTranslator(System.getProperty("user.home")+"/Research/Data/HumanGenome/splitted/Homo_sapiens.NCBI36.54.dna_rm."+splitNum+".fasta");
			translator.translate(System.getProperty("user.home")+"/Research/Data/HumanGenome/translated/HSRM.NCBI36.54.translation."+splitNum+".fasta");
		}
		System.out.println("Done");
		*/
	}
}
