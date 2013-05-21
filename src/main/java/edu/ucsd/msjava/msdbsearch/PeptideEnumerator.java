package edu.ucsd.msjava.msdbsearch;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Modification;
import edu.ucsd.msjava.sequences.Constants;

public class PeptideEnumerator {
	
	private static final int MIN_PEPTIDE_LENGTH = 6;
	private static final int MAX_PEPTIDE_LENGTH = 30;
	private static final int MAX_NUM_MODS = 0;
	private static final int MAX_NUM_MISSED_CLEAVAGES = 2;
	private static final int NTT = 2;
	
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit("Wrong parameter!");
		
		File fastaFile = new File(argv[0]);
		if(!fastaFile.exists())
			printUsageAndExit("File does not exist!");
		if(fastaFile.isDirectory())
			printUsageAndExit("File must not be a directory!");
		if(!fastaFile.getName().endsWith(".fasta") && !fastaFile.getName().endsWith(".fa"))
			printUsageAndExit("Not a fasta file!");
		
		File outputFile = new File(argv[1]);
		enumerate(fastaFile, outputFile);
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println(message);
		System.out.println("Usage: java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.msdbsearch.PeptideEnumerator FastaFile(*.fasta or *.fa) OutputFile");
		System.exit(-1);
	}
	
	public static void enumerate(File fastaFile, File outputFile) throws Exception
	{
		CompactFastaSequence fastaSequence = new CompactFastaSequence(fastaFile.getPath());
		CompactSuffixArray sa = new CompactSuffixArray(fastaSequence, MAX_PEPTIDE_LENGTH);

		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		
		DataInputStream indices = new DataInputStream(new BufferedInputStream(new FileInputStream(sa.getIndexFile())));
		indices.skip(CompactSuffixArray.INT_BYTE_SIZE*2);	// skip size and id
		
		DataInputStream nlcps = new DataInputStream(new BufferedInputStream(new FileInputStream(sa.getNeighboringLcpFile())));
		nlcps.skip(CompactSuffixArray.INT_BYTE_SIZE*2);
		CompactFastaSequence sequence = sa.getSequence();

		int i = Integer.MAX_VALUE - 1000;
		int size = sa.getSize();
		
//		ArrayList<Modification.Instance> mods = new ArrayList<Modification.Instance>();
//		mods.add(new Modification.Instance(Modification.get("Oxidation"), 'M'));
//		mods.add(new Modification.Instance(Modification.get("Carbamidomethyl"), 'C').fixedModification());
//		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(mods);
//		aaSet.setMaxNumberOfVariableModificationsPerPeptide(MAX_NUM_MODS);
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
		
		Enzyme enzyme = Enzyme.TRYPSIN;
		
		CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aaSet, MAX_PEPTIDE_LENGTH);
		int[] numMissedCleavages = new int[MAX_PEPTIDE_LENGTH+1];
		int nnet = 0;
		for(int bufferIndex=0; bufferIndex<size; bufferIndex++)
		{
			int index = indices.readInt();
			int lcp = nlcps.readByte();
			if(lcp >= i+1)		
			{
				continue;
			}
			else if(lcp == 0)	// preceding aa is changed
			{
				char precedingAA = sequence.getCharAt(index);
				if(precedingAA != Constants.TERMINATOR_CHAR && !enzyme.isCleavable(precedingAA))
				{
					i=0;
					nnet = 1;
					if(nnet > 2-NTT)
					{
						continue;
					}
				}
			}			
			if(lcp == 0)
				i = 1;
			else if(lcp < i+1)
				i = lcp;
			
			for(; i<MAX_PEPTIDE_LENGTH+1 && index+i<size-1; i++)	// ith character of a peptide
			{
				char residue = sequence.getCharAt(index+i);
				
				if(candidatePepGrid.addResidue(i, residue) == false)
					break;

				if(enzyme.isCleavable(residue))
					numMissedCleavages[i] = numMissedCleavages[i-1]+1;
				else
					numMissedCleavages[i] = numMissedCleavages[i-1];
				
				if(numMissedCleavages[i] > MAX_NUM_MISSED_CLEAVAGES+1)
					break;
				
				if(i < MIN_PEPTIDE_LENGTH)
					continue;
				
				char next = sequence.getCharAt(index+i+1);
				if(!enzyme.isCleavable(residue) && next != Constants.TERMINATOR_CHAR)
				{
					if(nnet+1 > 2-NTT)
						continue;
				}
				
				for(int j=0; j<candidatePepGrid.size(); j++)
				{
					char pre = sequence.getCharAt(index);
//					String pepSeq = candidatePepGrid.getPeptideSeq(j).replaceAll("m", "M@").replaceAll("C", "C!");
					String pepSeq = candidatePepGrid.getPeptideSeq(j);
					out.println(pre+"."+pepSeq+"."+next);
				}
				if(numMissedCleavages[i] == MAX_NUM_MISSED_CLEAVAGES+1)
					break;
			}
		}
		
		indices.close();
		nlcps.close();
		out.close();
		
		System.out.println("Done");
	}
}
