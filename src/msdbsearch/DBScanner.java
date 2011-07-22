package msdbsearch;

import java.io.FileNotFoundException;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Map.Entry;

import parser.BufferedLineReader;

import msgf.DeNovoGraph;
import msgf.FlexAminoAcidGraph;
import msgf.MSGFDBResultGenerator;
import msgf.GeneratingFunction;
import msgf.GeneratingFunctionGroup;
import msgf.NominalMass;
import msgf.Tolerance;
import msscorer.SimpleDBSearchScorer;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Enzyme;
import msutil.Peptide;
import msutil.SpecKey;
import msutil.Modification.Location;
import sequences.Constants;
import suffixarray.SuffixArraySequence;

public class DBScanner {

	private int minPeptideLength;
	private int maxPeptideLength;
	
	private AminoAcidSet aaSet;
	private double[] aaMass;
	private int[] intAAMass;
	
	private Enzyme enzyme;
	private int numPeptidesPerSpec;

	private final SuffixArrayForMSGFDB sa;
	private final int size;
	// to scan the database partially
	// Input spectra
	private final ScoredSpectraMap specScanner;
	
	// DB search results
	private Map<SpecKey,PriorityQueue<DatabaseMatch>> specKeyDBMatchMap;
	private Map<Integer,PriorityQueue<DatabaseMatch>> specIndexDBMatchMap;

	public DBScanner(
			ScoredSpectraMap specScanner,
			SuffixArrayForMSGFDB sa,
			Enzyme enzyme,
			AminoAcidSet aaSet,
			int numPeptidesPerSpec,
			int minPeptideLength,
			int maxPeptideLength
			) 
	{
		this.specScanner = specScanner;
		this.sa = sa;
		this.size = sa.getSize();
		this.aaSet = aaSet;
		this.enzyme = enzyme;
		this.numPeptidesPerSpec = numPeptidesPerSpec;
		this.minPeptideLength = minPeptideLength;
		this.maxPeptideLength = maxPeptideLength;
		
		// Initialize mass arrays for a faster search
		aaMass = new double[aaSet.getMaxResidue()];
		intAAMass = new int[aaSet.getMaxResidue()];
		for(int i=0; i<aaMass.length; i++)
		{
			aaMass[i] = -1;
			intAAMass[i] = -1;
		}
		for(AminoAcid aa : aaSet.getAllAminoAcidArr())
		{
			aaMass[aa.getResidue()] = aa.getAccurateMass();
			intAAMass[aa.getResidue()] = aa.getNominalMass();
		}
		
		specKeyDBMatchMap = Collections.synchronizedMap(new HashMap<SpecKey,PriorityQueue<DatabaseMatch>>());
		specIndexDBMatchMap = Collections.synchronizedMap(new HashMap<Integer,PriorityQueue<DatabaseMatch>>());
	}

	// builder
	public DBScanner maxPeptideLength(int maxPeptideLength)
	{
		this.maxPeptideLength = maxPeptideLength;
		return this;
	}

	// builder
	public DBScanner minPeptideLength(int minPeptideLength)
	{
		if(minPeptideLength > 1)
			this.minPeptideLength = minPeptideLength;
		else
			minPeptideLength = 1;
		return this;
	}
	
	public synchronized void addDBMatches(Map<SpecKey,PriorityQueue<DatabaseMatch>> map)
	{
		if(map == null)
			return;
		Iterator<Entry<SpecKey, PriorityQueue<DatabaseMatch>>> itr = map.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<SpecKey, PriorityQueue<DatabaseMatch>> entry = itr.next();
			SpecKey specKey = entry.getKey(); 
			PriorityQueue<DatabaseMatch> queue = specKeyDBMatchMap.get(entry.getKey());
			if(queue == null)
			{
				queue = new PriorityQueue<DatabaseMatch>();
				specKeyDBMatchMap.put(specKey, queue);
			}
			for(DatabaseMatch match : entry.getValue())
			{
				if(queue.size() < this.numPeptidesPerSpec)
				{
					queue.add(match);
				}
				else if(queue.size() >= this.numPeptidesPerSpec)
				{
					if(match.getScore() > queue.peek().getScore())
					{
						queue.poll();
						queue.add(match);
					}
				}
			}
		}
	}
	
	
	public void dbSearchCTermEnzymeNoMod(int numberOfAllowableNonEnzymaticTermini, boolean verbose)
	{
		dbSearchCTermEnzymeNoMod(numberOfAllowableNonEnzymaticTermini, 0, size, verbose);
	}

	public void dbSearchCTermEnzyme(int numberOfAllowableNonEnzymaticTermini, boolean verbose)
	{
		dbSearchCTermEnzyme(numberOfAllowableNonEnzymaticTermini, 0, size, verbose);
	}
	
	public void dbSearchNTermEnzyme(int numberOfAllowableNonEnzymaticTermini, boolean verbose)
	{
		dbSearchNTermEnzyme(numberOfAllowableNonEnzymaticTermini, 0, size, verbose);
	}
	
	public void dbSearchNoEnzyme(boolean verbose)
	{
		dbSearchNoEnzyme(0, size, verbose);
	}
	
	
	// duplicated for speeding-up the search
	public void dbSearchCTermEnzymeNoMod(
			int numberOfAllowableNonEnzymaticTermini, 
			int fromIndex, 
			int toIndex, 
			boolean verbose)
	{
		IntBuffer indexBuf = sa.getIndices().duplicate();
		ByteBuffer lcpBuf = sa.getNeighboringLcps().duplicate();
		SuffixArraySequence sequence = sa.getSequence();
		
		Map<SpecKey,PriorityQueue<DatabaseMatch>> curSpecKeyDBMatchMap = new HashMap<SpecKey,PriorityQueue<DatabaseMatch>>();
		double[] prm = new double[maxPeptideLength+2];
		prm[0] = 0;
		int[] intPRM = new int[maxPeptideLength+2];
		intPRM[0] = 0;
		
		int i = Integer.MAX_VALUE;
		
		boolean enzymaticSearch;
		if(numberOfAllowableNonEnzymaticTermini == 2)
			enzymaticSearch = false;
		else
			enzymaticSearch = true;
		
		int neighboringAACleavageCredit = aaSet.getNeighboringAACleavageCredit();
		int neighboringAACleavagePenalty = aaSet.getNeighboringAACleavagePenalty();
		int peptideCleavageCredit = aaSet.getPeptideCleavageCredit();
		int peptideCleavagePenalty = aaSet.getPeptideCleavagePenalty();
		
		Tolerance leftPMTolerance = specScanner.getLeftParentMassTolerance();
		Tolerance rightPMTolerance = specScanner.getRightParentMassTolerance();
		
		boolean isProteinNTerm = true;
		int nTermAAScore = neighboringAACleavageCredit;
		int nNET = 0;	// number of non-enzymatic termini
		for(int bufferIndex=fromIndex; bufferIndex<toIndex; bufferIndex++)
		{
			int index = indexBuf.get(bufferIndex);
			int lcp;
			if(bufferIndex == fromIndex)
				lcp = 0;
			else
				lcp = lcpBuf.get(bufferIndex);
			if(lcp > i)		// skip redundant peptide
				continue;
			else if(lcp == 0)	// preceding aa is changed
			{
				char nAA = sequence.getCharAt(index);
				if(nAA == Constants.TERMINATOR_CHAR)
					isProteinNTerm = true;
				else
				{	
					if(nAA == 'M' && sequence.getCharAt(index-1) == Constants.TERMINATOR_CHAR)
						isProteinNTerm = true;
					else
						isProteinNTerm = false;
				}
				if(isProteinNTerm || enzyme.isCleavable(nAA))
					nTermAAScore = neighboringAACleavageCredit;
				else
					nTermAAScore = neighboringAACleavagePenalty;
				if(enzymaticSearch)
				{
					if(!isProteinNTerm && !enzyme.isCleavable(nAA))
					{
						nNET = 1;
						if(nNET > numberOfAllowableNonEnzymaticTermini)
						{
							i=0;
							continue;
						}
					}
					else
						nNET = 0;
				}
			}
			
//			if(verbose && rank % 1000000 == 0)
//				System.out.println("DBSearch: " + rank/(float)size*100 + "%");
			
			for(i=lcp > 1 ? lcp : 1; i<maxPeptideLength+1 && index+i<size; i++)	// ith character of a peptide
			{
				char residue = sequence.getCharAt(index+i);
				double m = aaMass[residue];
				if(m <= 0)
					break;
				prm[i+1] = prm[i] + m;	// prm[i]: mass of peptide upto ith residue
				intPRM[i+1] = intPRM[i] + intAAMass[residue];
				boolean isProteinCTerm = index+i+1 == size || sequence.getCharAt(index+i+1) == Constants.TERMINATOR_CHAR;
				
				if(enzymaticSearch && !enzyme.isCleavable(residue))
				{
					if(!isProteinCTerm && nNET+1 > numberOfAllowableNonEnzymaticTermini)
						continue;
				}
				if(i < minPeptideLength)
					continue;
				float peptideMass = (float)prm[i+1];
				int nominalPeptideMass = intPRM[i+1];
				float tolDaLeft = leftPMTolerance.getToleranceAsDa(peptideMass);
				float tolDaRight = rightPMTolerance.getToleranceAsDa(peptideMass);
				
				double leftThr = (double)(peptideMass - tolDaRight);
				double rightThr = (double)(peptideMass + tolDaLeft);
				Collection<SpecKey> matchedSpecKeyList = specScanner.getPepMassSpecKeyMap().subMap(leftThr, rightThr).values();
				if(matchedSpecKeyList.size() > 0)
				{
//					String pepStr = sequence.getSubsequence(index+1, index+i+1);
//					if(pepStr.contains("X"))
//						System.out.println("*****************"+index);
//					if(index == 1)
//						System.out.println("Debug");
					int peptideCleavageScore;
					if(enzyme.isCleavable(sequence.getCharAt(index+i)))
						peptideCleavageScore = peptideCleavageCredit;
					else
						peptideCleavageScore = peptideCleavagePenalty;
					
					for(SpecKey specKey : matchedSpecKeyList)
					{
						SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(specKey);
						int score = nTermAAScore + scorer.getScore(prm, intPRM, 2, i+2) + peptideCleavageScore;
						PriorityQueue<DatabaseMatch> prevMatchQueue = curSpecKeyDBMatchMap.get(specKey);
						if(prevMatchQueue == null)
						{
							prevMatchQueue = new PriorityQueue<DatabaseMatch>();
							curSpecKeyDBMatchMap.put(specKey, prevMatchQueue);
						}
						if(prevMatchQueue.size() < this.numPeptidesPerSpec)
						{
							prevMatchQueue.add(new DatabaseMatch(index, i+2, score, nominalPeptideMass));
						}
						else if(prevMatchQueue.size() >= this.numPeptidesPerSpec)
						{
							if(score > prevMatchQueue.peek().getScore())
							{
								prevMatchQueue.poll();
								prevMatchQueue.add(new DatabaseMatch(index, i+2, score, nominalPeptideMass));
							}
						}
					}
				}
			}
		}
		this.addDBMatches(curSpecKeyDBMatchMap);
	}  	

	public void dbSearchCTermEnzyme(int numberOfAllowableNonEnzymaticTermini, int fromIndex, int toIndex, boolean verbose)
	{
		Map<SpecKey,PriorityQueue<DatabaseMatch>> curSpecKeyDBMatchMap = new HashMap<SpecKey,PriorityQueue<DatabaseMatch>>();
		
		CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aaSet, maxPeptideLength);
		
		int i = Integer.MAX_VALUE;
		
		boolean enzymaticSearch;
		if(numberOfAllowableNonEnzymaticTermini == 2)
			enzymaticSearch = false;
		else
			enzymaticSearch = true;
		
		int neighboringAACleavageCredit = aaSet.getNeighboringAACleavageCredit();
		int neighboringAACleavagePenalty = aaSet.getNeighboringAACleavagePenalty();
		int peptideCleavageCredit = aaSet.getPeptideCleavageCredit();
		int peptideCleavagePenalty = aaSet.getPeptideCleavagePenalty();
		
		boolean containsCTermMod = aaSet.containsCTermModification();
		
		IntBuffer indexBuf = sa.getIndices().duplicate();
		ByteBuffer lcpBuf = sa.getNeighboringLcps().duplicate();
		SuffixArraySequence sequence = sa.getSequence();
		
		boolean isProteinNTerm = true;
		int nTermAAScore = neighboringAACleavageCredit;
		boolean isExtensionAtTheSameIndex;
		int nNET = 0;	// number of non-enzymatic termini
		for(int bufferIndex=fromIndex; bufferIndex<toIndex; bufferIndex++) 
		{
			isExtensionAtTheSameIndex = false;
			int index = indexBuf.get(bufferIndex);
			int lcp;
			if(bufferIndex == fromIndex)
				lcp = 0;
			else
				lcp = lcpBuf.get(bufferIndex);
			
			if(lcp > i)		// skip redundant peptide
				continue;
			else if(lcp == 0)	// preceding aa is changed
			{
				char nAA = sequence.getCharAt(index);
				if(nAA == Constants.TERMINATOR_CHAR)
					isProteinNTerm = true;
				else
				{
					if(nAA == 'M' && sequence.getCharAt(index-1) == Constants.TERMINATOR_CHAR)
						isProteinNTerm = true;
					else
						isProteinNTerm = false;
				}
				if(isProteinNTerm || enzyme.isCleavable(nAA))
					nTermAAScore = neighboringAACleavageCredit;
				else
					nTermAAScore = neighboringAACleavagePenalty;
				if(enzymaticSearch)
				{
					if(!isProteinNTerm && !enzyme.isCleavable(nAA))
					{
						nNET = 1;
						if(nNET > numberOfAllowableNonEnzymaticTermini)
						{
							i=0;
							continue;
						}
					}
					else
						nNET = 0;
				}
			}
			
			
//			if(verbose && rank % 1000000 == 0)
//				System.out.println("DBSearch: " + rank/(float)size*100 + "%");

			for(i=lcp > 1 ? lcp : 1; i<maxPeptideLength+1 && index+i<size-1; i++)	// ith character of a peptide
			{
				char residue = sequence.getCharAt(index+i);
				boolean isProteinCTerm = false;
				if(i==1)
				{
					if(isProteinNTerm)
					{
						if(candidatePepGrid.addProtNTermResidue(residue) == false)
							break;
					}
					else
					{
						if(candidatePepGrid.addNTermResidue(residue) == false)
							break;
					}
				}
				else
				{
					if(!containsCTermMod)
					{
						if(candidatePepGrid.addResidue(i, residue) == false)
							break;
					}
					else
					{
						if(i < minPeptideLength)
						{
							if(candidatePepGrid.addResidue(i, residue) == false)
								break;
							else
								continue;
						}
						else
						{
							if(isExtensionAtTheSameIndex && i > minPeptideLength)
								candidatePepGrid.addResidue(i-1, sequence.getCharAt(index+i-1));
							boolean success;
							if(isProteinCTerm = (sequence.getCharAt(index+i+1) == Constants.TERMINATOR_CHAR))	// protein C-term
								success = candidatePepGrid.addProtCTermResidue(i, residue);
							else	// peptide C-term
								success = candidatePepGrid.addCTermResidue(i, residue);
							if(!success)
								break;
						}
					}
				}
				
//				if(sequence.getSubsequence(index, index+i+1).equalsIgnoreCase("KTQDAHFQR"))
//					System.out.println("Debug");

				if(i < minPeptideLength)
					continue;
				
				if(enzymaticSearch && !enzyme.isCleavable(residue))
				{
					if(!isProteinCTerm && nNET+1 > numberOfAllowableNonEnzymaticTermini)
							continue;
				}
				
				int peptideCleavageScore;
				if(enzyme.isCleavable(sequence.getCharAt(index+i)))
					peptideCleavageScore = peptideCleavageCredit;
				else
					peptideCleavageScore = peptideCleavagePenalty;
				
				int enzymeScore = nTermAAScore + peptideCleavageScore;
				
				for(int j=0; j<candidatePepGrid.size(); j++)
				{
					float peptideMass = candidatePepGrid.getPeptideMass(j);
					int nominalPeptideMass = candidatePepGrid.getNominalPeptideMass(j);
					float tolDaLeft = specScanner.getLeftParentMassTolerance().getToleranceAsDa(peptideMass);
					float tolDaRight = specScanner.getRightParentMassTolerance().getToleranceAsDa(peptideMass);
					
					double leftThr = (double)(peptideMass - tolDaRight);
					double rightThr = (double)(peptideMass + tolDaLeft);
					Collection<SpecKey> matchedSpecKeyList = specScanner.getPepMassSpecKeyMap().subMap(leftThr, rightThr).values();
					if(matchedSpecKeyList.size() > 0)
					{
						for(SpecKey specKey : matchedSpecKeyList)
						{
							SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(specKey);
							int score = enzymeScore + scorer.getScore(candidatePepGrid.getPRMGrid()[j], candidatePepGrid.getNominalPRMGrid()[j], 1, i+1); 
							PriorityQueue<DatabaseMatch> prevMatchQueue = curSpecKeyDBMatchMap.get(specKey);
							if(prevMatchQueue == null)
							{
								prevMatchQueue = new PriorityQueue<DatabaseMatch>();
								curSpecKeyDBMatchMap.put(specKey, prevMatchQueue);
							}
							if(prevMatchQueue.size() < this.numPeptidesPerSpec)
							{
								prevMatchQueue.add(new DatabaseMatch(index, i+2, score, nominalPeptideMass).pepSeq(candidatePepGrid.getPeptideSeq(j)).setProteinNTerm(isProteinNTerm).setProteinCTerm(isProteinCTerm));
							}
							else if(prevMatchQueue.size() >= this.numPeptidesPerSpec)
							{
								if(score > prevMatchQueue.peek().getScore())
								{
									prevMatchQueue.poll();
									prevMatchQueue.add(new DatabaseMatch(index, i+2, score, nominalPeptideMass).pepSeq(candidatePepGrid.getPeptideSeq(j)).setProteinNTerm(isProteinNTerm).setProteinCTerm(isProteinCTerm));
								}
							}
						}
					}					
				}
				isExtensionAtTheSameIndex = true;
			}
		}
//		indices.rewind();
//		neighboringLcps.rewind();
		this.addDBMatches(curSpecKeyDBMatchMap);
	}  		
	
	// dupulicated to speed-up the search
	public void dbSearchNoEnzyme(int fromIndex, int toIndex, boolean verbose)
	{
		Map<SpecKey,PriorityQueue<DatabaseMatch>> curSpecKeyDBMatchMap = new HashMap<SpecKey,PriorityQueue<DatabaseMatch>>();
		CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aaSet, maxPeptideLength);
		
		int i = 0;
		
		IntBuffer indexBuf = sa.getIndices().duplicate();
		ByteBuffer lcpBuf = sa.getNeighboringLcps().duplicate();
		SuffixArraySequence sequence = sa.getSequence();
		
		boolean containsCTermMod = aaSet.containsCTermModification();
		
		boolean isExtensionAtTheSameIndex;
		boolean isProteinNTerm = true;
		for(int bufferIndex=fromIndex; bufferIndex<toIndex; bufferIndex++)
		{
			int index = indexBuf.get(bufferIndex);
			int lcp;
			if(bufferIndex == fromIndex)
				lcp = 0;
			else
				lcp = lcpBuf.get(bufferIndex);

			if(lcp > i)		// skip redundant peptide
				continue;
			else if(lcp == 0)
			{
				if(index != 0)
				{
					char nAA = sequence.getCharAt(index-1);
					if(nAA == Constants.TERMINATOR_CHAR)
						isProteinNTerm = true;
					else
					{
						if(nAA == 'M' && sequence.getCharAt(index-2) == Constants.TERMINATOR_CHAR)
							isProteinNTerm = true;
						else
							isProteinNTerm = false;
					}
				}
				else
					continue;
			}
			
			isExtensionAtTheSameIndex = false;
//			if(verbose && rank % 1000000 == 0)
//				System.out.println("DBSearch: " + rank/(float)size*100 + "%");

			for(i=lcp; i<maxPeptideLength && index+i<size-1; i++)	// (i+1)th character of a peptide, length=i+1
			{
				char residue = sequence.getCharAt(index+i);
				if(i==0)
				{
					if(isProteinNTerm)
					{
						if(candidatePepGrid.addProtNTermResidue(residue) == false)
							break;
					}
					else
					{
						if(candidatePepGrid.addNTermResidue(residue) == false)
							break;
					}
				}
				else
				{
					int length = i+1;
					if(!containsCTermMod)
					{
						if(candidatePepGrid.addResidue(length, residue) == false)
							break;
					}
					else
					{
						if(length < minPeptideLength)
						{
							if(candidatePepGrid.addResidue(length, residue) == false)
								break;
							else
								continue;
						}
						else
						{
							if(isExtensionAtTheSameIndex && length > minPeptideLength)
								candidatePepGrid.addResidue(length-1, sequence.getCharAt(index+i-1));
							boolean success;
							if(sequence.getCharAt(index+i+1) == Constants.TERMINATOR_CHAR)	// protein C-term
								success = candidatePepGrid.addProtCTermResidue(length, residue);
							else	// peptide C-term
								success = candidatePepGrid.addCTermResidue(length, residue);
							if(!success)
								break;
						}
					}
				}
				
				if(i+1 < minPeptideLength)
					continue;
//				if(sequence.getSubsequence(index, index+i+1).equalsIgnoreCase("TPEVTCVVVDVSHEDPEVQFK"))
//					System.out.println("Debug");
				char nextAA = sequence.getCharAt(index+i+1);
				boolean isProteinCTerm = nextAA == Constants.TERMINATOR_CHAR;
				
				for(int j=0; j<candidatePepGrid.size(); j++)
				{
					float peptideMass = candidatePepGrid.getPeptideMass(j);
					int nominalPeptideMass = candidatePepGrid.getNominalPeptideMass(j);
					float tolDaLeft = specScanner.getLeftParentMassTolerance().getToleranceAsDa(peptideMass);
					float tolDaRight = specScanner.getRightParentMassTolerance().getToleranceAsDa(peptideMass);
					
					double leftThr = (double)(peptideMass - tolDaRight);
					double rightThr = (double)(peptideMass + tolDaLeft);
					Collection<SpecKey> matchedSpecKeyList = specScanner.getPepMassSpecKeyMap().subMap(leftThr, rightThr).values();
					if(matchedSpecKeyList.size() > 0)
					{
						for(SpecKey specKey : matchedSpecKeyList)
						{
							SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(specKey);
							int score =  scorer.getScore(candidatePepGrid.getPRMGrid()[j], candidatePepGrid.getNominalPRMGrid()[j], 1, i+2); 
							PriorityQueue<DatabaseMatch> prevMatchQueue = curSpecKeyDBMatchMap.get(specKey);
							if(prevMatchQueue == null)
							{
								prevMatchQueue = new PriorityQueue<DatabaseMatch>();
								curSpecKeyDBMatchMap.put(specKey, prevMatchQueue);
							}
							if(prevMatchQueue.size() < this.numPeptidesPerSpec)
							{
								prevMatchQueue.add(new DatabaseMatch(index-1, i+2, score, nominalPeptideMass).pepSeq(candidatePepGrid.getPeptideSeq(j)).setProteinNTerm(isProteinNTerm).setProteinCTerm(isProteinCTerm));
							}
							else if(prevMatchQueue.size() >= this.numPeptidesPerSpec)
							{
								if(score > prevMatchQueue.peek().getScore())
								{
									prevMatchQueue.poll();
									prevMatchQueue.add(new DatabaseMatch(index-1, i+2, score, nominalPeptideMass).pepSeq(candidatePepGrid.getPeptideSeq(j)).setProteinNTerm(isProteinNTerm).setProteinCTerm(isProteinCTerm));
								}
							}
						}
					}					
				}
				isExtensionAtTheSameIndex = true;
			}
		}
//		indices.rewind();
//		neighboringLcps.rewind();
		this.addDBMatches(curSpecKeyDBMatchMap);
	}
	
	public void dbSearchNTermEnzyme(int numberOfAllowableNonEnzymaticTermini, int fromIndex, int toIndex, boolean verbose)
	{
		Map<SpecKey,PriorityQueue<DatabaseMatch>> curSpecKeyDBMatchMap = new HashMap<SpecKey,PriorityQueue<DatabaseMatch>>();
		CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aaSet, maxPeptideLength);
		
		int i = Integer.MAX_VALUE;
		
		IntBuffer indexBuf = sa.getIndices().duplicate();
		ByteBuffer lcpBuf = sa.getNeighboringLcps().duplicate();
		SuffixArraySequence sequence = sa.getSequence();

		boolean enzymaticSearch;
		if(enzyme == null || numberOfAllowableNonEnzymaticTermini == 2)
			enzymaticSearch = false;	// consider all positions in the db
		else
			enzymaticSearch = true;
		
		int neighboringAACleavageCredit = aaSet.getNeighboringAACleavageCredit();
		int neighboringAACleavagePenalty = aaSet.getNeighboringAACleavagePenalty();
		int peptideCleavageCredit = aaSet.getPeptideCleavageCredit();
		int peptideCleavagePenalty = aaSet.getPeptideCleavagePenalty();
		
		boolean containsCTermMod = aaSet.containsCTermModification();
		
		boolean isProteinNTerm = true;
		boolean isExtensionAtTheSameIndex;
		int peptideCleavageScore = 0;	// N-term residue, when i=0
		
		int nNET = 0;	// number of non-enzymatic termini
		for(int bufferIndex=fromIndex; bufferIndex<toIndex; bufferIndex++)
		{
			int index = indexBuf.get(bufferIndex);
			int lcp;
			if(bufferIndex == fromIndex)
				lcp = 0;
			else
				lcp = lcpBuf.get(bufferIndex);
			isExtensionAtTheSameIndex = false;

			if(lcp > i)		// skip redundant peptide
				continue;
			else if(lcp == 0)	
			{
				i=0;
				char nAA = sequence.getCharAt(index);
				if(enzyme.isCleavable(nAA))
				{
					peptideCleavageScore = peptideCleavageCredit;
					nNET = 0;
				}
				else
				{
					if(!Character.isLetter(nAA))
						continue;
					if(enzymaticSearch)
					{
						nNET = 1;
						if(nNET > numberOfAllowableNonEnzymaticTermini)
							continue;
					}
					peptideCleavageScore = peptideCleavagePenalty;
					
				}
				char preAA = sequence.getCharAt(index-1);
				if(preAA == Constants.TERMINATOR_CHAR)
					isProteinNTerm = true;
				else
				{
					if(preAA == 'M' && sequence.getCharAt(index-2) == Constants.TERMINATOR_CHAR)
						isProteinNTerm = true;
					else
						isProteinNTerm = false;
				}
			}
			
//			if(verbose && rank % 1000000 == 0)
//				System.out.println("DBSearch: " + rank/(float)size*100 + "%");

			for(i=lcp; i<maxPeptideLength+1 && index+i<size-1; i++)	// (i+1)th character of a peptide, length=i+1
			{
				char residue = sequence.getCharAt(index+i);
				if(i==0)
				{
					if(isProteinNTerm)
					{
						if(candidatePepGrid.addProtNTermResidue(residue) == false)
							break;
					}
					else
					{
						if(candidatePepGrid.addNTermResidue(residue) == false)
							break;
					}
				}
				else
				{
					int length = i+1;
					if(!containsCTermMod)
					{
						if(candidatePepGrid.addResidue(length, residue) == false)
							break;
					}
					else
					{
						if(length < minPeptideLength)
						{
							if(candidatePepGrid.addResidue(length, residue) == false)
								break;
							else
								continue;
						}
						else
						{
							if(isExtensionAtTheSameIndex && length > minPeptideLength)
								candidatePepGrid.addResidue(length-1, sequence.getCharAt(index+i-1));
							boolean success;
							if(sequence.getCharAt(index+i+1) == Constants.TERMINATOR_CHAR)	// protein C-term
								success = candidatePepGrid.addProtCTermResidue(length, residue);
							else	// peptide C-term
								success = candidatePepGrid.addCTermResidue(length, residue);
							if(!success)
								break;
						}
					}
				}
				
				if(i+1 < minPeptideLength)
					continue;
				
				int cTermAAScore = 0;
				char nextAA = sequence.getCharAt(index+i+1);
				boolean isProteinCTerm = nextAA == Constants.TERMINATOR_CHAR;
				if(isProteinCTerm || enzyme.isCleavable(nextAA))
				{
					cTermAAScore = neighboringAACleavageCredit;
				}
				else
				{
					if(nNET+1 > numberOfAllowableNonEnzymaticTermini)
						continue;
					cTermAAScore = neighboringAACleavagePenalty;
				}
				
				int enzymeScore = cTermAAScore + peptideCleavageScore;
				
				for(int j=0; j<candidatePepGrid.size(); j++)
				{
					float peptideMass = candidatePepGrid.getPeptideMass(j);
					int nominalPeptideMass = candidatePepGrid.getNominalPeptideMass(j);
					float tolDaLeft = specScanner.getLeftParentMassTolerance().getToleranceAsDa(peptideMass);
					float tolDaRight = specScanner.getRightParentMassTolerance().getToleranceAsDa(peptideMass);
					
					double leftThr = (double)(peptideMass - tolDaRight);
					double rightThr = (double)(peptideMass + tolDaLeft);
					Collection<SpecKey> matchedSpecKeyList = specScanner.getPepMassSpecKeyMap().subMap(leftThr, rightThr).values();
					if(matchedSpecKeyList.size() > 0)
					{
						for(SpecKey specKey : matchedSpecKeyList)
						{
							SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(specKey);
							int score = enzymeScore + scorer.getScore(candidatePepGrid.getPRMGrid()[j], candidatePepGrid.getNominalPRMGrid()[j], 1, i+2); 
							PriorityQueue<DatabaseMatch> prevMatchQueue = curSpecKeyDBMatchMap.get(specKey);
							if(prevMatchQueue == null)
							{
								prevMatchQueue = new PriorityQueue<DatabaseMatch>();
								curSpecKeyDBMatchMap.put(specKey, prevMatchQueue);
							}
							if(prevMatchQueue.size() < this.numPeptidesPerSpec)
							{
								prevMatchQueue.add(new DatabaseMatch(index-1, i+2, score, nominalPeptideMass).pepSeq(candidatePepGrid.getPeptideSeq(j)).setProteinNTerm(isProteinNTerm).setProteinCTerm(isProteinCTerm));
							}
							else if(prevMatchQueue.size() >= this.numPeptidesPerSpec)
							{
								if(score > prevMatchQueue.peek().getScore())
								{
									prevMatchQueue.poll();
									prevMatchQueue.add(new DatabaseMatch(index-1, i+2, score, nominalPeptideMass).pepSeq(candidatePepGrid.getPeptideSeq(j)).setProteinNTerm(isProteinNTerm).setProteinCTerm(isProteinCTerm));
								}
							}
						}
					}					
				}
				isExtensionAtTheSameIndex = true;
			}
		}
//		indices.rewind();
//		neighboringLcps.rewind();
		this.addDBMatches(curSpecKeyDBMatchMap);
	}
	
	public void computeSpecProb(boolean storeScoreDist)
	{
		computeSpecProb(storeScoreDist, 0, specScanner.getSpecKeyList().size());
	}
	
	public void computeSpecProb(boolean storeScoreDist, int fromIndex, int toIndex)
	{
		List<SpecKey> specKeyList = specScanner.getSpecKeyList().subList(fromIndex, toIndex);
		
		for(SpecKey specKey : specKeyList)
		{
			PriorityQueue<DatabaseMatch> matchQueue = specKeyDBMatchMap.get(specKey);
			if(matchQueue == null)
				continue;

			int specIndex = specKey.getSpecIndex();
			
			boolean useProtNTerm = false;
			boolean useProtCTerm = false;
			int minScore = Integer.MAX_VALUE;
//			int minMatchedNominalMass = Integer.MAX_VALUE;
//			int maxMatchedNominalMass = Integer.MIN_VALUE;
			for(DatabaseMatch m : matchQueue)
			{
				if(m.isProteinNTerm())
					useProtNTerm = true;
				if(m.isProteinCTerm())
					useProtCTerm = true;
				if(m.getScore() < minScore)
					minScore = m.getScore();
//				if(m.getNominalPeptideMass() < minMatchedNominalMass)
//					minMatchedNominalMass = m.getNominalPeptideMass();
//				if(m.getNominalPeptideMass() > maxMatchedNominalMass)
//					maxMatchedNominalMass = m.getNominalPeptideMass();
			}
			
			GeneratingFunctionGroup<NominalMass> gf = new GeneratingFunctionGroup<NominalMass>();
			SimpleDBSearchScorer<NominalMass> scoredSpec = specScanner.getSpecKeyScorerMap().get(specKey);
//			float peptideMass = spec.getParentMass() - (float)Composition.H2O;
			float peptideMass = scoredSpec.getPrecursorPeak().getMass() - (float)Composition.H2O;
			int nominalPeptideMass = NominalMass.toNominalMass(peptideMass);
			float tolDaLeft = specScanner.getLeftParentMassTolerance().getToleranceAsDa(peptideMass);
			float tolDaRight = specScanner.getRightParentMassTolerance().getToleranceAsDa(peptideMass);
			int maxPeptideMassIndex, minPeptideMassIndex;
			maxPeptideMassIndex = nominalPeptideMass + Math.round(tolDaLeft-0.4999f);
			minPeptideMassIndex = nominalPeptideMass - Math.round(tolDaRight-0.4999f);
			
			if(tolDaRight < 0.5f)
				minPeptideMassIndex -= specScanner.getNumAllowedC13();

//			if(maxPeptideMassIndex < maxMatchedNominalMass)
//				maxPeptideMassIndex = maxMatchedNominalMass;
//			if(minPeptideMassIndex > minMatchedNominalMass)
//				minPeptideMassIndex = minMatchedNominalMass;
			
			for(int peptideMassIndex = minPeptideMassIndex; peptideMassIndex<=maxPeptideMassIndex; peptideMassIndex++)
			{
				DeNovoGraph<NominalMass> graph = new FlexAminoAcidGraph(
						aaSet, 
						peptideMassIndex,
						enzyme,
						scoredSpec,
						useProtNTerm,
						useProtCTerm
						);
				
				GeneratingFunction<NominalMass> gfi = new GeneratingFunction<NominalMass>(graph)
				.doNotBacktrack()
				.doNotCalcNumber();
				gfi.setUpScoreThreshold(minScore);
				gf.registerGF(graph.getPMNode(), gfi);
			}

			boolean isGFComputed = gf.computeGeneratingFunction();
			
			for(DatabaseMatch match : matchQueue)
			{
				if(!isGFComputed || match.getNominalPeptideMass() < minPeptideMassIndex || match.getNominalPeptideMass() > maxPeptideMassIndex)
				{
					match.setDeNovoScore(Integer.MIN_VALUE);
					match.setSpecProb(1);
				}
				else
				{
					match.setDeNovoScore(gf.getMaxScore()-1);
					int score = match.getScore();
					double specProb = gf.getSpectralProbability(score);
					assert(specProb > 0): specIndex + ": " + match.getDeNovoScore()+" "+match.getScore()+" "+specProb; 
					match.setSpecProb(specProb);
					if(storeScoreDist)
						match.setScoreDist(gf.getScoreDist());
				}
			}
		}
	}
	
	public synchronized void addDBSearchResults(MSGFDBResultGenerator gen, String specFileName)
	{
		Iterator<Entry<SpecKey, PriorityQueue<DatabaseMatch>>> itr = specKeyDBMatchMap.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<SpecKey, PriorityQueue<DatabaseMatch>> entry = itr.next();
			SpecKey specKey = entry.getKey();
			PriorityQueue<DatabaseMatch> matchQueue = entry.getValue();
			if(matchQueue == null || matchQueue.size() == 0)
				continue;
			
			int specIndex = specKey.getSpecIndex();
			int charge = specKey.getCharge();
			PriorityQueue<DatabaseMatch> existingQueue = specIndexDBMatchMap.get(specIndex);
			if(existingQueue == null)
			{
				existingQueue = new PriorityQueue<DatabaseMatch>(this.numPeptidesPerSpec, new DatabaseMatch.SpecProbComparator());
				specIndexDBMatchMap.put(specIndex, existingQueue);
			}
			
			for(DatabaseMatch match : matchQueue)
			{
				match.setCharge(charge);
				if(existingQueue.size() < this.numPeptidesPerSpec)
				{
					existingQueue.add(match);
				}
				else if(existingQueue.size() >= this.numPeptidesPerSpec)
				{
					if(match.getSpecProb() < existingQueue.peek().getSpecProb())
					{
						existingQueue.poll();
						existingQueue.add(match);
					}
				}
			}
		}		
		
		Iterator<Entry<Integer, PriorityQueue<DatabaseMatch>>> itr2 = specIndexDBMatchMap.entrySet().iterator();
		while(itr2.hasNext())
		{
			Entry<Integer, PriorityQueue<DatabaseMatch>> entry = itr2.next();
			int specIndex = entry.getKey();
			PriorityQueue<DatabaseMatch> matchQueue = entry.getValue();
			if(matchQueue == null)
				continue;

			ArrayList<DatabaseMatch> matchList = new ArrayList<DatabaseMatch>(matchQueue);
			if(matchList.size() == 0)
				continue;

			for(int i=matchList.size()-1; i>=0; --i)
			{
				DatabaseMatch match = matchList.get(i);
				
				if(match.getDeNovoScore() < 0)
					continue;
				
				int index = match.getIndex();
				int length = match.getLength();
				
				String peptideStr = match.getPepSeq();
				if(peptideStr == null)
					peptideStr = sa.getSequence().getSubsequence(index+1, index+length-1);
//				if(peptideStr.contains("X"))
//					System.out.println("Debug");
				Peptide pep = aaSet.getPeptide(peptideStr);
				String annotationStr = sa.getSequence().getCharAt(index)+"."+pep+"."+sa.getSequence().getCharAt(index+length-1);
				SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(new SpecKey(specIndex, match.getCharge()));
//				float expMass = scorer.getParentMass();
				float expMass = scorer.getPrecursorPeak().getMass();
				float theoMass = pep.getParentMass();
				float pmError = Float.MAX_VALUE;
				float peptideMass = expMass - (float)Composition.H2O;
				
				float tolDaRight = specScanner.getRightParentMassTolerance().getToleranceAsDa(peptideMass);
				int nC13 = tolDaRight >= 0.5f ? 0 : specScanner.getNumAllowedC13();
				for(int numC13=0; numC13<=nC13; numC13++)
				{
					float error = expMass-theoMass-(float)(Composition.ISOTOPE)*numC13; 
					if(Math.abs(error) < Math.abs(pmError))
						pmError = error;
				}
				if(specScanner.getRightParentMassTolerance().isTolerancePPM())
					pmError = pmError/theoMass*1e6f;
				
				String protein = sa.getAnnotation(index+1);
				
				int score = match.getScore();
				double specProb = match.getSpecProb();
				int numPeptides = sa.getNumDistinctPeptides(peptideStr.length()+1);
				double pValue = MSGFDBResultGenerator.DBMatch.getPValue(specProb, numPeptides);
				String specProbStr;
				if(specProb < Float.MIN_NORMAL)
					specProbStr = String.valueOf(specProb);
				else
					specProbStr = String.valueOf((float)specProb);
				String pValueStr;
				if(specProb < Float.MIN_NORMAL)
					pValueStr = String.valueOf(pValue);
				else
					pValueStr = String.valueOf((float)pValue);
				String actMethodStr = scorer.getActivationMethodName();
				
				String resultStr =
					specFileName+"\t"
					+specIndex+"\t"
					+scorer.getScanNum()+"\t"
					+actMethodStr+"\t" 
					+scorer.getPrecursorPeak().getMz()+"\t"
					+pmError+"\t"
					+match.getCharge()+"\t"
					+annotationStr+"\t"
					+protein+"\t"
					+match.getDeNovoScore()+"\t"
					+score+"\t"
					+specProbStr+"\t"
					+pValueStr;
				MSGFDBResultGenerator.DBMatch dbMatch = new MSGFDBResultGenerator.DBMatch(specProb, numPeptides, resultStr, match.getScoreDist());		
				gen.add(dbMatch);				
			}
		}
	}
	
	public static void setAminoAcidProbabilities(String databaseFileName, AminoAcidSet aaSet)
	{
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(databaseFileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		long[] aaCount = new long[128];
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith(">"))	// annotation
				continue;
			for(int i=0; i<s.length(); i++)
			{
				char residue = s.charAt(i);
				if(aaSet.getAminoAcid(residue) != null)
					aaCount[residue]++;
			}
		}
		long totalAACount = 0;
		for(AminoAcid aa : aaSet.getAAList(Location.Anywhere))
			if(!aa.isModified())
				totalAACount += aaCount[aa.getResidue()];
		
		boolean success = true;
		for(AminoAcid aa : aaSet.getAllAminoAcidArr())
		{
			long count = aaCount[aa.getUnmodResidue()];
			if(count == 0)
			{
				success = false;
				break;
			}
			aa.setProbability(count/(float)totalAACount);
		}
		
		if(!success)
		{
			System.out.println("Warning: database does not contain all standard amino acids. Probability 0.05 will be used for all amino acids.");
			for(AminoAcid aa : aaSet.getAllAminoAcidArr())
				aa.setProbability(0.05f);
		}
	}
}
