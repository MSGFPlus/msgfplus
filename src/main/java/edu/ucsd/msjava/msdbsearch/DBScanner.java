package edu.ucsd.msjava.msdbsearch;

import java.io.DataInputStream;
import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Map.Entry;

import edu.ucsd.msjava.msgf.DeNovoGraph;
import edu.ucsd.msjava.msgf.FlexAminoAcidGraph;
import edu.ucsd.msjava.msgf.GeneratingFunction;
import edu.ucsd.msjava.msgf.GeneratingFunctionGroup;
import edu.ucsd.msjava.msgf.MSGFDBResultGenerator;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msscorer.SimpleDBSearchScorer;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.SpecKey;
import edu.ucsd.msjava.msutil.Modification.Location;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.sequences.Constants;

public class DBScanner {

	private int minPeptideLength;
	private int maxPeptideLength;
	
	private AminoAcidSet aaSet;
	private double[] aaMass;
	private int[] intAAMass;
	
	private Enzyme enzyme;
	private int numPeptidesPerSpec;

	private final CompactSuffixArray sa;
	private final int size;
	// to scan the database partially
	// Input spectra
	private final ScoredSpectraMap specScanner;
	
	// DB search results
	private Map<SpecKey,PriorityQueue<DatabaseMatch>> specKeyDBMatchMap;
//	private Map<Integer,PriorityQueue<DatabaseMatch>> specIndexDBMatchMap;

	// For output
	private String threadName = "";
	public DBScanner(
			ScoredSpectraMap specScanner,
			CompactSuffixArray sa,
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
//		specIndexDBMatchMap = Collections.synchronizedMap(new HashMap<Integer,PriorityQueue<DatabaseMatch>>());
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
	
	public DBScanner setThreadName(String threadName)
	{
		this.threadName = threadName;
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
			PriorityQueue<DatabaseMatch> queue = entry.getValue();
			
			PriorityQueue<DatabaseMatch> existingQueue = specKeyDBMatchMap.get(entry.getKey());
			if(existingQueue == null)
			{
				existingQueue = new PriorityQueue<DatabaseMatch>();
				specKeyDBMatchMap.put(specKey, existingQueue);
			}
			existingQueue.addAll(queue);
		}
	}
	
	public Map<SpecKey,PriorityQueue<DatabaseMatch>> getSpecKeyDBMatchMap()
	{
		return specKeyDBMatchMap;
	}
	
	public void dbSearchCTermEnzymeNoMod(int numberOfAllowableNonEnzymaticTermini, boolean verbose)
	{
		dbSearch(numberOfAllowableNonEnzymaticTermini, 0, size, verbose);
	}

	public void dbSearchCTermEnzyme(int numberOfAllowableNonEnzymaticTermini, boolean verbose)
	{
		dbSearch(numberOfAllowableNonEnzymaticTermini, 0, size, verbose);
	}
	
	public void dbSearchNTermEnzyme(int numberOfAllowableNonEnzymaticTermini, boolean verbose)
	{
		dbSearch(numberOfAllowableNonEnzymaticTermini, 0, size, verbose);
	}
	
	public void dbSearchNoEnzyme(boolean verbose)
	{
		dbSearch(2, 0, size, verbose);
	}

	public void dbSearch(int numberOfTolerableTermini)
	{
		dbSearch(numberOfTolerableTermini, 0, size, true);
	}
	
	public void dbSearch(int numberOfTolerableTermini, int fromIndex, int toIndex, boolean verbose)
	{
		Map<SpecKey,PriorityQueue<DatabaseMatch>> curSpecKeyDBMatchMap = new HashMap<SpecKey,PriorityQueue<DatabaseMatch>>();
		
		CandidatePeptideGrid candidatePepGrid;
		if(enzyme != null)
			candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aaSet, maxPeptideLength);
		else
			candidatePepGrid = new CandidatePeptideGrid(aaSet, maxPeptideLength);
		
		int i = Integer.MAX_VALUE - 1000;
		
		boolean enzymaticSearch;
		if(numberOfTolerableTermini == 2)
			enzymaticSearch = false;
		else
			enzymaticSearch = true;
		
		int neighboringAACleavageCredit = aaSet.getNeighboringAACleavageCredit();
		int neighboringAACleavagePenalty = aaSet.getNeighboringAACleavagePenalty();
		int peptideCleavageCredit = aaSet.getPeptideCleavageCredit();
		int peptideCleavagePenalty = aaSet.getPeptideCleavagePenalty();
		
		boolean containsCTermMod = aaSet.containsCTermModification();
		
		try {
			DataInputStream indices = new DataInputStream(new BufferedInputStream(new FileInputStream(sa.getIndexFile())));
			indices.skip(CompactSuffixArray.INT_BYTE_SIZE*2+CompactSuffixArray.INT_BYTE_SIZE*fromIndex);	// skip size and id
			
			DataInputStream nlcps = new DataInputStream(new BufferedInputStream(new FileInputStream(sa.getNeighboringLcpFile())));
			nlcps.skip(CompactSuffixArray.INT_BYTE_SIZE*2+fromIndex);
			CompactFastaSequence sequence = sa.getSequence();
		
			boolean isProteinNTerm = true;
			int nTermCleavageScore = 0;
			
			boolean isExtensionAtTheSameIndex;
			int numNonEnzTermini = 0;	// number of non-enzymatic termini
			int numIndices = toIndex-fromIndex;
			
			DatabaseMatch[] prevMatch = new DatabaseMatch[maxPeptideLength+2];
			
			for(int bufferIndex=0; bufferIndex<numIndices; bufferIndex++)
			{
				// Print out the progress
				if(verbose && bufferIndex % 2000000 == 0)
				{
					System.out.print(threadName + ": Database search progress... "); 
					System.out.format("%.1f%% complete\n", bufferIndex/(float)numIndices*100);
				}
				isExtensionAtTheSameIndex = false;
				int index = indices.readInt();
				int lcp = nlcps.readByte();
				if(bufferIndex == 0)
					lcp = 0;

				// For debugging
//				System.out.println(sequence.getSubsequence(index, sequence.getSize()));
//				if(index == 1)
//					System.out.println("Debug");
				// skip redundant peptides
				
				for(int prevMatchIndex=minPeptideLength; prevMatchIndex<prevMatch.length; prevMatchIndex++)
				{
					if(prevMatchIndex<lcp)
					{
						if(prevMatch[prevMatchIndex] != null)
							prevMatch[prevMatchIndex].addIndex(index);
					}
					else
						prevMatch[prevMatchIndex] = null;
					
				}

				if(lcp > i+1 ||
						lcp == i+1 && (enzyme == null || enzyme.isCTerm()))		
				{
					continue;
				}
				else if(lcp == 0)	// preceding aa is changed
				{
					char precedingAA = sequence.getCharAt(index);
					if(precedingAA == Constants.TERMINATOR_CHAR)
						isProteinNTerm = true;
					else
						isProteinNTerm = false;
					
					// determine neighboring N-term score
					if(enzyme == null || enzyme.isNTerm())
					{
						nTermCleavageScore = 0;
					}
					else if(enzyme.isCTerm())
					{
						if(isProteinNTerm || enzyme.isCleavable(precedingAA))// || precedingAA == Constants.INVALID_CHAR)
						{
							nTermCleavageScore = neighboringAACleavageCredit;
							if(enzymaticSearch)
								numNonEnzTermini = 0;
						}
						else
						{
							nTermCleavageScore = neighboringAACleavagePenalty;
							if(enzymaticSearch)
							{
								numNonEnzTermini = 1;
								if(numNonEnzTermini > numberOfTolerableTermini)
								{
									i=0;
									continue;
								}
							}
						}
					}
				}	// end lcp=0
				
				if(lcp == 0)
					i = 1;
				else if(lcp < i+1)
					i = lcp;
					
				for(; i<maxPeptideLength+1 && index+i<size-1; i++)	// ith character of a peptide
				{
					char residue = sequence.getCharAt(index+i);
					boolean isProteinCTerm = false;
					if(i==1)	// N-term residue
					{
						if(enzyme != null && enzyme.isNTerm())
						{
							if(isProteinNTerm || enzyme.isCleavable(residue))	// || sequence.getCharAt(index) == Constants.INVALID_CHAR)
							{
								nTermCleavageScore = peptideCleavageCredit;
								if(enzymaticSearch)
									numNonEnzTermini = 0;
							}
							else
							{
								nTermCleavageScore = peptideCleavagePenalty;
								if(enzymaticSearch)
								{
									numNonEnzTermini = 1;
									if(numNonEnzTermini > numberOfTolerableTermini)
										break;
								}
							}
						}
						
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
					
					if(i < minPeptideLength)
						continue;
					
//					System.out.println(sequence.getSubsequence(index+1, index+i+1));
//					if(sequence.getSubsequence(index, index+i+1).equalsIgnoreCase("RIGAYLFVDMAHVAGLIAAGVYPNPVPHAHVVTSTTHK"))
//						System.out.println("Debug");
					
					int cTermCleavageScore = 0;
					if(enzyme != null)
					{
						char cTermNeighboringResidue = sequence.getCharAt(index+i+1);
						isProteinCTerm = (cTermNeighboringResidue == Constants.TERMINATOR_CHAR);						
						if(enzyme.isCTerm())
						{
//							if(isProteinCTerm || enzyme.isCleavable(residue)) // || cTermNeighboringResidue == Constants.INVALID_CHAR)
							if(enzyme.isCleavable(residue)) // || cTermNeighboringResidue == Constants.INVALID_CHAR)	// changed by Sangtae to avoid SpecProb=0
								cTermCleavageScore = peptideCleavageCredit;
							else
							{
								cTermCleavageScore = peptideCleavagePenalty;
								if(!isProteinCTerm && numNonEnzTermini+1 > numberOfTolerableTermini)
								{
									isExtensionAtTheSameIndex = true;
									continue;
								}
							}
						}
						else if(enzyme.isNTerm())
						{
							if(isProteinCTerm || enzyme.isCleavable(cTermNeighboringResidue)) // || cTermNeighboringResidue == Constants.INVALID_CHAR)
								cTermCleavageScore = neighboringAACleavageCredit;
							else
							{
								cTermCleavageScore = neighboringAACleavagePenalty;
								if(numNonEnzTermini+1 > numberOfTolerableTermini)
								{
									isExtensionAtTheSameIndex = true;
									continue;
								}
							}
						}
					}
					
					int cleavageScore = nTermCleavageScore + cTermCleavageScore;
					
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
							boolean isNTermMetCleaved = candidatePepGrid.isNTermMetCleaved(j);
							int pepLength;
							if(!isNTermMetCleaved)
								pepLength = i;
							else
								pepLength = i-1;
							
							for(SpecKey specKey : matchedSpecKeyList)
							{
								SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(specKey);
//								if(sequence.getSubsequence(index, index+i+1).equalsIgnoreCase("SRDTAIKT"))
//									System.out.println("Debug");
								int score = cleavageScore + scorer.getScore(candidatePepGrid.getPRMGrid(j), candidatePepGrid.getNominalPRMGrid(j), 1, pepLength+1, candidatePepGrid.getNumMods(j)); 
								PriorityQueue<DatabaseMatch> prevMatchQueue = curSpecKeyDBMatchMap.get(specKey);
								if(prevMatchQueue == null)
								{
									prevMatchQueue = new PriorityQueue<DatabaseMatch>();
									curSpecKeyDBMatchMap.put(specKey, prevMatchQueue);
								}
								
								if(prevMatchQueue.size() < this.numPeptidesPerSpec || score == prevMatchQueue.peek().getScore())
								{
									DatabaseMatch dbMatch = new DatabaseMatch(index, (byte)(i+2), score, peptideMass, nominalPeptideMass, specKey.getCharge(), candidatePepGrid.getPeptideSeq(j)).setProteinNTerm(isProteinNTerm).setProteinCTerm(isProteinCTerm);
									prevMatchQueue.add(dbMatch);
									prevMatch[i] = dbMatch;
								}
								else if(prevMatchQueue.size() >= this.numPeptidesPerSpec)
								{
									if(score > prevMatchQueue.peek().getScore())
									{
										while(!prevMatchQueue.isEmpty() && prevMatchQueue.peek().getScore() < score)
											prevMatchQueue.poll();
										DatabaseMatch dbMatch = new DatabaseMatch(index, (byte)(i+2), score, peptideMass, nominalPeptideMass, specKey.getCharge(), candidatePepGrid.getPeptideSeq(j)).setProteinNTerm(isProteinNTerm).setProteinCTerm(isProteinCTerm);
										prevMatchQueue.add(dbMatch);
										prevMatch[i] = dbMatch;
									}
								}
							}
						}					
					}
					isExtensionAtTheSameIndex = true;
				}
			}
			this.addDBMatches(curSpecKeyDBMatchMap);
			indices.close();
			nlcps.close();
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}  		
	
	public void computeSpecProb(boolean storeScoreDist)
	{
		computeSpecProb(storeScoreDist, 0, specScanner.getSpecKeyList().size());
	}
	
	public void computeSpecProb(boolean storeScoreDist, int fromIndex, int toIndex)
	{
		List<SpecKey> specKeyList = specScanner.getSpecKeyList().subList(fromIndex, toIndex);
		
		int numSpecs = toIndex-fromIndex;
		int numProcessedSpecs = 0;
		for(SpecKey specKey : specKeyList)
		{
			numProcessedSpecs++;
			if(numProcessedSpecs % 1000 == 0)
			{
				System.out.print(threadName + ": Computing spectral E-values... "); 
				System.out.format("%.1f%% complete\n", numProcessedSpecs/(float)numSpecs*100);
			}
			
			PriorityQueue<DatabaseMatch> matchQueue = specKeyDBMatchMap.get(specKey);
			if(matchQueue == null)
				continue;

			int specIndex = specKey.getSpecIndex();
			
			boolean useProtNTerm = false;
			boolean useProtCTerm = false;
			int minScore = Integer.MAX_VALUE;
			for(DatabaseMatch m : matchQueue)
			{
				if(m.isProteinNTerm())
					useProtNTerm = true;
				if(m.isProteinCTerm())
					useProtCTerm = true;
				if(m.getScore() < minScore)
					minScore = m.getScore();
			}
			
			GeneratingFunctionGroup<NominalMass> gf = new GeneratingFunctionGroup<NominalMass>();
			SimpleDBSearchScorer<NominalMass> scoredSpec = specScanner.getSpecKeyScorerMap().get(specKey);
			float peptideMass = scoredSpec.getPrecursorPeak().getMass() - (float)Composition.H2O;
			int nominalPeptideMass = NominalMass.toNominalMass(peptideMass);
			int minNominalPeptideMass = nominalPeptideMass - specScanner.getMaxNum13C();
			int maxNominalPeptideMass = nominalPeptideMass + specScanner.getMinNum13C();
			
			float tolDaLeft = specScanner.getLeftParentMassTolerance().getToleranceAsDa(peptideMass);
			float tolDaRight = specScanner.getRightParentMassTolerance().getToleranceAsDa(peptideMass);
			int maxPeptideMassIndex, minPeptideMassIndex;
			
			maxPeptideMassIndex = maxNominalPeptideMass + Math.round(tolDaLeft-0.4999f);
			minPeptideMassIndex = minNominalPeptideMass - Math.round(tolDaRight-0.4999f);
			
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
	
	// for MS-GFDB
	public synchronized void addDBSearchResults(List<MSGFDBResultGenerator.DBMatch> gen, String specFileName, boolean replicateMergedResults)
	{
		Map<Integer,PriorityQueue<DatabaseMatch>> specIndexDBMatchMap = new HashMap<Integer,PriorityQueue<DatabaseMatch>>();
		
		Iterator<Entry<SpecKey, PriorityQueue<DatabaseMatch>>> itr = specKeyDBMatchMap.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<SpecKey, PriorityQueue<DatabaseMatch>> entry = itr.next();
			SpecKey specKey = entry.getKey();
			PriorityQueue<DatabaseMatch> matchQueue = entry.getValue();
			if(matchQueue == null || matchQueue.size() == 0)
				continue;
			
			int specIndex = specKey.getSpecIndex();
			PriorityQueue<DatabaseMatch> existingQueue = specIndexDBMatchMap.get(specIndex);
			if(existingQueue == null)
			{
				existingQueue = new PriorityQueue<DatabaseMatch>(this.numPeptidesPerSpec, new DatabaseMatch.SpecProbComparator());
				specIndexDBMatchMap.put(specIndex, existingQueue);
			}
			
			for(DatabaseMatch match : matchQueue)
			{
				if(existingQueue.size() < this.numPeptidesPerSpec)
				{
					existingQueue.add(match);
				}
				else if(existingQueue.size() >= this.numPeptidesPerSpec)
				{
					if(match.getSpecEValue() < existingQueue.peek().getSpecEValue())
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
				int charge = match.getCharge();
				
				String peptideStr = match.getPepSeq();
				if(peptideStr == null)
					peptideStr = sa.getSequence().getSubsequence(index+1, index+length-1);
				Peptide pep = aaSet.getPeptide(peptideStr);
				String annotationStr = sa.getSequence().getCharAt(index)+"."+pep+"."+sa.getSequence().getCharAt(index+length-1);
				SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(new SpecKey(specIndex, charge));
				ArrayList<Integer> specIndexList = specScanner.getSpecKey(specIndex, charge).getSpecIndexList();
				if(specIndexList == null)
				{
					specIndexList = new ArrayList<Integer>();
					specIndexList.add(specIndex);
				}
				
				float expMass = scorer.getPrecursorPeak().getMass();
//				float theoMass = pep.getParentMass();
				float peptideMass = match.getPeptideMass();
				float pmError = Float.MAX_VALUE;
				float theoMass = peptideMass + (float)Composition.H2O;

				int deltaNominalMass = 0;
				for(int delta=specScanner.getMinNum13C(); delta<=specScanner.getMaxNum13C(); delta++)
				{
					float error = expMass-theoMass-(float)(Composition.ISOTOPE)*delta; 
					if(Math.abs(error) < Math.abs(pmError))
					{
						pmError = error;
						deltaNominalMass = delta;
					}
				}
				if(specScanner.getRightParentMassTolerance().isTolerancePPM())
					pmError = pmError/theoMass*1e6f;
				
				String protein = sa.getAnnotation(index+1);
				
				int score = match.getScore();
				double specProb = match.getSpecEValue();
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

				if(!replicateMergedResults)
				{
					StringBuffer specIndexStrBuf = new StringBuffer();
					StringBuffer scanNumStrBuf = new StringBuffer();
					StringBuffer actMethodStrBuf = new StringBuffer();
					specIndexStrBuf.append(specIndexList.get(0));
					actMethodStrBuf.append(scorer.getActivationMethodArr()[0]);
					scanNumStrBuf.append(scorer.getScanNumArr()[0]);
					for(int j=1; j<scorer.getActivationMethodArr().length; j++)
					{
						specIndexStrBuf.append("/"+specIndexList.get(j));
						scanNumStrBuf.append("/"+scorer.getScanNumArr()[j]);
						actMethodStrBuf.append("/"+scorer.getActivationMethodArr()[j]);
					}
					
					String resultStr =
						specFileName+"\t"
						+specIndexStrBuf.toString()+"\t"
						+scanNumStrBuf.toString()+"\t"
						+actMethodStrBuf.toString()+"\t" 
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
				else
				{
					for(int j=0; j<scorer.getActivationMethodArr().length; j++)
					{
						String resultStr =
							specFileName+"\t"
							+specIndexList.get(j)+"\t"
							+scorer.getScanNumArr()[j]+"\t"
							+scorer.getActivationMethodArr()[j]+"\t" 
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
