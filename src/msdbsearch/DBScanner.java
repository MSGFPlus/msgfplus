package msdbsearch;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.PriorityQueue;
import java.util.TreeMap;

import parser.BufferedLineReader;

import msgf.DeNovoGraph;
import msgf.GF;
import msgf.GenericDeNovoGraph;
import msgf.MSGFDBResultGenerator;
import msgf.GeneratingFunction;
import msgf.GeneratingFunctionGroup;
import msgf.NominalMass;
import msgf.NominalMassFactory;
import msgf.ScoredSpectrum;
import msgf.Tolerance;
import msscorer.DBScanScorer;
import msscorer.FastScorer;
import msscorer.NewRankScorer;
import msscorer.NewScoredSpectrum;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Enzyme;
import msutil.ModifiedAminoAcid;
import msutil.Peptide;
import msutil.Spectrum;
import msutil.SpectrumAccessorByScanNum;
import msutil.Modification.Location;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class DBScanner extends SuffixArray {

	private int minPeptideLength;
	private int maxPeptideLength;
	private NominalMassFactory factory;
	
	private AminoAcidSet aaSet;
	private double[] aaMass;
	private int[] intAAMass;
	private Enzyme enzyme;
	private int numAllowedC13;
	private int numPeptidesPerSpec;
	
	private SpectrumAccessorByScanNum specMap;
	private Tolerance parentMassTolerance;
	
	private List<Integer> scanNumList;
	private TreeMap<Float,Integer> pepMassScanNumMap;
	private HashMap<Integer, FastScorer> scanNumScorerMap;
	private int[] numDisinctPeptides;
	
	// DB search results
	private HashMap<Integer,PriorityQueue<DatabaseMatch>> scanNumDBMatchMap;

	public DBScanner(
			NominalMassFactory factory,
			SuffixArraySequence sequence, 
			AminoAcidSet aaSet,
			SpectrumAccessorByScanNum specMap,
			List<Integer> scanNumList,
			Tolerance parentMassTolerance,
			int numAllowedC13,
			NewRankScorer scorer,
			Enzyme enzyme,
			int numPeptidesPerSpec,
			int minPeptideLength,
			int maxPeptideLength
			) 
	{
		super(sequence);
		this.factory = factory;
		aaMass = new double[128];
		intAAMass = new int[128];
		for(int i=0; i<aaMass.length; i++)
		{
			aaMass[i] = -1;
			intAAMass[i] = -1;
		}
		for(AminoAcid aa : aaSet)
		{
			aaMass[aa.getResidue()] = aa.getAccurateMass();
			intAAMass[aa.getResidue()] = factory.getMassIndex(aa.getMass());
		}
		
		this.aaSet = aaSet;
		this.parentMassTolerance = parentMassTolerance;
		this.enzyme = enzyme;
		this.numAllowedC13 = numAllowedC13;
		this.scanNumList = scanNumList;
		this.specMap = specMap;
		this.numPeptidesPerSpec = numPeptidesPerSpec;
		this.minPeptideLength = minPeptideLength;
		this.maxPeptideLength = maxPeptideLength;
		
		// set amino acid probabilities
		setAAProbabilities();
		
		pepMassScanNumMap = new TreeMap<Float,Integer>();	// pepMass -> scanNum
		scanNumScorerMap = new HashMap<Integer, FastScorer>();	// scanNumber -> peptideScorer
		
		for(int scanNum : scanNumList)
		{
			Spectrum spec = specMap.getSpectrumByScanNum(scanNum);
			NewScoredSpectrum<NominalMass> scoredSpec = scorer.getScoredSpectrum(spec);
			float peptideMass = spec.getParentMass() - (float)Composition.H2O;
			int peptideMassIndex = factory.getMassIndex(peptideMass);
			if(scorer.supportEdgeScores())
				scanNumScorerMap.put(scanNum, new DBScanScorer(factory, scoredSpec, peptideMassIndex));
			else
				scanNumScorerMap.put(scanNum, new FastScorer(factory, scoredSpec, peptideMassIndex));
			while(pepMassScanNumMap.get(peptideMass) != null)
				peptideMass = Math.nextUp(peptideMass);
			pepMassScanNumMap.put(peptideMass, scanNum);
			
			if(numAllowedC13 > 0 && parentMassTolerance.getToleranceAsDa(peptideMass) <= 2)
			{
				if(numAllowedC13 >= 1)
				{
					float mass1 = peptideMass-(float)Composition.ISOTOPE;
					while(pepMassScanNumMap.get(mass1) != null)
						mass1 = Math.nextUp(mass1);
					pepMassScanNumMap.put(mass1, scanNum);
				}
				
				if(numAllowedC13 >= 2)
				{
					float mass2 = peptideMass-2*(float)Composition.ISOTOPE;
					while(pepMassScanNumMap.get(mass2) != null)
						mass2 = Math.nextUp(mass2);
					pepMassScanNumMap.put(mass2, scanNum);
				}
			}
		}		
		
		// compute the number of distinct peptides
		numDisinctPeptides = new int[maxPeptideLength+2];
		for(int length=minPeptideLength; length<=maxPeptideLength+1; length++)
			numDisinctPeptides[length] = getNumDistinctSeq(length);
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
	
	// duplicated for speeding-up the search
	public void dbSearchCTermEnzyme(int numberOfAllowableNonEnzymaticTermini, boolean verbose)
	{
		scanNumDBMatchMap = new HashMap<Integer,PriorityQueue<DatabaseMatch>>();

		double[] prm = new double[maxPeptideLength+2];
		prm[0] = 0;
		int[] intPRM = new int[maxPeptideLength+2];
		intPRM[0] = 0;
		
		int i = Integer.MAX_VALUE;
		
		int rank = -1;
		boolean enzymaticSearch;
		if(numberOfAllowableNonEnzymaticTermini == 2)
			enzymaticSearch = false;
		else
			enzymaticSearch = true;
		while(indices.hasRemaining()) {
			int index = indices.get();
			rank++;
			if(verbose && rank % 1000000 == 0)
				System.out.println("DBSearch: " + rank/(float)size*100 + "%");
			char nAA = sequence.getCharAt(index);
			int nNET = 0;	// number of non-enzymatic termini
			int lcp = this.neighboringLcps.get(rank);
			if(lcp > i)
				continue;
			if(enzymaticSearch)
			{
				if(Character.isLetter(nAA) && !enzyme.isCleavable(nAA))
				{
					nNET++;
					if(nNET > numberOfAllowableNonEnzymaticTermini)
					{
						i=0;
						continue;
					}
				}
			}
			int nTermAAScore;
			if(nAA == '_' || enzyme.isCleavable(nAA))
				nTermAAScore = enzyme.getNeighboringAACleavageCredit();
			else
				nTermAAScore = enzyme.getNeighboringAACleavagePenalty();
			
			for(i=lcp > 1 ? lcp : 1; i<maxPeptideLength+1 && index+i<size; i++)	// ith character of a peptide
			{
				char residue = sequence.getCharAt(index+i);
				double m = aaMass[residue];
				if(m <= 0)
					break;
				prm[i+1] = prm[i] + m;	// prm[i]: mass of peptide upto ith residue
				intPRM[i+1] = intPRM[i] + intAAMass[residue];
				if(enzymaticSearch && !enzyme.isCleavable(residue))
				{
					if(nNET+1 > numberOfAllowableNonEnzymaticTermini)
						continue;
				}
				if(i < minPeptideLength)
					continue;
				
				float peptideMass = (float)prm[i+1];
				float tolDa = parentMassTolerance.getToleranceAsDa(peptideMass);
				
				Collection<Integer> matchedScanNumList = pepMassScanNumMap.subMap(peptideMass-tolDa, peptideMass+tolDa).values();
				if(matchedScanNumList.size() > 0)
				{
					int peptideCleavageScore;
					if(enzyme.isCleavable(sequence.getCharAt(index+i)))
						peptideCleavageScore = enzyme.getPeptideCleavageCredit();
					else
						peptideCleavageScore = enzyme.getPeptideCleavagePenalty();
					
//					if(sequence.getSubsequence(index+1, index+i+1).equalsIgnoreCase("GNTIEIQGDDAPSLWVYGFSDR"))
//						System.out.println("Debug");
					
					for(Integer scanNum : matchedScanNumList)
					{
						FastScorer scorer = scanNumScorerMap.get(scanNum);
						int score = nTermAAScore + scorer.getScore(prm, intPRM, 2, i+2) + peptideCleavageScore;
						PriorityQueue<DatabaseMatch> prevMatchQueue = scanNumDBMatchMap.get(scanNum);
						if(prevMatchQueue == null)
						{
							prevMatchQueue = new PriorityQueue<DatabaseMatch>();
							scanNumDBMatchMap.put(scanNum, prevMatchQueue);
						}
						else if(prevMatchQueue.size() == 0 || score > prevMatchQueue.peek().getScore())
						{
							if(prevMatchQueue.size() >= this.numPeptidesPerSpec)
								prevMatchQueue.poll();
//							String peptide = sequence.getSubsequence(index+1, index+i+1);
//							char nTermAA = sequence.getCharAt(index);
//							char cTermAA = sequence.getCharAt(index+i+1);
//							annotationStr = nTermAA+"."+peptide+"."+cTermAA;
							prevMatchQueue.add(new DatabaseMatch(index, i+2, score));
						}
					}
				}
			}
		}
		indices.rewind();
		neighboringLcps.rewind();
	}  	
	
	public void dbSearchCTermEnzymeWithMods(int numberOfAllowableNonEnzymaticTermini, boolean verbose)
	{
		scanNumDBMatchMap = new HashMap<Integer,PriorityQueue<DatabaseMatch>>();	// scanNum -> dbHits

		CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aaSet, maxPeptideLength);
		
		int i = Integer.MAX_VALUE;
		
		int rank = -1;
		boolean enzymaticSearch;
		if(numberOfAllowableNonEnzymaticTermini == 2)
			enzymaticSearch = false;
		else
			enzymaticSearch = true;
		while(indices.hasRemaining()) {
			int index = indices.get();
			rank++;
			if(verbose && rank % 1000000 == 0)
				System.out.println("DBSearch: " + rank/(float)size*100 + "%");
			char nAA = sequence.getCharAt(index);
			int nNET = 0;	// number of non-enzymatic termini
			int lcp = this.neighboringLcps.get(rank);
			if(lcp > i)
				continue;
			
			boolean isProteinNTerm; 
//			if(nAA == '_' || (nAA == 'M' && sequence.getCharAt(index-1) == '_'))
			if(nAA == '_')
				isProteinNTerm = true;
			else
				isProteinNTerm = false;
			
			if(enzymaticSearch)
			{
				if(!isProteinNTerm && Character.isLetter(nAA) && !enzyme.isCleavable(nAA))
				{
					nNET++;
					if(nNET > numberOfAllowableNonEnzymaticTermini)
					{
						i=0;
						continue;
					}
				}
			}

			int nTermAAScore;
			if(isProteinNTerm || enzyme.isCleavable(nAA))
				nTermAAScore = enzyme.getNeighboringAACleavageCredit();
			else
				nTermAAScore = enzyme.getNeighboringAACleavagePenalty();
			
			for(i=lcp > 1 ? lcp : 1; i<maxPeptideLength+1 && index+i<size; i++)	// ith character of a peptide
			{
				// TODO: apply C-term mods
				char residue = sequence.getCharAt(index+i);
				if(isProteinNTerm)
				{
					if(candidatePepGrid.addProtNTermResidue(residue) == false)
						break;
				}
				else if(i == 1)
				{
					if(candidatePepGrid.addNTermResidue(residue) == false)
						break;
				}
				else
				{
					if(candidatePepGrid.addResidue(i, residue) == false)
						break;
				}
				if(enzymaticSearch && !enzyme.isCleavable(residue))
				{
					if(nNET+1 > numberOfAllowableNonEnzymaticTermini)
						continue;
				}
				if(i < minPeptideLength)
					continue;
				
				int peptideCleavageScore;
				if(enzyme.isCleavable(sequence.getCharAt(index+i)))
					peptideCleavageScore = enzyme.getPeptideCleavageCredit();
				else
					peptideCleavageScore = enzyme.getPeptideCleavagePenalty();
//				if(sequence.getSubsequence(index+1, index+i+1).equalsIgnoreCase("KALYNET"))
//					System.out.println("Debug");
				
				int enzymeScore = nTermAAScore + peptideCleavageScore;
				
				for(int j=0; j<candidatePepGrid.size(); j++)
				{
					float peptideMass = candidatePepGrid.getPeptideMass(j);
					float tolDa = parentMassTolerance.getToleranceAsDa(peptideMass);
					Collection<Integer> matchedScanNumList = pepMassScanNumMap.subMap(peptideMass-tolDa, peptideMass+tolDa).values();
					if(matchedScanNumList.size() > 0)
					{
						for(Integer scanNum : matchedScanNumList)
						{
							FastScorer scorer = scanNumScorerMap.get(scanNum);
							int score = enzymeScore + scorer.getScore(candidatePepGrid.getPRMGrid()[j], candidatePepGrid.getNominalPRMGrid()[j], 1, i+1); 
							PriorityQueue<DatabaseMatch> prevMatchQueue = scanNumDBMatchMap.get(scanNum);
							if(prevMatchQueue == null)
							{
								prevMatchQueue = new PriorityQueue<DatabaseMatch>();
								scanNumDBMatchMap.put(scanNum, prevMatchQueue);
							}
							if(prevMatchQueue.size() == 0 || score > prevMatchQueue.peek().getScore())
							{
								if(prevMatchQueue.size() >= this.numPeptidesPerSpec)
									prevMatchQueue.poll();
								prevMatchQueue.add(new DatabaseMatch(index, i+2, score).pepSeq(candidatePepGrid.getPeptideSeq(j)));
							}
						}
					}					
				}
			}
		}
		indices.rewind();
		neighboringLcps.rewind();
	}  	
	
	// dupulicated to speed-up the search
	public void dbSearchNoEnzyme()
	{
		assert(false);
	}
	
	// dupulicated to speed-up the search
	public void dbSearchNTermEnzyme(int numberOfAllowableNonEnzymaticTermini)
	{
		assert(false);
	}
	
	public void computeSpecProb(boolean grouped)
	{
		for(int scanNum : scanNumList)
		{
			PriorityQueue<DatabaseMatch> matchQueue = scanNumDBMatchMap.get(scanNum);
			if(matchQueue == null)
				continue;

			Spectrum spec = specMap.getSpectrumByScanNum(scanNum);
			ScoredSpectrum<NominalMass> scoredSpec = scanNumScorerMap.get(scanNum);
			GF<NominalMass> gf = null;
			if(grouped)
			{
				GeneratingFunctionGroup<NominalMass> gfGroup = new GeneratingFunctionGroup<NominalMass>();
				float peptideMass = spec.getParentMass() - (float)Composition.H2O;
				float tolDa = parentMassTolerance.getToleranceAsDa(peptideMass);
				int maxPeptideMassIndex = factory.getMassIndex(peptideMass + tolDa);
				int minPeptideMassIndex = factory.getMassIndex(Math.min(peptideMass-tolDa, peptideMass-numAllowedC13*(float)Composition.ISOTOPE));
				for(int peptideMassIndex = minPeptideMassIndex; peptideMassIndex<=maxPeptideMassIndex; peptideMassIndex++)
				{
					NominalMass pmNode = factory.getInstance(peptideMassIndex);
					if(pmNode != null)
					{
						DeNovoGraph<NominalMass> graph = new GenericDeNovoGraph<NominalMass>(
								factory, 
								factory.getMassFromIndex(peptideMassIndex+factory.getMassIndex((float)Composition.H2O)), 
								Tolerance.ZERO_TOLERANCE, 
								enzyme,
								scoredSpec);
						GeneratingFunction<NominalMass> gfi = new GeneratingFunction<NominalMass>(graph)
						.doNotBacktrack()
						.doNotCalcNumber();
						gfGroup.registerGF(pmNode, gfi);
					}
				}			
				gf = gfGroup;
			}
			else
			{
				DeNovoGraph<NominalMass> graph = new GenericDeNovoGraph<NominalMass>(
						factory, 
						spec.getParentMass(), 
						parentMassTolerance, 
						enzyme,
						scoredSpec);
				gf = new GeneratingFunction<NominalMass>(graph)
					.doNotBacktrack()
					.doNotCalcNumber();
			}
			
			gf.computeGeneratingFunction();
			
			for(DatabaseMatch match : matchQueue)
			{
				match.setDeNovoScore(gf.getMaxScore()-1);
				int score = match.getScore();
				float specProb = gf.getSpectralProbability(score);
				assert(specProb > 0): scanNum + ": " + match.getDeNovoScore()+" "+match.getScore()+" "+specProb; 
				match.setSpecProb(specProb);
				match.setScoreDist(gf.getScoreDist());
			}
		}
	}
	
	public void addDBSearchResults(MSGFDBResultGenerator gen, String specFileName, boolean showTitle)
	{
		for(int scanNum : scanNumList)
		{
			Spectrum spec = specMap.getSpectrumByScanNum(scanNum);
			PriorityQueue<DatabaseMatch> matchQueue = scanNumDBMatchMap.get(scanNum);
			if(matchQueue == null)
				continue;

			ArrayList<DatabaseMatch> matchList = new ArrayList<DatabaseMatch>(matchQueue);
			if(matchList.size() == 0)
				continue;

			DatabaseMatch match = matchList.get(matchList.size()-1);
			// parent mass error
			
//			String peptide = sequence.getSubsequence(index+1, index+i+1);
//			char nTermAA = sequence.getCharAt(index);
//			char cTermAA = sequence.getCharAt(index+i+1);
//			annotationStr = nTermAA+"."+peptide+"."+cTermAA;
			
			int index = match.getIndex();
			int length = match.getLength();
			
//			String peptideStr = sequence.getSubsequence(index+1, index+length-1);
			String peptideStr = match.getPepSeq();
			if(peptideStr == null)
				peptideStr = sequence.getSubsequence(index+1, index+length-1);
			Peptide pep = aaSet.getPeptide(peptideStr);
//			String annotationStr = sequence.getCharAt(index)+"."+peptideStr+"."+sequence.getCharAt(index+length-1);
			String annotationStr = sequence.getCharAt(index)+"."+pep+"."+sequence.getCharAt(index+length-1);
			float expMass = spec.getParentMass();
			float theoMass = pep.getParentMass();
			float pmError = Float.MAX_VALUE;
			for(int numC13=0; numC13<=numAllowedC13; numC13++)
			{
				float error = expMass-theoMass-(float)(Composition.ISOTOPE)*numC13; 
				if(Math.abs(error) < Math.abs(pmError))
					pmError = error;
			}
			if(parentMassTolerance.isTolerancePPM())
				pmError = pmError/theoMass*1e6f;
			
			String protein = getAnnotation(index+1);
			
			int score = match.getScore();
			float specProb = match.getSpecProb();
			int numPeptides = this.numDisinctPeptides[peptideStr.length()+1];
			float pValue = MSGFDBResultGenerator.DBMatch.getPValue(specProb, numPeptides);
			String resultStr =
				specFileName+"\t"
				+scanNum+"\t"
				+(showTitle ? spec.getTitle()+"\t" : "")
				+spec.getPrecursorPeak().getMz()+"\t"
				+pmError+"\t"
				+spec.getCharge()+"\t"
				+annotationStr+"\t"
				+protein+"\t"
				+match.getDeNovoScore()+"\t"
				+score+"\t"
				+specProb+"\t"
				+pValue;
			MSGFDBResultGenerator.DBMatch eFDRMatch = new MSGFDBResultGenerator.DBMatch(specProb, numPeptides, resultStr, match.getScoreDist());		
			gen.add(eFDRMatch);
		}
	}
	
	private void setAAProbabilities()
	{
		long[] aaCount = new long[128];
		for(long i=0; i<super.sequence.getSize(); i++)
		{
			char residue = sequence.getCharAt(i);
			aaCount[residue]++;
		}
		long sum = 0;
		for(AminoAcid aa : aaSet)
			sum += aaCount[aa.getResidue()];
		for(AminoAcid aa : aaSet)
			aa.setProbability(aaCount[aa.getResidue()]/(float)sum);
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
			if(!aa.isVariableModification())
				totalAACount += aaCount[aa.getResidue()];
		for(AminoAcid aa : aaSet.getAllAminoAcidArr())
		{
			aa.setProbability(aaCount[aa.getUnmodResidue()]/(float)totalAACount);
		}
	}
}
