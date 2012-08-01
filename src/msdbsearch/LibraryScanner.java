package msdbsearch;

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

import parser.BufferedLineReader;

import msgf.DeNovoGraph;
import msgf.FlexAminoAcidGraph;
import msgf.MSGFDBResultGenerator;
import msgf.GeneratingFunction;
import msgf.GeneratingFunctionGroup;
import msgf.NominalMass;
import msscorer.SimpleDBSearchScorer;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Modification;
import msutil.SpecKey;
import msutil.Modification.Location;

public class LibraryScanner {

	private final int MAX_LIBRARY_PEPTIDE_LENGTH = 100;

	private double[] aaMass;
	private int[] intAAMass;

	private int numPeptidesPerSpec;

	// Input spectra
	private final ScoredSpectraMap specScanner;

	// DB search results
	private Map<SpecKey,PriorityQueue<LibraryMatch>> specKeyDBMatchMap;
	private Map<Integer,PriorityQueue<LibraryMatch>> specIndexDBMatchMap;
	private int numPeptidesInLib = 0;

	// For output
	private String threadName = "";
	public LibraryScanner(
			ScoredSpectraMap specScanner,
			int numPeptidesPerSpec
	) 
	{
		this.specScanner = specScanner;
		this.numPeptidesPerSpec = numPeptidesPerSpec;

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

		specKeyDBMatchMap = Collections.synchronizedMap(new HashMap<SpecKey,PriorityQueue<LibraryMatch>>());
		specIndexDBMatchMap = Collections.synchronizedMap(new HashMap<Integer,PriorityQueue<LibraryMatch>>());
	}

	public LibraryScanner setThreadName(String threadName)
	{
		this.threadName = threadName;
		return this;
	}

	public synchronized void addDBMatches(Map<SpecKey,PriorityQueue<LibraryMatch>> map)
	{
		if(map == null)
			return;
		Iterator<Entry<SpecKey, PriorityQueue<LibraryMatch>>> itr = map.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<SpecKey, PriorityQueue<LibraryMatch>> entry = itr.next();
			SpecKey specKey = entry.getKey(); 
			PriorityQueue<LibraryMatch> queue = specKeyDBMatchMap.get(entry.getKey());
			if(queue == null)
			{
				queue = new PriorityQueue<LibraryMatch>();
				specKeyDBMatchMap.put(specKey, queue);
			}
			for(LibraryMatch match : entry.getValue())
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

	public void libSearch(String libFilePath, boolean verbose)
	{
//		Map<SpecKey,PriorityQueue<LibraryMatch>> targetSpecKeyDBMatchMap = libSearch(libFilePath, false, true);
//		Map<SpecKey,PriorityQueue<LibraryMatch>> decoySpecKeyDBMatchMap = libSearch(libFilePath, true, true);
//		this.addDBMatches(targetSpecKeyDBMatchMap);
//		this.addDBMatches(decoySpecKeyDBMatchMap);
		
		this.addDBMatches(libSearchPlain(libFilePath, true));
	}
	
	// Reads peptide variants from sptxt file
	private Map<SpecKey,PriorityQueue<LibraryMatch>> libSearchPlain(String libFilePath, boolean verbose)
	{
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(libFilePath);
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		
		Map<SpecKey,PriorityQueue<LibraryMatch>> curSpecKeyDBMatchMap = new HashMap<SpecKey,PriorityQueue<LibraryMatch>>();

		String s;

		int numPeptides = 0;
		
		String pepStr = null;
		int pepLength = 0;
		int charge = -1;
		
		while((s=in.readLine()) != null)
		{
			if(s.trim().length() == 0)
				continue;
			else if(s.startsWith("Name:"))
			{
				numPeptides++;
				// Print out the progress
				if(numPeptides % 100000 == 100000-1)
				{
					System.out.print(threadName + ": Database search progress... "); 
					System.out.format("%dE5 peptides complete\n", numPeptides/100000);
				}
				// Name: AAAAA...GAK/2
				String[] token = s.split("\\s+");
				String name = token[1];
				charge = Integer.parseInt(name.substring(name.lastIndexOf('/')+1));
				StringBuffer pepBuf = new StringBuffer();
				for(int i=0; i<name.length(); i++)
				{
					if(Character.isUpperCase(name.charAt(i)))
					{
						pepBuf.append(name.charAt(i));
					}
				}
				pepLength = pepBuf.length();
				pepStr = pepBuf.toString();
			}
			else if(s.startsWith("Comment:"))
			{
				int numMods = -1;
				double[] modMass = new double[MAX_LIBRARY_PEPTIDE_LENGTH]; // 1-based
				int[] nominalModMass = new int[MAX_LIBRARY_PEPTIDE_LENGTH]; // 1-based
				String[] modResidues = new String[MAX_LIBRARY_PEPTIDE_LENGTH]; // 1-based
				String protein = null;
				
				// Comment:
				String[] token = s.split("\\s+");
				for(int i=0; i<token.length; i++)
				{
					String curToken = token[i]; 
					
					// modification
					if(curToken.startsWith("Mods="))
					{
						String[] modToken = curToken.split("[=/]");
						numMods = Integer.parseInt(modToken[1]);
						for(int j=2; j<modToken.length; j++)
						{
							String[] mod = modToken[j].split(",");
							int location = Integer.parseInt(mod[0]);	// 0-base
							if(location == -1)
								location = 0;
							
							String modName = mod[2];
							double deltaMass = modTable.get(modName);
							modMass[location+1] = deltaMass;
							nominalModMass[location+1] = NominalMass.toNominalMass((float)deltaMass);
							modResidues[location+1] = modResidueTable.get(modName);
						}					
					}
					// protein
					else if(curToken.startsWith("Protein="))
					{
						String[] protToken = curToken.split("[=/]");
						protein = protToken[2];
					}
				}				
				
				// always 0 at index 0, mass of ith prefix at index i
				int[] nominalPRM = new int[MAX_LIBRARY_PEPTIDE_LENGTH];
				double[] prm = new double[MAX_LIBRARY_PEPTIDE_LENGTH];

				nominalPRM[0] = 0;
				prm[0] = 0;
				StringBuffer peptideOutput = new StringBuffer();
				for(int i=0; i<pepLength; i++)	// ith character of a peptide (base 0)
				{
					char residue = pepStr.charAt(i);
					nominalPRM[i+1] = nominalPRM[i] + intAAMass[residue] + nominalModMass[i+1];
					prm[i+1] = prm[i] + aaMass[residue] + modMass[i+1];
					peptideOutput.append(pepStr.charAt(i)+(modResidues[i+1] == null ? "" : modResidues[i+1]));
				}

				float peptideMass = (float)prm[pepLength];
				int nominalPeptideMass = nominalPRM[pepLength];
				float tolDaLeft = specScanner.getLeftParentMassTolerance().getToleranceAsDa(peptideMass);
				float tolDaRight = specScanner.getRightParentMassTolerance().getToleranceAsDa(peptideMass);

				double leftThr = (double)(peptideMass - tolDaRight);
				double rightThr = (double)(peptideMass + tolDaLeft);
				Collection<SpecKey> matchedSpecKeyList = specScanner.getPepMassSpecKeyMap().subMap(leftThr, rightThr).values();
				for(SpecKey specKey : matchedSpecKeyList)
				{
					if(charge != specKey.getCharge())
						continue;
					SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(specKey);
					int score = scorer.getScore(prm, nominalPRM, 1, pepLength+1, numMods); 
					PriorityQueue<LibraryMatch> prevMatchQueue = curSpecKeyDBMatchMap.get(specKey);
					if(prevMatchQueue == null)
					{
						prevMatchQueue = new PriorityQueue<LibraryMatch>();
						curSpecKeyDBMatchMap.put(specKey, prevMatchQueue);
					}
					if(prevMatchQueue.size() < this.numPeptidesPerSpec)
					{
						prevMatchQueue.add(new LibraryMatch(score, peptideMass, nominalPeptideMass, charge, peptideOutput.toString(), protein));
					}
					else if(prevMatchQueue.size() >= this.numPeptidesPerSpec)
					{
						if(score > prevMatchQueue.peek().getScore())
						{
							prevMatchQueue.poll();
							prevMatchQueue.add(new LibraryMatch(score, peptideMass, nominalPeptideMass, charge, peptideOutput.toString(), protein));
						}
					}
				}				
			}
		}

		if(in != null)
		{
			try {
				in.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return curSpecKeyDBMatchMap;
	}  		
	
	// Reads peptide variants from sptxt file
	private Map<SpecKey,PriorityQueue<LibraryMatch>> libSearch(String libFilePath, boolean isDecoy, boolean verbose)
	{
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(libFilePath);
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		
		Map<SpecKey,PriorityQueue<LibraryMatch>> curSpecKeyDBMatchMap = new HashMap<SpecKey,PriorityQueue<LibraryMatch>>();

		String s;

		int numPeptides = 0;
		while((s=in.readLine()) != null)
		{
			if(!s.startsWith("Comment:"))
				continue;

			// Print out the progress
			if(verbose && numPeptides > 0 && numPeptides % 100000 == 0)
			{
				System.out.print(threadName + ": Database search progress... "); 
				System.out.format("%dE5 peptides complete\n", numPeptides/100000);
			}

			// these should be filled by parsing the file
			String pepStr = null;
			int pepLength = 0;
			int charge = -1;
			int numMods = -1;
			double[] modMass = new double[MAX_LIBRARY_PEPTIDE_LENGTH]; // 1-based
			int[] nominalModMass = new int[MAX_LIBRARY_PEPTIDE_LENGTH]; // 1-based
			String[] modResidues = new String[MAX_LIBRARY_PEPTIDE_LENGTH]; // 1-based
			String protein = null;

			String[] token = s.split("\\s+");
			for(int i=0; i<token.length; i++)
			{
				String curToken = token[i]; 
				if(curToken.startsWith("Fullname="))
				{
					String[] pepToken = curToken.split("[=./]");
					pepStr = pepToken[2];
					pepStr = pepStr.replaceAll("M\\(O\\)","M");
					pepLength = pepStr.length();
					charge = Integer.parseInt(pepToken[4]);
					
					if(isDecoy)	
					{
						// e.g. QGACK -> QCAGK
						StringBuffer reversePepStr = new StringBuffer();
						reversePepStr.append(pepStr.charAt(0));
						for(int j=pepLength-2; j>=1; j--)
							reversePepStr.append(pepStr.charAt(j));
						reversePepStr.append(pepStr.charAt(pepLength-1));
						pepStr = reversePepStr.toString();
					}
				}
				
				// modification
				else if(curToken.startsWith("Mods="))
				{
					String[] modToken = curToken.split("[=/]");
					numMods = Integer.parseInt(modToken[1]);
					for(int j=2; j<modToken.length; j++)
					{
						String[] mod = modToken[j].split(",");
						int location = Integer.parseInt(mod[0]);	// 0-base
						if(location == -1)
							location = 0;
						
						if(isDecoy)
						{
							if(location > 0 && location < pepLength-1)
								location = pepLength-1-location;
						}
						
						String modName = mod[2];
						double deltaMass = modTable.get(modName);
						modMass[location+1] = deltaMass;
						nominalModMass[location+1] = NominalMass.toNominalMass((float)deltaMass);
						modResidues[location+1] = modResidueTable.get(modName);
					}					
				}
				// protein
				else if(curToken.startsWith("Protein="))
				{
					String[] protToken = curToken.split("[=/]");
					protein = protToken[2];
					if(isDecoy)
						protein = "DECOY_" + protein;
				}
			}

			numPeptides++;
			
			// always 0 at index 0, mass of ith prefix at index i
			int[] nominalPRM = new int[MAX_LIBRARY_PEPTIDE_LENGTH];
			double[] prm = new double[MAX_LIBRARY_PEPTIDE_LENGTH];

			nominalPRM[0] = 0;
			prm[0] = 0;
			StringBuffer peptideOutput = new StringBuffer();
			for(int i=0; i<pepLength; i++)	// ith character of a peptide (base 0)
			{
				char residue = pepStr.charAt(i);
				nominalPRM[i+1] = nominalPRM[i] + intAAMass[residue] + nominalModMass[i+1];
				prm[i+1] = prm[i] + aaMass[residue] + modMass[i+1];
				peptideOutput.append(pepStr.charAt(i)+(modResidues[i+1] == null ? "" : modResidues[i+1]));
			}

			float peptideMass = (float)prm[pepLength];
			int nominalPeptideMass = nominalPRM[pepLength];
			float tolDaLeft = specScanner.getLeftParentMassTolerance().getToleranceAsDa(peptideMass);
			float tolDaRight = specScanner.getRightParentMassTolerance().getToleranceAsDa(peptideMass);

			double leftThr = (double)(peptideMass - tolDaRight);
			double rightThr = (double)(peptideMass + tolDaLeft);
			Collection<SpecKey> matchedSpecKeyList = specScanner.getPepMassSpecKeyMap().subMap(leftThr, rightThr).values();
			for(SpecKey specKey : matchedSpecKeyList)
			{
				if(charge != specKey.getCharge())
					continue;
				SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(specKey);
				int score = scorer.getScore(prm, nominalPRM, 1, pepLength+1, numMods); 
				PriorityQueue<LibraryMatch> prevMatchQueue = curSpecKeyDBMatchMap.get(specKey);
				if(prevMatchQueue == null)
				{
					prevMatchQueue = new PriorityQueue<LibraryMatch>();
					curSpecKeyDBMatchMap.put(specKey, prevMatchQueue);
				}
				if(prevMatchQueue.size() < this.numPeptidesPerSpec)
				{
					prevMatchQueue.add(new LibraryMatch(score, peptideMass, nominalPeptideMass, charge, peptideOutput.toString(), protein));
				}
				else if(prevMatchQueue.size() >= this.numPeptidesPerSpec)
				{
					if(score > prevMatchQueue.peek().getScore())
					{
						prevMatchQueue.poll();
						prevMatchQueue.add(new LibraryMatch(score, peptideMass, nominalPeptideMass, charge, peptideOutput.toString(), protein));
					}
				}
			}
		}

		if(in != null)
		{
			try {
				in.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return curSpecKeyDBMatchMap;
	}  		

	public void computeSpecProb()
	{
		computeSpecProb(0, specScanner.getSpecKeyList().size());
	}

	public void computeSpecProb(int fromIndex, int toIndex)
	{
		List<SpecKey> specKeyList = specScanner.getSpecKeyList().subList(fromIndex, toIndex);

		int numSpecs = toIndex-fromIndex;
		int numProcessedSpecs = 0;
		for(SpecKey specKey : specKeyList)
		{
			numProcessedSpecs++;
			if(numProcessedSpecs % 1000 == 0)
			{
				System.out.print(threadName + ": Computing spectral probabilities... "); 
				System.out.format("%.1f%% complete\n", numProcessedSpecs/(float)numSpecs*100);
			}

			PriorityQueue<LibraryMatch> matchQueue = specKeyDBMatchMap.get(specKey);
			if(matchQueue == null)
				continue;

			int specIndex = specKey.getSpecIndex();
			int minScore = Integer.MAX_VALUE;
			for(LibraryMatch m : matchQueue)
			{
				if(m.getScore() < minScore)
					minScore = m.getScore();
			}

			GeneratingFunctionGroup<NominalMass> gf = new GeneratingFunctionGroup<NominalMass>();
			SimpleDBSearchScorer<NominalMass> scoredSpec = specScanner.getSpecKeyScorerMap().get(specKey);
			float peptideMass = scoredSpec.getPrecursorPeak().getMass() - (float)Composition.H2O;
			int nominalPeptideMass = NominalMass.toNominalMass(peptideMass);
			int minNominalPeptideMass = nominalPeptideMass + specScanner.getMinNum13C();
			int maxNominalPeptideMass = nominalPeptideMass + specScanner.getMaxNum13C();
			
			float tolDaLeft = specScanner.getLeftParentMassTolerance().getToleranceAsDa(peptideMass);
			float tolDaRight = specScanner.getRightParentMassTolerance().getToleranceAsDa(peptideMass);
			int maxPeptideMassIndex, minPeptideMassIndex;
			
			maxPeptideMassIndex = minNominalPeptideMass + Math.round(tolDaLeft-0.4999f);
			minPeptideMassIndex = maxNominalPeptideMass - Math.round(tolDaRight-0.4999f);

			for(int peptideMassIndex = minPeptideMassIndex; peptideMassIndex<=maxPeptideMassIndex; peptideMassIndex++)
			{
				DeNovoGraph<NominalMass> graph = new FlexAminoAcidGraph(
						aaSet, 
						peptideMassIndex,
						null,
						scoredSpec,
						true,
						false
				);

				GeneratingFunction<NominalMass> gfi = new GeneratingFunction<NominalMass>(graph)
				.doNotBacktrack()
				.doNotCalcNumber();
				gfi.setUpScoreThreshold(minScore);
				gf.registerGF(graph.getPMNode(), gfi);
			}

			boolean isGFComputed = gf.computeGeneratingFunction();

			for(LibraryMatch match : matchQueue)
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
				}
			}
		}
	}

	public synchronized void addLibSearchResults(List<MSGFDBResultGenerator.DBMatch> gen, String specFileName)
	{
		Iterator<Entry<SpecKey, PriorityQueue<LibraryMatch>>> itr = specKeyDBMatchMap.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<SpecKey, PriorityQueue<LibraryMatch>> entry = itr.next();
			SpecKey specKey = entry.getKey();
			PriorityQueue<LibraryMatch> matchQueue = entry.getValue();
			if(matchQueue == null || matchQueue.size() == 0)
				continue;

			int specIndex = specKey.getSpecIndex();
			PriorityQueue<LibraryMatch> existingQueue = specIndexDBMatchMap.get(specIndex);
			if(existingQueue == null)
			{
				existingQueue = new PriorityQueue<LibraryMatch>(this.numPeptidesPerSpec, new Match.SpecProbComparator());
				specIndexDBMatchMap.put(specIndex, existingQueue);
			}

			for(LibraryMatch match : matchQueue)
			{
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

		Iterator<Entry<Integer, PriorityQueue<LibraryMatch>>> itr2 = specIndexDBMatchMap.entrySet().iterator();
		while(itr2.hasNext())
		{
			Entry<Integer, PriorityQueue<LibraryMatch>> entry = itr2.next();
			int specIndex = entry.getKey();
			PriorityQueue<LibraryMatch> matchQueue = entry.getValue();
			if(matchQueue == null)
				continue;

			ArrayList<LibraryMatch> matchList = new ArrayList<LibraryMatch>(matchQueue);
			if(matchList.size() == 0)
				continue;

			for(int i=matchList.size()-1; i>=0; --i)
			{
				LibraryMatch match = matchList.get(i);

				if(match.getDeNovoScore() < 0)
					continue;

				int charge = match.getCharge();

				String annotationStr = match.getPepSeq();
				SimpleDBSearchScorer<NominalMass> scorer = specScanner.getSpecKeyScorerMap().get(new SpecKey(specIndex, charge));
				ArrayList<Integer> specIndexList = specScanner.getSpecKey(specIndex, charge).getSpecIndexList();
				if(specIndexList == null)
				{
					specIndexList = new ArrayList<Integer>();
					specIndexList.add(specIndex);
				}

				float expMass = scorer.getPrecursorPeak().getMass();
				float peptideMass = match.getPeptideMass();
				float theoMass = peptideMass + (float)Composition.H2O;
				float pmError = Float.MAX_VALUE;

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

				String protein = match.getProtein();	// current no protein id is assigned

				int score = match.getScore();
				double specProb = match.getSpecProb();
				double pValue = MSGFDBResultGenerator.DBMatch.getPValue(specProb, numPeptidesInLib);
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
				MSGFDBResultGenerator.DBMatch dbMatch = new MSGFDBResultGenerator.DBMatch(specProb, numPeptidesInLib, resultStr, match.getScoreDist());		
				gen.add(dbMatch);				
			}
		}
	}

	private static HashMap<String,Double> modTable;
	private static HashMap<String,String> modResidueTable;
	private static AminoAcidSet aaSet;

	static {
		modTable = new HashMap<String,Double>();
		//		modTable.put("Carbamidomethyl", Modification.get("Carbamidomethylation").getAccurateMass());
		modTable.put("Carbamidomethyl", 0.);
		modTable.put("Pyro-carbamidomethyl", Modification.get("PyroCarbamidomethyl").getAccurateMass());
		modTable.put("Oxidation", Modification.get("Oxidation").getAccurateMass());
		modTable.put("Acetyl", Modification.get("Acetylation").getAccurateMass());
		modTable.put("Gln->pyro-Glu", Modification.get("PyrogluQ").getAccurateMass());
		modTable.put("Glu->pyro-Glu", Modification.get("PyrogluE").getAccurateMass());

		modResidueTable = new HashMap<String,String>();
		//		modResidueTable.put("Carbamidomethyl", String.format("%.3f", "+"+Modification.get("Carbamidomethylation").getMass()));
		modResidueTable.put("Carbamidomethyl", "");
		modResidueTable.put("Pyro-carbamidomethyl", String.format("%.3f", Modification.get("PyroCarbamidomethyl").getMass()));
		modResidueTable.put("Oxidation", String.format("+%.3f", Modification.get("Oxidation").getMass()));
		modResidueTable.put("Acetyl", String.format("+%.3f", Modification.get("Acetylation").getMass()));
		modResidueTable.put("Gln->pyro-Glu", String.format("%.3f", Modification.get("PyrogluQ").getMass()));
		modResidueTable.put("Glu->pyro-Glu", String.format("%.3f", Modification.get("PyrogluE").getMass()));

		// set up aaSet
		ArrayList<Modification.Instance> mods = new ArrayList<Modification.Instance>();
		mods.add(new Modification.Instance(Modification.get("Carbamidomethylation"), 'C').fixedModification());
		mods.add(new Modification.Instance(Modification.get("PyroCarbamidomethyl"), 'C', Location.N_Term));
		mods.add(new Modification.Instance(Modification.get("Oxidation"), 'M', Location.Anywhere));
		mods.add(new Modification.Instance(Modification.get("Acetylation"), '*', Location.N_Term));
		mods.add(new Modification.Instance(Modification.get("PyrogluQ"), 'Q', Location.N_Term));
		mods.add(new Modification.Instance(Modification.get("PyrogluE"), 'E', Location.N_Term));

		aaSet = AminoAcidSet.getAminoAcidSet(mods);
	}

}
