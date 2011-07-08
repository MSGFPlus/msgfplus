package msscorer;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.TreeMap;
import java.util.TreeSet;

import msgf.Histogram;
import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.IonType;
import msutil.Pair;
import msutil.Peak;
import msutil.Peptide;
import msutil.SpectraContainer;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

/**
 * This supports both low and high accuracy fragment ions.
 * @author sangtaekim
 * @deprecated
 */
public class NewScoringParameterGenerator extends NewRankScorer {
	private static final float MIN_PRECURSOR_OFFSET = -300f;	// for precursors
	private static final float MAX_PRECURSOR_OFFSET = 30f;
	private static final int MIN_NUM_SPECTRA_PER_PARTITION = 400;	// 400
	private static final int MIN_NUM_SPECTRA_FOR_PRECURSOR_OFF = 150;
	
	private static final float SIGNAL_TO_NOISE_RATIO_PRECURSOR_FILTER = 20f;
	private static final float SIGNAL_TO_NOISE_RATIO = 10f;
	private static final float NOISE_LEVEL_RATIO = 3f;
	private static final int NUM_NOISE_IONS = 10;
	
	private static final int MAX_RANK = 150;
	private static final int NUM_SEGMENTS_PER_SPECTRUM = 2;	// 2
	
    private static final int[] smoothingRanks = {3, 5, 10, 20, 50, Integer.MAX_VALUE}; //Ranks around which smoothing occurs
    private static final int[] smoothingWindowSize = {0, 1, 2, 3, 4, 5}; //Smoothing windows for each smoothing rank
	
	protected static final int MAX_CHARGE = 20;
 	
	public static void main(String argv[])
	{
		File specFile = null;
		File outputFile = null;
		boolean isText = false;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();;
		Tolerance mme = new Tolerance(0.5f);
		boolean useError = false;
		int numSpecsPerPeptide = 1;
		
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length)
				printUsageAndExit("Illegal parameter!");
			if(argv[i].equalsIgnoreCase("-i"))
			{
				specFile = new File(argv[i+1]);
				if(!specFile.exists())
				{
					printUsageAndExit(argv[i+1]+" doesn't exist.");
				}
				int posDot = specFile.getName().lastIndexOf('.');
				if(posDot >= 0)
				{
					String extension = specFile.getName().substring(posDot);
					if(!extension.equalsIgnoreCase(".mgf"))
						printUsageAndExit("Illegal spectrum format: " + argv[i+1]);
				}
				else
					printUsageAndExit("Illegal spectrum format: " + argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-o"))
			{
				outputFile = new File(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-t"))
			{
				outputFile = new File(argv[i+1]);
				isText = true;
			}
			else if(argv[i].equalsIgnoreCase("-fixMod"))
			{
				// 0: No mod, 1: Carbamidomethyl C, 2: Carboxymethyl C
				if(argv[i+1].equalsIgnoreCase("0"))
					aaSet = AminoAcidSet.getStandardAminoAcidSet();
				else if(argv[i+1].equalsIgnoreCase("1"))
					aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
				else if(argv[i+1].equalsIgnoreCase("2"))
					aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarboxymethylatedCys();
				else
					printUsageAndExit("Illigal -fixMod parameter: " + argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-mme"))
			{
				mme = Tolerance.parseToleranceStr(argv[i+1]);
				if(mme == null)
				{
					printUsageAndExit("Illegal mme value: " + argv[i+1]);
				}
			}
			else if(argv[i].equalsIgnoreCase("-pep"))
			{
				numSpecsPerPeptide = Integer.parseInt(argv[i+1]);
			}
			else
				printUsageAndExit("Illegal parameters!");
		}
		if(specFile == null)
			printUsageAndExit("missing annotatedMgfFileName!");
		if(outputFile == null)
			printUsageAndExit("missing outputFileName!");
		
		generateParameters(specFile, numSpecsPerPeptide, outputFile, aaSet, isText, false, mme, useError);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.err.println(message);
		System.out.println("usage: java -Xmx2000M -cp MSGF.jar msscorer.ScoringParameterGenerator\n" +
				"\t-i annotatedMgfFileName (*.mgf)\n" +
				"\t-o outputFileName (e.g. CID_Tryp.param)\n" +
				"\t[-mme maximum mass error] (e.g. 0.5Da or 30ppm, default: 0.5Da)\n" +
				"\t[-fixMod 0/1/2] (0: NoCysteineProtection, 1: CarbamidomethyC (default), 2: CarboxymethylC)\n" +
				"\t[-pep numPeptidesPerSpec]  (default: 1)" 
		);
		
		System.exit(-1);
	}
	
	public static void generateParameters(File specFile, int numSpecsPerPeptide, File outputFile, AminoAcidSet aaSet, boolean isText, boolean verbose, Tolerance mme, boolean useError)
	{
		SpectraContainer container = new SpectraContainer(specFile.getPath(), new MgfSpectrumParser());
		
		NewScoringParameterGenerator gen = null;
		// multiple spectra with the same peptide -> one spec per peptide
		HashMap<String,ArrayList<Spectrum>> pepSpecMap = new HashMap<String,ArrayList<Spectrum>>();
		SpectraContainer specContOnePerPep = new SpectraContainer();
		for(Spectrum spec : container)
		{
			String pep = spec.getAnnotationStr()+":"+spec.getCharge();
			if(pep != null && pep.length() > 0)
			{
				ArrayList<Spectrum> specList = pepSpecMap.get(pep);
				if(specList == null)
				{
					specList = new ArrayList<Spectrum>();
					pepSpecMap.put(pep, specList);
				}
				if(specList.size() < numSpecsPerPeptide)
					specList.add(spec);
			}
		}
		for(ArrayList<Spectrum> specList : pepSpecMap.values())
			for(Spectrum spec : specList)
				specContOnePerPep.add(spec);
		gen = new NewScoringParameterGenerator(specContOnePerPep);
		
		
		// set up the tolerance
		gen.tolerance(mme);
		
		// Step 1: partition spectra
		gen.partition(NUM_SEGMENTS_PER_SPECTRUM);
		if(verbose)
			System.out.println("Partition: " + gen.partitionSet.size());
		
		// Step 2: compute offset frequency functions of precursor peaks and their neutral losses
		gen.precursorOFF();
		if(verbose)
			System.out.println("PrecursorOFF Done.");

		// Step 3: filter out "significant" precursor offsets
		gen.filterPrecursorPeaks();
		if(verbose)
			System.out.println("Filtering Done.");
		
		// Step 4: compute offset frequency fnction of fragment peaks and determine ion types to be considered for scoring
		gen.selectIonTypes();
		if(verbose)
			System.out.println("Ion types selected.");
		
		// Step 5: compute rank distributions
		gen.generateRankDist(MAX_RANK);
		if(verbose)
			System.out.println("Rank distribution computed.");

		// Step 7: smoothing parameters
		gen.smoothing();
		if(verbose)
			System.out.println("Smoothing complete.");
		
		// output
		if(!isText)
			gen.writeParameters(outputFile);
		else
			gen.writeParametersPlainText(outputFile);

		if(verbose)
			System.out.println("Writing Done.");
	}
	
	// Required
	private SpectraContainer specContainer;
	
	public NewScoringParameterGenerator(SpectraContainer specContainer)
	{
		this.specContainer = specContainer;
	}
	
	public void partition(int numSegments)
	{
		super.numSegments = numSegments;
		chargeHist = new Histogram<Integer>();
		partitionSet = new TreeSet<Partition>();
		
		Hashtable<Integer,ArrayList<Float>> parentMassMap = new Hashtable<Integer,ArrayList<Float>>();
		for(Spectrum spec : specContainer)
		{
			int charge = spec.getCharge();
			if(charge <= 0)
				continue;
			chargeHist.add(charge);
			if(spec.getAnnotation() != null)
			{
				ArrayList<Float> precursorList = parentMassMap.get(charge);
				if(precursorList == null)
				{
					precursorList = new ArrayList<Float>();
					parentMassMap.put(charge, precursorList);
				}
				precursorList.add(spec.getParentMass());
			}
		}
		
		for(int c=chargeHist.minKey(); c<=chargeHist.maxKey(); c++)
		{
			ArrayList<Float> parentMassList = parentMassMap.get(c);
			if(parentMassList == null)
				continue;
			
			int numSpec = parentMassList.size();
			if(numSpec < Math.round(MIN_NUM_SPECTRA_PER_PARTITION*0.9f))	// to few spectra
				continue;
			
			Collections.sort(parentMassList);
			int bestSetSize = 0;
			int smallestRemainder = MIN_NUM_SPECTRA_PER_PARTITION;
			for(int i=Math.round(MIN_NUM_SPECTRA_PER_PARTITION*0.9f); i<=Math.round(MIN_NUM_SPECTRA_PER_PARTITION*1.1f); i++)
			{
				int remainder = numSpec % i;
				if(i-remainder < remainder)
					remainder = i-remainder;
				if(remainder < smallestRemainder || (remainder==smallestRemainder && Math.abs(MIN_NUM_SPECTRA_PER_PARTITION-i) < Math.abs(MIN_NUM_SPECTRA_PER_PARTITION-bestSetSize)))
				{
					bestSetSize = i;
					smallestRemainder = remainder;
				}
			}
			int num=0;
			for(int i=0; i==0 || i<Math.round(numSpec/(float)bestSetSize); i++)
			{
				if(num != 0)
				{
					for(int seg=0; seg<numSegments; seg++)
						partitionSet.add(new Partition(c, parentMassList.get(num), seg));
				}
				else
				{
					for(int seg=0; seg<numSegments; seg++)
						partitionSet.add(new Partition(c, 0f, seg));
				}
				num += bestSetSize;
			}
		}
	}	
	
	private void precursorOFF()
	{
		float granularity = 0.001f;
		if(chargeHist == null)
		{
			assert(false): "partition() must have been called before";
			return;
		}
		precursorOFFMap = new TreeMap<Integer,ArrayList<PrecursorOffsetFrequency>>();
		numPrecurOFF = 0;
		
		for(int charge=chargeHist.minKey(); charge<=chargeHist.maxKey(); charge++)
		{
			if(chargeHist.get(charge) < MIN_NUM_SPECTRA_FOR_PRECURSOR_OFF)
				continue;
			ArrayList<PrecursorOffsetFrequency> precursorOffsetList = new ArrayList<PrecursorOffsetFrequency>();
			int numSpecs = 0;
			Hashtable<Integer, Histogram<Float>> histList = new Hashtable<Integer, Histogram<Float>>();
			for(int c=charge; c>=2; c--)
				histList.put(c, new Histogram<Float>());
			
	 		for(Spectrum spec : specContainer)
			{
				if(spec.getAnnotation() == null)
					continue;
				if(spec.getCharge() != charge)
					continue;
				numSpecs++;
				spec = filter.apply(spec);
				float precursorNeutralMass = spec.getParentMass();
				for(int c=charge; c>=2; c--)
				{
					float precursorMz = (precursorNeutralMass+c*(float)Composition.NEUTRON)/c;
					ArrayList<Peak> peakList = spec.getPeakListByMassRange(
							precursorMz+MIN_PRECURSOR_OFFSET/(float)c-mme.getToleranceAsDa(precursorMz+MIN_PRECURSOR_OFFSET/(float)c)/2, 
							precursorMz+MAX_PRECURSOR_OFFSET/(float)c+mme.getToleranceAsDa(precursorMz+MAX_PRECURSOR_OFFSET/(float)c)/2);

					float prevMassDiff = Float.MIN_VALUE;
					for(Peak p : peakList)
					{
						float peakMass = p.getMz();
						float massDiff = Math.round((peakMass-precursorMz)/granularity)*granularity;
						if(massDiff > prevMassDiff)
						{
							histList.get(c).add(massDiff);
							prevMassDiff = massDiff;
						}
					}
				}
			}
	 		
			for(int c=charge; c>=2; c--)
			{
				ArrayList<Float> keyList = new ArrayList<Float>(histList.get(c).keySet());
				ArrayList<Float> probList = new ArrayList<Float>();
				
				Collections.sort(keyList);
				for(Float key : keyList)
				{
					float prob = (histList.get(c).get(key))/(float)numSpecs;
					probList.add(prob);
				}
				
				float minProbThreshold = new ListStat(probList).getSignalThreshold(SIGNAL_TO_NOISE_RATIO_PRECURSOR_FILTER);
				int index = -1;
				for(Float key : keyList)
				{
					float prob = probList.get(++index);
					if(prob > minProbThreshold)
					{
						precursorOffsetList.add(new PrecursorOffsetFrequency((charge-c), key, prob));
					}
				}
			}
			ArrayList<PrecursorOffsetFrequency> clusteredOFF = PrecursorOffsetFrequency.getClusteredOFF(precursorOffsetList, granularity);
			precursorOFFMap.put(charge, clusteredOFF);
			numPrecurOFF += clusteredOFF.size();
		}
	}
	
	private void filterPrecursorPeaks()
	{
		if(this.precursorOFFMap == null)
			return;
		for(Spectrum spec : specContainer)
		{
			for(PrecursorOffsetFrequency off : this.getPrecursorOFF(spec.getCharge()))
				spec.filterPrecursorPeaks(off.getTolerance(), off.getReducedCharge(), off.getOffset());
		}
	}
	
	private Pair<Float,Float> getParentMassRange(Partition partition)
	{
		float minParentMass = partition.getParentMass();
		float maxParentMass = Float.MAX_VALUE;
		Partition higherPartition = partitionSet.higher(partition);
		if(higherPartition != null)
		{
			if(higherPartition.getCharge() == partition.getCharge() && higherPartition.getSegNum() == partition.getSegNum())
			{
				maxParentMass = higherPartition.getParentMass();
			}
		}
		return new Pair<Float,Float>(minParentMass, maxParentMass);
	}
	
	private void selectIonTypes()
	{
		if(partitionSet == null)
		{
			assert(false) : "partition() must have been called before!";
			return;
		}

		fragOFFTable = new Hashtable<Partition,ArrayList<FragmentOffsetFrequency>>();
		insignificantFragOFFTable = new Hashtable<Partition,ArrayList<FragmentOffsetFrequency>>();
		
		for(Partition partition : partitionSet)
		{
			int charge = partition.getCharge();
			// parent mass range check
			Pair<Float,Float> parentMassRange = getParentMassRange(partition);

			SpectraContainer curPartContainer = new SpectraContainer();
			for(Spectrum spec : specContainer)
			{
				if(spec.getAnnotation() == null)
					continue;
				if(spec.getCharge() != charge)
					continue;
				
				float curParentMass = spec.getParentMass();
				if(curParentMass < parentMassRange.getFirst() || curParentMass >= parentMassRange.getSecond())
					continue;
				
				curPartContainer.add(spec);
			}
			
			ArrayList<FragmentOffsetFrequency> signalFragmentOffsetFrequencyList = new ArrayList<FragmentOffsetFrequency>();
			ArrayList<FragmentOffsetFrequency> insignificantFragmentOffsetFrequencyList = new ArrayList<FragmentOffsetFrequency>();
			
			int seg = partition.getSegNum();
			IonType[] allIonTypes = IonType.getAllKnownIonTypes(Math.min(charge, 4), true).toArray(new IonType[0]);
			IonProbability probGen = new IonProbability(
					curPartContainer.iterator(), 
					allIonTypes,
					mme)
					.filter(filter)
					.segment(seg, numSegments);
	 		
			float[] ionProb = probGen.getIonProb();
			ArrayList<Float> ionProbWithoutZero = new ArrayList<Float>();
			for(float prob : ionProb)
				if(prob > 0)
					ionProbWithoutZero.add(prob);
			
			ListStat stat = new ListStat(ionProbWithoutZero);
			float signalThreshold = stat.getSignalThreshold(SIGNAL_TO_NOISE_RATIO);
			float noiseThreshold = stat.getSignalThreshold(NOISE_LEVEL_RATIO);
			for(int i=0; i<allIonTypes.length; i++)
			{
				if(ionProb[i] >= signalThreshold)
					signalFragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(allIonTypes[i], ionProb[i]));
				else if(ionProb[i] > noiseThreshold/NOISE_LEVEL_RATIO/NOISE_LEVEL_RATIO && ionProb[i] < noiseThreshold)
					insignificantFragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(allIonTypes[i], ionProb[i]));
			}
			
	 		if(signalFragmentOffsetFrequencyList.size() == 0)
	 		{
				int maxIndex = -1;
				float maxIonProb = Float.MIN_VALUE;
				for(int i=0; i<allIonTypes.length; i++)
				{
					if(ionProb[i] > maxIonProb)
					{
						maxIndex = i;
						maxIonProb = ionProb[i];
					}
				}
 				signalFragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(allIonTypes[maxIndex], maxIonProb));
	 		}
	 		Collections.sort(insignificantFragmentOffsetFrequencyList);
	 		
	 		ArrayList<FragmentOffsetFrequency> noiseFragmentOffsetFrequencyList = new ArrayList<FragmentOffsetFrequency>();
	 		for(int i=0; i<NUM_NOISE_IONS && i < insignificantFragmentOffsetFrequencyList.size(); i++)
	 			noiseFragmentOffsetFrequencyList.add(insignificantFragmentOffsetFrequencyList.get(i));
	 		fragOFFTable.put(partition, signalFragmentOffsetFrequencyList);
	 		insignificantFragOFFTable.put(partition, noiseFragmentOffsetFrequencyList);
//	 		System.out.println("**************Partition " + partition.getCharge()+","+partition.getParentMass()+","+partition.getSegNum());
//	 		System.out.println("Signal:"+"\t"+signalThreshold);
//	 		for(FragmentOffsetFrequency off : signalFragmentOffsetFrequencyList)
//	 			System.out.println(off.getIonType().getName()+"\t"+off.getFrequency());
//	 		System.out.println("Noise:"+"\t"+noiseThreshold/NOISE_LEVEL_RATIO);
//	 		for(FragmentOffsetFrequency off : noiseFragmentOffsetFrequencyList)
//	 			System.out.println(off.getIonType().getName()+"\t"+off.getFrequency());
		}
	}	
	
	private void generateRankDist(int maxRank)
	{
		if(partitionSet == null)
		{
			assert(false): "partition() must have been called!";
			return;
		}
		
		rankDistTable = new Hashtable<Partition,Hashtable<IonType,Float[]>>();
		this.maxRank = maxRank;
		
		for(Partition partition : partitionSet)
		{
			int charge = partition.getCharge();
			IonType[] ionTypes = getIonTypes(partition);
			IonType[] noiseIons = getNoiseIonTypes(partition);
			IonType[] allIons = Arrays.copyOf(ionTypes, ionTypes.length+noiseIons.length);
			for(int i=0; i<noiseIons.length; i++)
				allIons[ionTypes.length+i] = noiseIons[i];
			
			Pair<Float,Float> parentMassRange = getParentMassRange(partition);
			int seg = partition.getSegNum();
			
			int numSpec = 0;
			Hashtable<IonType, Histogram<Integer>> rankDist = new Hashtable<IonType, Histogram<Integer>>();
			
			for(IonType ion : allIons)
				rankDist.put(ion, new Histogram<Integer>());
			rankDist.put(IonType.NOISE, new Histogram<Integer>());
			
			int numMaxRankPeaks = 0;
			int totalCleavageSites = 0;
			
			for(Spectrum spec : specContainer)
			{
				int numExplainedPeaks = 0;
				if(spec.getAnnotation() == null)
					continue;
				if(spec.getCharge() != charge)
					continue;
				float curParentMass = spec.getParentMass();
				if(curParentMass < parentMassRange.getFirst() || curParentMass >= parentMassRange.getSecond())
					continue;
				
				Peptide annotation = spec.getAnnotation();
				spec.setRanksOfPeaks();
				numSpec++;
				numMaxRankPeaks += spec.size()-maxRank+1;
				totalCleavageSites += annotation.size()-1;
				float prm = 0;
				float srm = 0;
				
				int numSignalBinsAtThisSegment = 0;
				for(int i=0; i<annotation.size()-1; i++)
				{
					prm += annotation.get(i).getMass();
					srm += annotation.get(annotation.size()-1-i).getMass();
					
					for(IonType ion : allIons)
					{
						float theoMass;
						if(ion instanceof IonType.PrefixIon)
							theoMass = ion.getMz(prm);
						else
							theoMass = ion.getMz(srm);
						
						int segNum = super.getSegmentNum(theoMass, curParentMass);
						if(segNum == seg)
						{
							numSignalBinsAtThisSegment++;
							Peak p = spec.getPeakByMass(theoMass, mme);
							if(p != null)
							{
								numExplainedPeaks++;
								int rank = p.getRank();
								if(rank >= maxRank)
								{
									rank = maxRank;
								}
								rankDist.get(ion).add(rank);
							}
							else
							{
								rankDist.get(ion).add(maxRank+1);	// maxRank+1: missing ion
							}							
						}
					}
				}
			}
			
			Hashtable<IonType,Float[]> freqDist = new Hashtable<IonType,Float[]>();
			for(IonType ion : ionTypes)
			{
				Float[] dist = new Float[maxRank+1];
				Histogram<Integer> hist = rankDist.get(ion);
				for(int i=1; i<=maxRank+1; i++)
				{
					Integer num = hist.get(i);
					dist[i-1] = (num/(float)numSpec);
				}
				freqDist.put(ion, dist);
			}
			
			// noise
			Float[] noiseDist = new Float[maxRank+1];
			for(int i=0; i<noiseDist.length; i++)
				noiseDist[i] = 0f;
			for(IonType ion : noiseIons)
			{
				Histogram<Integer> hist = rankDist.get(ion);
				for(int i=1; i<=maxRank+1; i++)
				{
					Integer num = hist.get(i);
					noiseDist[i-1] += num;
				}
			}
			for(int i=0; i<noiseDist.length; i++)
				noiseDist[i] = noiseDist[i]/numSpec/noiseIons.length;
			freqDist.put(IonType.NOISE, noiseDist);
			
			rankDistTable.put(partition, freqDist);
		}
	}
	
	protected void smoothing()
	{
		smoothingRankDistTable();
	}
	
	protected void smoothingRankDistTable()
	{
		if(rankDistTable == null)
			return;
		assert(smoothingRanks.length == smoothingWindowSize.length);
		for(Partition partition : rankDistTable.keySet())
		{
			Hashtable<IonType,Float[]> table = this.rankDistTable.get(partition);
			for(IonType ion : table.keySet())
			{
				Float[] freq = table.get(ion);
				Float[] smoothedFreq = new Float[freq.length];
				int smoothingIndex = 0;
				for(int i=0; i<freq.length-2; i++)	// last 2 columns: maxRank, unexplained
				{
					if(smoothingIndex < smoothingRanks.length-1 &&
							i == smoothingRanks[smoothingIndex])
						smoothingIndex++;
					int windowSize = smoothingWindowSize[smoothingIndex];
					float sumFrequencies = 0;
					int numIndicesSummed = 0;
					for(int d=-windowSize; d<=windowSize; d++)
					{
						int index = i+d;
						if(index < 0 || index > freq.length-3)
							continue;
						sumFrequencies += freq[index];
						numIndicesSummed++;
					}
					while(sumFrequencies == 0 && windowSize < freq.length-4)
					{
						windowSize++;
						int index = i-windowSize;
						if(index >= 0)
						{
							sumFrequencies += freq[index];
							numIndicesSummed++;
						}
						index = i+windowSize;
						if(index <= freq.length-3)
						{
							sumFrequencies += freq[index];
							numIndicesSummed++;
						}
					}
					if(sumFrequencies != 0)
						smoothedFreq[i] = sumFrequencies/numIndicesSummed;
					else
						assert(false);
				}
				for(int i=0; i<freq.length-2; i++)
					freq[i] = smoothedFreq[i];
				if(freq[freq.length-1] == 0)
					freq[freq.length-1] = Float.MIN_VALUE;
				if(freq[freq.length-2] == 0)	
					freq[freq.length-2] = freq[freq.length-3];
			}
		}
	}
}

