package msscorer;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.TreeMap;
import java.util.TreeSet;

import msgf.Histogram;
import msgf.IntHistogram;
import msgf.NominalMass;
import msgf.Tolerance;
import msscorer.NewScorerFactory.SpecDataType;
import msutil.ActivationMethod;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Enzyme;
import msutil.InstrumentType;
import msutil.IonType;
import msutil.Pair;
import msutil.Peak;
import msutil.Peptide;
import msutil.Protocol;
import msutil.SpectraContainer;
import msutil.Spectrum;
import msutil.IonType.PrefixIon;
import parser.MgfSpectrumParser;

/**
 * This only supports low accuracy fragment ions.
 * @author sangtaekim
 *
 */
public class ScoringParameterGeneratorWithErrors extends NewRankScorer {
	private static final float MIN_PRECURSOR_OFFSET = -300f;	// for precursors
	private static final float MAX_PRECURSOR_OFFSET = 30f;
	private static final int MIN_NUM_SPECTRA_PER_PARTITION = 400;	// 400
	private static final int MIN_NUM_SPECTRA_FOR_PRECURSOR_OFF = 150;

	private static final float MIN_PRECURSOR_OFFSET_PROBABILITY = 0.15f;	// 0.15
	private static final float MIN_ION_OFFSET_PROBABILITY = 0.15f;	// 0.15, for ion types
	private static final float MIN_MAIN_ION_OFFSET_PROBABILITY = 0.01f;	// ion of probability below this number will be ignored 

	private static final int MAX_RANK = 150;
	private static final int NUM_SEGMENTS_PER_SPECTRUM = 2;	// 2

	private static final int[] smoothingRanks = {3, 5, 10, 20, 50, Integer.MAX_VALUE}; //Ranks around which smoothing occurs
	private static final int[] smoothingWindowSize = {0, 1, 2, 3, 4, 5}; //Smoothing windows for each smoothing rank

	private static final float DECONVOLUTION_MASS_TOLERANCE = 0.02f;
	protected static final int MAX_CHARGE = 20;

	public static void generateParameters(
			File specFile,
			SpecDataType dataType,
			AminoAcidSet aaSet, 
			File outputDir,
			boolean isText, 
			boolean verbose,
			boolean singlePartition
			)
	{
		SpectraContainer container = new SpectraContainer(specFile.getPath(), new MgfSpectrumParser().aaSet(aaSet));
		generateParameters(container, dataType, aaSet, outputDir, isText, verbose, singlePartition);
	}

	public static void generateParameters(
			SpectraContainer container, 
			SpecDataType dataType,
			AminoAcidSet aaSet,
			File outputDir,
			boolean isText, 
			boolean verbose)
	{
		generateParameters(container, dataType, aaSet, outputDir, isText, verbose, false);
	}
	
	public static void generateParameters(
			SpectraContainer container, 
			SpecDataType dataType,
			AminoAcidSet aaSet,
			File outputDir,
			boolean isText, 
			boolean verbose,
			boolean singlePartition)
	{
		if(verbose)
			System.out.println("Number of annotated PSMs: " + container.size());
		
		String paramFileName = dataType.toString()+".param";
		
		File outputFile;
		if(outputDir != null)
			outputFile = new File(outputDir, paramFileName);
		else
			outputFile = new File(paramFileName);
		
		if(verbose)
			System.out.println("Output file name: " + outputFile.getAbsolutePath());
		int errorScalingFactor = 0;
		boolean applyDeconvolution = false;
		
		if(dataType.getInstrumentType() == InstrumentType.HIGH_RESOLUTION_LTQ 
				|| dataType.getInstrumentType() == InstrumentType.TOF)
		{
			errorScalingFactor = 100;
			applyDeconvolution = true;
			if(verbose)
				System.out.println("High-precision MS/MS data: " +
						"errorScalingFactor("+errorScalingFactor+") " +
						"chargeDeconvolution("+applyDeconvolution+")");
		}
		
		boolean considerPhosLoss = false;
		if(dataType.getProtocol().getName().equals("Phosphorylation"))
		{
			considerPhosLoss = true;
			if(verbose)
				System.out.println("Consider H3PO4 loss.");
		}
		
		HashSet<String> pepSet = new HashSet<String>();
		for(Spectrum spec : container)
			pepSet.add(spec.getAnnotationStr());

		if(verbose)
			System.out.println("Number of unique peptides: " + pepSet.size());
		int numSpecsPerPeptide;
		if(pepSet.size() < 2000)
		{
			numSpecsPerPeptide = 3;
		}
		else
		{
			numSpecsPerPeptide = 1;
		}
		if(verbose)
			System.out.println("Consider " + numSpecsPerPeptide + " per spectrum.");
		
		// multiple spectra with the same peptide -> one spec per peptide
		HashMap<String,ArrayList<Spectrum>> pepSpecMap = new HashMap<String,ArrayList<Spectrum>>();
		for(Spectrum spec : container)
		{
			if(spec.getAnnotationStr() == null)
				continue;
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

		SpectraContainer specContOnePerPep = new SpectraContainer();
		for(ArrayList<Spectrum> specList : pepSpecMap.values())
		{
			for(Spectrum spec : specList)
			{
				specContOnePerPep.add(spec);
			}
		}

		ScoringParameterGeneratorWithErrors gen = new ScoringParameterGeneratorWithErrors(
				specContOnePerPep, 
				dataType,
				considerPhosLoss, 
				applyDeconvolution);
		
		// set up the tolerance
		gen.tolerance(new Tolerance(0.5f));

		// Step 1: partition spectra
		if(singlePartition)
			gen.partition(2, true);
		else
			gen.partition(NUM_SEGMENTS_PER_SPECTRUM, false);
		if(verbose)
			System.out.println("Partition: " + gen.partitionSet.size());

		// Step 2: compute offset frequency functions of precursor peaks and their neutral losses
		gen.precursorOFF(MIN_PRECURSOR_OFFSET_PROBABILITY);
		if(verbose)
			System.out.println("PrecursorOFF Done.");

		// Step 3: filter out "significant" precursor offsets
		gen.filterPrecursorPeaks();
		if(verbose)
			System.out.println("Filtering Done.");

		if(applyDeconvolution)
		{
			gen.deconvoluteSpectra();
			if(verbose)
				System.out.println("Deconvolution Done.");
		}

		// Step 4: compute offset frequency function of fragment peaks and determine ion types to be considered for scoring
		gen.selectIonTypes();
		if(verbose)
			System.out.println("Ion types selected.");

		// Step 5: compute rank distributions
		gen.generateRankDist(MAX_RANK);
		if(verbose)
			System.out.println("Rank distribution computed.");

		// Step 6 (optional): generate error distribution
		gen.generateErrorDist(errorScalingFactor);
		if(verbose)
			System.out.println("Error disbribution computed");

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
	private final boolean considerPhosLoss;

	public ScoringParameterGeneratorWithErrors(SpectraContainer specContainer, SpecDataType dataType, boolean considerPhosLoss, boolean applyDeconvolution)
	{
		this.specContainer = specContainer;
		this.considerPhosLoss = considerPhosLoss;
		super.dataType = dataType;
		super.applyDeconvolution = applyDeconvolution;
		super.deconvolutionErrorTolerance = DECONVOLUTION_MASS_TOLERANCE;
	}

	public void partition(int numSegments, boolean singlePartition)
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
			
			if(singlePartition)
				bestSetSize = numSpec;
			else
			{
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

	private void precursorOFF(float minProbThreshold)
	{
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
			Hashtable<Integer, Histogram<Integer>> histList = new Hashtable<Integer, Histogram<Integer>>();
			for(int c=charge; c>=2; c--)
				histList.put(c, new Histogram<Integer>());

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
					float precursorMz = (precursorNeutralMass+c*(float)Composition.H)/c;
					ArrayList<Peak> peakList = spec.getPeakListByMassRange(
							precursorMz+MIN_PRECURSOR_OFFSET/(float)c-mme.getToleranceAsDa(precursorMz+MIN_PRECURSOR_OFFSET/(float)c)/2, 
							precursorMz+MAX_PRECURSOR_OFFSET/(float)c+mme.getToleranceAsDa(precursorMz+MAX_PRECURSOR_OFFSET/(float)c)/2);

					int prevMassIndexDiff = Integer.MIN_VALUE;
					for(Peak p : peakList)
					{
						float peakMass = p.getMz();
						int massIndexDiff = NominalMass.toNominalMass(peakMass - precursorMz);
						if(massIndexDiff > prevMassIndexDiff)
						{
							histList.get(c).add(massIndexDiff);
							prevMassIndexDiff = massIndexDiff;
						}
					}
				}
			}

			for(int c=charge; c>=2; c--)
			{
				ArrayList<Integer> keyList = new ArrayList<Integer>(histList.get(c).keySet());
				Collections.sort(keyList);
				for(Integer key : keyList)
				{
					float prob = (histList.get(c).get(key))/(float)numSpecs;
					if(prob > minProbThreshold)
					{
						precursorOffsetList.add(new PrecursorOffsetFrequency((charge-c), NominalMass.getMassFromNominalMass(key), prob));
					}
				}
			}
			precursorOFFMap.put(charge, precursorOffsetList);
			numPrecurOFF += precursorOffsetList.size();
		}
	}

	private void filterPrecursorPeaks()
	{
		if(this.precursorOFFMap == null)
			return;
		for(Spectrum spec : specContainer)
		{
			for(PrecursorOffsetFrequency off : this.getPrecursorOFF(spec.getCharge()))
				spec.filterPrecursorPeaks(mme, off.getReducedCharge(), off.getOffset());
		}
	}

	private void deconvoluteSpectra()
	{
		SpectraContainer newSpecContainer = new SpectraContainer();
		for(Spectrum spec : specContainer)
		{
			newSpecContainer.add(spec.getDeconvolutedSpectrum(DECONVOLUTION_MASS_TOLERANCE));
		}
		specContainer = newSpecContainer;
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

			int seg = partition.getSegNum();
			IonType[] allIonTypes = IonType.getAllKnownIonTypes(Math.min(charge, 4), true, considerPhosLoss).toArray(new IonType[0]);
			IonProbability probGen = new IonProbability(
					curPartContainer.iterator(), 
					allIonTypes,
					mme)
			.filter(filter)
			.segment(seg, numSegments);

			float[] ionProb = probGen.getIonProb();

			float signalThreshold = MIN_ION_OFFSET_PROBABILITY;
			for(int i=0; i<allIonTypes.length; i++)
			{
				if(ionProb[i] >= signalThreshold)
					signalFragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(allIonTypes[i], ionProb[i]));
			}

			if(signalFragmentOffsetFrequencyList.size() == 0)
			{
				int maxIndex = -1;
				float maxIonProb = Float.MIN_VALUE;
				for(int i=0; i<allIonTypes.length; i++)
				{
					if(ionProb[i] > MIN_MAIN_ION_OFFSET_PROBABILITY && ionProb[i] > maxIonProb)
					{
						maxIndex = i;
						maxIonProb = ionProb[i];
					}
				}
				if(maxIndex >= 0)
					signalFragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(allIonTypes[maxIndex], maxIonProb));
			}

			Collections.sort(signalFragmentOffsetFrequencyList, Collections.reverseOrder());
			fragOFFTable.put(partition, signalFragmentOffsetFrequencyList);
		}
		super.determineIonTypes();
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
			if(ionTypes == null || ionTypes.length == 0)
				continue;
			Pair<Float,Float> parentMassRange = getParentMassRange(partition);
			int seg = partition.getSegNum();

			int numSpec = 0;
			Hashtable<IonType, Histogram<Integer>> rankDist = new Hashtable<IonType, Histogram<Integer>>();
			Hashtable<IonType, Float> rankDistMaxRank = new Hashtable<IonType, Float>();
			Hashtable<IonType, Float> rankDistUnexplained = new Hashtable<IonType, Float>();

			for(IonType ion : ionTypes)
			{
				rankDist.put(ion, new Histogram<Integer>());
				rankDistMaxRank.put(ion, 0f);
				rankDistUnexplained.put(ion, 0f);
			}
			rankDist.put(IonType.NOISE, new Histogram<Integer>());

			float[] noiseDist = new float[maxRank+2];
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
				int prmMassIndex = 0;
				int srmMassIndex = 0;

				HashSet<Peak> explainedPeakSet = new HashSet<Peak>();
				Hashtable<IonType, Integer> numExplainedMaxRankPeaks = new Hashtable<IonType, Integer>();
				for(IonType ion : ionTypes)
				{
					numExplainedMaxRankPeaks.put(ion, 0);
				}

				int numSignalBinsAtThisSegment = 0;
				for(int i=0; i<annotation.size()-1; i++)
				{
					prmMassIndex += NominalMass.toNominalMass(annotation.get(i).getMass());
					srmMassIndex += NominalMass.toNominalMass(annotation.get(annotation.size()-1-i).getMass());

					float prm = NominalMass.getMassFromNominalMass(prmMassIndex);
					float srm = NominalMass.getMassFromNominalMass(srmMassIndex);
					for(IonType ion : ionTypes)
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
									numExplainedMaxRankPeaks.put(ion, numExplainedMaxRankPeaks.get(ion)+1);
								}
								explainedPeakSet.add(p);
								rankDist.get(ion).add(rank);
							}
							else
							{
								rankDist.get(ion).add(maxRank+1);	// maxRank+1: missing ion
							}							
						}
					}
				}

				ArrayList<Peak> unexplainedPeaksAtThisSegment = new ArrayList<Peak>();
				int numPeaksAtThisSegment = 0;
				int numMaxRankPeaksAtThisSegment = 0;
				for(Peak p : spec)
				{
					if(super.getSegmentNum(p.getMz(), curParentMass) == seg)
					{
						numPeaksAtThisSegment++;
						if(p.getRank() >= maxRank)
							numMaxRankPeaksAtThisSegment++;
						if(!explainedPeakSet.contains(p))
							unexplainedPeaksAtThisSegment.add(p);
					}
				}

				float midMassThisSegment = (1f/numSegments*seg+1f/numSegments/2)*annotation.getParentMass();
				float numBinsAtThisSegment = annotation.getParentMass()/numSegments/mme.getToleranceAsDa(midMassThisSegment)/2;

				for(Peak p : unexplainedPeaksAtThisSegment)
				{
					int rank = p.getRank();
					//					float noiseFreq = (float)(annotation.size()-1)/(annotation.getParentMass()/(mme.getToleranceAsDa(midMassThisSegment)*2));
					float noiseFreq = (annotation.size()-1)/numSegments/numBinsAtThisSegment;
					if(rank >= maxRank)
						noiseDist[maxRank] += noiseFreq/numMaxRankPeaksAtThisSegment;
					else
						noiseDist[rank] += noiseFreq;
				}

				for(IonType ion : ionTypes)
				{
					if(numMaxRankPeaksAtThisSegment > 0)
					{
						Float prevSumFreq = rankDistMaxRank.get(ion);
						float curFreq = numExplainedMaxRankPeaks.get(ion)/(float)numMaxRankPeaksAtThisSegment;
						rankDistMaxRank.put(ion, prevSumFreq+curFreq);
					}
				}

				noiseDist[maxRank+1] += (numBinsAtThisSegment-numPeaksAtThisSegment)*(annotation.size()-1)/numSegments/numBinsAtThisSegment;
			}

			Hashtable<IonType,Float[]> freqDist = new Hashtable<IonType,Float[]>();
			for(IonType ion : ionTypes)
			{
				Float[] dist = new Float[maxRank+1];
				Histogram<Integer> hist = rankDist.get(ion);
				for(int i=1; i<=maxRank-1; i++)
				{
					Integer num = hist.get(i);
					dist[i-1] = (num/(float)numSpec);
				}
				dist[maxRank-1] = rankDistMaxRank.get(ion)/numSpec;
				dist[maxRank] = hist.get(maxRank+1)/(float)numSpec;
				freqDist.put(ion, dist);
			}

			// noise
			Float[] dist = new Float[maxRank+1];
			for(int i=1; i<=maxRank+1; i++)
				dist[i-1] = noiseDist[i]/numSpec;
			freqDist.put(IonType.NOISE, dist);

			rankDistTable.put(partition, freqDist);
		}
	}

	private void generateErrorDist(int errorScalingFactor)
	{
		this.errorScalingFactor = errorScalingFactor;
		if(errorScalingFactor > 0)
		{
			generateIonErrorDist();
			generateNoiseErrorDist();
		}
	}

	private void generateIonErrorDist()
	{
		ionErrDistTable = new Hashtable<Partition,Float[]>();
		ionExistenceTable = new Hashtable<Partition,Float[]>();
		for(Partition partition : partitionSet)
		{
			int charge = partition.getCharge();
			Pair<Float,Float> parentMassRange = getParentMassRange(partition);
			int seg = partition.getSegNum();
			if(seg != super.getNumSegments()-1)
				continue;
			IonType mainIon = this.getMainIonType(partition);
			IntHistogram errHist = new IntHistogram();
			int[] edgeCount = new int[4];
			int numSpecs = 0;
			for(Spectrum spec : specContainer)
			{
				if(spec.getAnnotation() == null)
					continue;
				if(spec.getCharge() != charge)
					continue;

				float curParentMass = spec.getParentMass();
				if(curParentMass < parentMassRange.getFirst() || curParentMass >= parentMassRange.getSecond())
					continue;

				numSpecs++;
				Peptide peptide;

				peptide = spec.getAnnotation();

				int intResidueMass = 0;
				float[] obsMass = new float[peptide.size()+1];

				obsMass[0] = mainIon.getOffset(); 
				for(int i=0; i<peptide.size()-1; i++)
				{
					if(mainIon instanceof PrefixIon)
						intResidueMass += peptide.get(i).getNominalMass();
					else
						intResidueMass += peptide.get(peptide.size()-1-i).getNominalMass();

					float theoMass = mainIon.getMz(NominalMass.getMassFromNominalMass(intResidueMass));
					Peak p = spec.getPeakByMass(theoMass, mme);
					if(p != null)
						obsMass[i+1] = p.getMz();
					else
						obsMass[i+1] = -1;
				}

				obsMass[peptide.size()] = mainIon.getMz(peptide.getMass());
				for(int i=1; i<=peptide.size(); i++)
				{
					if(obsMass[i] >= 0)
					{
						if(obsMass[i-1] >=0)		// yy
						{
							AminoAcid aa;
							if(mainIon instanceof PrefixIon)
								aa = peptide.get(i-1);
							else
								aa = peptide.get(peptide.size()-i);

							float expMass = obsMass[i]-obsMass[i-1];
							float theoMass = aa.getMass()/mainIon.getCharge();
							float diff = expMass-theoMass;
							int diffIndex = Math.round(diff*errorScalingFactor);
							if(diffIndex > errorScalingFactor)
								diffIndex = errorScalingFactor;
							else if(diffIndex < -errorScalingFactor)
								diffIndex = -errorScalingFactor;
							errHist.add(diffIndex);
							edgeCount[3]++;
						}
						else	// ny
							edgeCount[1]++;
					}
					else
					{
						if(obsMass[i-1] >=0)		// yn
							edgeCount[2]++;
						else						// nn
							edgeCount[0]++;
					}
				}	
			}

			Float[] ionErrHist = new Float[2*errorScalingFactor+1];
			// smoothing
			float[] smoothedHist = errHist.getSmoothedHist(errorScalingFactor);
			for(int i=-errorScalingFactor; i<=errorScalingFactor; i++)
				ionErrHist[i+errorScalingFactor] = smoothedHist[i+errorScalingFactor]/(float)errHist.totalCount();

			Float[] ionExistence = new Float[edgeCount.length];
			int sumEdgeCount = 0;
			for(int i=0; i<edgeCount.length; i++)
				sumEdgeCount += edgeCount[i];
			for(int i=0; i<edgeCount.length; i++)
				ionExistence[i] = edgeCount[i]/(float)sumEdgeCount;

			for(int i=0; i<this.numSegments; i++)
			{
				Partition part = new Partition(partition.getCharge(), partition.getParentMass(), i);
				if(partitionSet.contains(part))
				{
					ionErrDistTable.put(part, ionErrHist);
					ionExistenceTable.put(part, ionExistence);
				}
			}
			//			if(partition.getCharge() == 2 && partition.getParentMass() > 1000 && partition.getParentMass() < 1110)
			//			{
			//				System.out.println("Partition\t"+partition.getCharge()+"\t"+partition.getParentMass());
			//				System.out.println("ErrorHist:");
			//				for(int i=0; i<errorScalingFactor*2+1; i++)
			//					System.out.println((i-errorScalingFactor)+"\t"+errHist.get(i-errorScalingFactor)+"\t"+ionErrHist[i]);
			//				System.out.println("IonExistence:");
			//				for(int i=0;i<ionExistence.length; i++)
			//					System.out.println(i+"\t"+ionExistence[i]);
			//			}
		}
	}

	private void generateNoiseErrorDist()
	{
		this.noiseErrDistTable = new Hashtable<Partition,Float[]>();
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		AminoAcid aaK = aaSet.getAminoAcid('K');
		AminoAcid aaQ = aaSet.getAminoAcid('Q');
		int heaviestAANominalMass = aaSet.getMaxNominalMass();
		float[] nominalMass = new float[heaviestAANominalMass+1];
		for(AminoAcid aa : aaSet)
			nominalMass[aa.getNominalMass()] = aa.getMass();

		for(Partition partition : partitionSet)
		{
			int charge = partition.getCharge();
			Pair<Float,Float> parentMassRange = getParentMassRange(partition);
			int seg = partition.getSegNum();
			if(seg != super.getNumSegments()-1)
				continue;

			IntHistogram errHist = new IntHistogram();
			int numSpecs = 0;
			for(Spectrum spec : specContainer)
			{
				if(spec.getAnnotation() == null)
					continue;
				if(spec.getCharge() != charge)
					continue;

				float curParentMass = spec.getParentMass();
				if(curParentMass < parentMassRange.getFirst() || curParentMass >= parentMassRange.getSecond())
					continue;

				Spectrum noiseSpec = (Spectrum)spec.clone();

				numSpecs++;

				for(int i=0; i<noiseSpec.size()-1; i++)
				{
					Peak p1 = noiseSpec.get(i);
					float p1Mass = p1.getMz();
					int nominalP1 = NominalMass.toNominalMass(p1.getMz());
					for(int j=i+1; j<noiseSpec.size(); j++)
					{
						Peak p2 = noiseSpec.get(j);
						float p2Mass = p2.getMz();
						int nominalP2 = NominalMass.toNominalMass(p2.getMz());
						int nominalDiff = nominalP2-nominalP1;
						if(nominalDiff > heaviestAANominalMass)
							break;
						if(nominalMass[nominalDiff] == 0)
							continue;

						float diff = p2Mass-p1Mass;
						float aaMass = nominalMass[nominalDiff];
						if(nominalDiff == 128)	// K or Q
						{
							if(Math.abs(diff-aaQ.getMass()) > Math.abs(diff-aaK.getMass()))
								aaMass = aaK.getMass();
							else
								aaMass = aaQ.getMass();
						}
						float err = diff-aaMass;
						errHist.add(Math.round(err*errorScalingFactor));
					}
				}
			}
			Float[] noiseErrHist = new Float[2*errorScalingFactor+1];
			// smoothing
			float[] smoothedHist = errHist.getSmoothedHist(errorScalingFactor);
			for(int i=-errorScalingFactor; i<=errorScalingFactor; i++)
				noiseErrHist[i+errorScalingFactor] = smoothedHist[i+errorScalingFactor]/(float)errHist.totalCount();

			for(int i=0; i<this.numSegments; i++)
			{
				Partition part = new Partition(partition.getCharge(), partition.getParentMass(), i);
				if(partitionSet.contains(part))
				{
					noiseErrDistTable.put(part, noiseErrHist);
				}
			}
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

