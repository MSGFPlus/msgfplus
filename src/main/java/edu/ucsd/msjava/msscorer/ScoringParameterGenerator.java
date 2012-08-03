package edu.ucsd.msjava.msscorer;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.ucsd.msjava.msgf.Histogram;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Constants;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.IonType;
import edu.ucsd.msjava.msutil.Pair;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.SpectraContainer;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.WindowFilter;
import edu.ucsd.msjava.msutil.IonType.PrefixIon;
import edu.ucsd.msjava.parser.MgfSpectrumParser;

/**
 * This only supports low accuracy fragment ions.
 * @author sangtaekim
 *
 */
public class ScoringParameterGenerator extends NewRankScorer {
	private static final float MIN_OFFSET_MASS = -120f;	// for ion types
	private static final float MAX_OFFSET_MASS = 38f;	
	private static final float MIN_PRECURSOR_OFFSET = -300f;	// for precursors
	private static final float MAX_PRECURSOR_OFFSET = 30f;
	private static final int MIN_NUM_SPECTRA_PER_PARTITION = 400;	// 400
	private static final int MIN_NUM_SPECTRA_FOR_PRECURSOR_OFF = 150;
	
	private static final float MIN_PRECURSOR_OFFSET_PROBABILITY = 0.15f;	// 0.15
	private static final float MIN_ION_OFFSET_PROBABILITY = 0.15f;	// 0.15, for ion types
	private static final int MAX_RANK = 150;
	private static final int NUM_SEGMENTS_PER_SPECTRUM = 2;	// 2
	
	
    private static final int[] smoothingRanks = {3, 5, 10, 20, 50, Integer.MAX_VALUE}; //Ranks around which smoothing occurs
    private static final int[] smoothingWindowSize = {0, 1, 2, 3, 4, 5}; //Smoothing windows for each smoothing rank
	
	private static final int NUM_NOISE_IONS = 10;
	protected static final int MAX_CHARGE = 20;
 	
	public static void main(String argv[])
	{
		File specFile = null;
		File outputFile = null;
		boolean isText = false;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		int numSpecsPerPeptide = 1;
		int errorScalingFactor = 10;
		
		ActivationMethod activationMethod = null;
		InstrumentType instType = null;
		Enzyme enzyme = null;
		
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
			else if(argv[i].equalsIgnoreCase("-pep"))
			{
				numSpecsPerPeptide = Integer.parseInt(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-err"))
			{
				errorScalingFactor = Integer.parseInt(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-m"))	// Fragmentation method
			{
				// (0: written in the spectrum, 1: CID , 2: ETD, 3: HCD)
				if(argv[i+1].equalsIgnoreCase("1"))
				{
					activationMethod = ActivationMethod.CID;
				}
				else if(argv[i+1].equalsIgnoreCase("2"))
				{
					activationMethod = ActivationMethod.ETD;
				}
				else if(argv[i+1].equalsIgnoreCase("3"))
				{
					activationMethod = ActivationMethod.HCD;
				}
				else
				{
					printUsageAndExit("Illegal activation method: " + argv[i+1]);
				}
			}			
			else if(argv[i].equalsIgnoreCase("-inst"))	// Instrument type
			{
				if(argv[i+1].equalsIgnoreCase("0"))
				{
					instType = InstrumentType.LOW_RESOLUTION_LTQ;
				}
				else if(argv[i+1].equalsIgnoreCase("1"))
				{
					instType = InstrumentType.TOF;
				}
				else if(argv[i+1].equalsIgnoreCase("2"))
				{
					instType = InstrumentType.HIGH_RESOLUTION_LTQ;
				}
				else
				{
					printUsageAndExit("Illegal instrument type: " + argv[i+1]);
				}
			}		
			else if(argv[i].equalsIgnoreCase("-e"))	// Enzyme
			{
				// 0: No enzyme, 1: Trypsin, 2: Chymotrypsin, 3: LysC, 4: LysN, 5: GluC, 6: ArgC, 7: AspN
				if(argv[i+1].equalsIgnoreCase("0"))
					enzyme = null;
				else if(argv[i+1].equalsIgnoreCase("1"))
					enzyme = Enzyme.TRYPSIN;
				else if(argv[i+1].equalsIgnoreCase("2"))
					enzyme = Enzyme.CHYMOTRYPSIN;
				else if(argv[i+1].equalsIgnoreCase("3"))
					enzyme = Enzyme.LysC;
				else if(argv[i+1].equalsIgnoreCase("4"))
					enzyme = Enzyme.LysN;
				else if(argv[i+1].equalsIgnoreCase("5"))
					enzyme = Enzyme.GluC;
				else if(argv[i+1].equalsIgnoreCase("6"))
					enzyme = Enzyme.ArgC;
				else if(argv[i+1].equalsIgnoreCase("7"))
					enzyme = Enzyme.AspN;
				else
					printUsageAndExit("Illegal enzyme: " + argv[i+1]);
			}			
			else
				printUsageAndExit("Illegal parameters!");
		}
		if(specFile == null)
			printUsageAndExit("missing annotatedMgfFileName!");
		if(outputFile == null)
			printUsageAndExit("missing outputFileName!");
		if(activationMethod == null)
			printUsageAndExit("missing activationMethod!");
		if(instType == null)
			printUsageAndExit("missing instrumentType!");
		
		generateParameters(specFile, activationMethod, instType, enzyme, Protocol.NOPROTOCOL, numSpecsPerPeptide, errorScalingFactor, outputFile, aaSet, isText, false);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.err.println(message);
		System.out.println("usage: java -Xmx2000M -cp MSGF.jar msscorer.ScoringParameterGenerator\n" +
				"\t-i annotatedMgfFileName (*.mgf)\n" +
				"\t-o outputFileName (e.g. CID_Tryp.param)\n" +
				"\t-m FragmentationMethodID (1: CID, 2: ETD, 3: HCD)\n" +
				"\t-inst InstrumentID (0: Low-res LCQ/LTQ, 1: TOF , 2: High-res LTQ)\n" +
				"\t-e EnzymeID (0: No enzyme, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N)\n" +
				"\t[-fixMod 0/1/2] (0: NoCysteineProtection, 1: CarbamidomethyC (default), 2: CarboxymethylC)\n" +
				"\t[-pep numPeptidesPerSpec]  (default: 1)\n" +
				"\t[-err errorScalingFactor]  (default: 10)"
		);
		System.exit(0);
	}
	
	public static void generateParameters(
			File specFile, 
			ActivationMethod activationMethod,
			InstrumentType instType,
			Enzyme enzyme,
			Protocol protocol,
			int numSpecsPerPeptide, 
			int errorScalingFactor,
			File outputFile, 
			AminoAcidSet aaSet, 
			boolean isText, 
			boolean verbose)
	{
		SpectraContainer container = new SpectraContainer(specFile.getPath(), new MgfSpectrumParser().aaSet(aaSet));
		
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
		
		SpecDataType dataType = new SpecDataType(activationMethod, instType, enzyme, protocol);
		ScoringParameterGenerator gen = new ScoringParameterGenerator(specContOnePerPep, dataType);
		
		// set up the tolerance
		gen.tolerance(new Tolerance(1/Constants.INTEGER_MASS_SCALER/2));
		
		// Step 1: partition spectra
		gen.partition(NUM_SEGMENTS_PER_SPECTRUM);
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
		
		// Step 4: compute offset frequency fnction of fragment peaks and determine ion types to be considered for scoring
		gen.selectIonTypes(MIN_ION_OFFSET_PROBABILITY);
		if(verbose)
			System.out.println("Ion types selected.");
		
		// Step 5: compute rank distributions
		gen.generateRankDist(MAX_RANK);
		if(verbose)
			System.out.println("Rank distribution computed.");

		// Step 6 (optional): generate error distribution, currently not in use

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
	
	public ScoringParameterGenerator(SpectraContainer specContainer, SpecDataType dataType)
	{
		this.specContainer = specContainer;
		super.dataType = dataType;
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
	
	private void selectIonTypes(float minProbThreshold)
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
			int seg = partition.getSegNum();

			int numSpec = 0;
			Hashtable<Integer, Histogram<Integer>> prefixIonFreq = new Hashtable<Integer, Histogram<Integer>>();
			Hashtable<Integer, Histogram<Integer>> suffixIonFreq = new Hashtable<Integer, Histogram<Integer>>();
			for(int c=1; c<=charge; c++)
			{
				prefixIonFreq.put(c, new Histogram<Integer>());
				suffixIonFreq.put(c, new Histogram<Integer>());
			}
			
			int numCleavages = 0;
			for(Spectrum spec : specContainer)
			{
				if(spec.getAnnotation() == null)
					continue;
				if(spec.getCharge() != charge)
					continue;
				
				float curParentMass = spec.getParentMass();
				if(curParentMass < parentMassRange.getFirst() || curParentMass >= parentMassRange.getSecond())
					continue;
				
				Peptide annotation = spec.getAnnotation();
				numCleavages += annotation.size()-1;
				numSpec++;
				spec = filter.apply(spec);
				
				for(int c=1; c<=charge; c++)
				{
					for(int direction=0; direction<2; direction++)
					{
						double accurateMass = 0;
						Hashtable<Integer, Histogram<Integer>> ionFreq = null;
						for(int i=0; i<annotation.size()-1; i++)
						{
							if(direction == 0)
							{
								accurateMass += annotation.get(i).getAccurateMass();
								ionFreq = prefixIonFreq;
							}
							else if(direction == 1)
							{
								accurateMass += annotation.get(annotation.size()-1-i).getAccurateMass();
								ionFreq = suffixIonFreq;
							}
							float mass = (float)(accurateMass/c);
							ArrayList<Peak> peakList = spec.getPeakListByMassRange(
									mass+MIN_OFFSET_MASS/(float)c-mme.getToleranceAsDa(mass), 
									mass+MAX_OFFSET_MASS/(float)c+mme.getToleranceAsDa(mass));
							int prevIntOffset = Integer.MIN_VALUE;
							for(Peak p : peakList)
							{
								float peakMz = p.getMz();
								int segNum = getSegmentNum(peakMz, curParentMass);
								if(segNum != seg)
									continue;
								float offset = peakMz - mass;
								int intOffset = NominalMass.toNominalMass(offset);
								if(intOffset > prevIntOffset)
								{
									ionFreq.get(c).add(intOffset);
									prevIntOffset = intOffset;
								}
							}
						}
					}
				}
			}
	 		
			float maxProb = 0;
			int maxCharge = 0;
			int maxDirection = 0;
			float maxOffset = 0;
			
			ArrayList<FragmentOffsetFrequency> fragmentOffsetFrequencyList = new ArrayList<FragmentOffsetFrequency>();
			ArrayList<FragmentOffsetFrequency> insignificantFragmentOffsetFrequencyList = new ArrayList<FragmentOffsetFrequency>();
	 		for(int c=1; c<=charge; c++)
	 		{
	 			for(int direction=0; direction<2; direction++)
	 			{
		 	 		ArrayList<Integer> keyList;
		 	 		if(direction == 0)
		 	 			keyList = new ArrayList<Integer>(prefixIonFreq.get(c).keySet());
		 	 		else
		 	 			keyList = new ArrayList<Integer>(suffixIonFreq.get(c).keySet());
		 	 		
		 	 		Collections.sort(keyList);
		 	 		for(Integer key : keyList)
		 	 		{
		 	 			float offset = NominalMass.getMassFromNominalMass(key);
			 	 		int freq;
			 	 		if(direction == 0)
			 	 			freq = prefixIonFreq.get(c).get(key);
			 	 		else
			 	 			freq = suffixIonFreq.get(c).get(key);
		 	 			float prob = freq/(float)numCleavages*numSegments;
		 	 			if(prob > maxProb)
		 	 			{
		 	 				maxProb = prob;
		 	 				maxCharge = c;
		 	 				maxDirection = direction;
		 	 				maxOffset = offset;
		 	 			}
		 	 			if(prob > minProbThreshold)
		 	 			{
		 	 				if(direction == 0)
		 	 					fragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(new IonType.PrefixIon(c, offset), prob));
		 	 				else
		 	 					fragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(new IonType.SuffixIon(c, offset), prob));
		 	 			}
		 	 			else
		 	 			{
		 	 				if(direction == 0)
		 	 					insignificantFragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(new IonType.PrefixIon(c, offset), prob));
		 	 				else
		 	 					insignificantFragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(new IonType.SuffixIon(c, offset), prob));
		 	 			}
		 	 		}
	 			}
	 		}			

	 		if(fragmentOffsetFrequencyList.size() == 0)
	 		{
	 			if(maxDirection == 0)
	 				fragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(new IonType.PrefixIon(maxCharge, maxOffset), maxProb));
	 			else
	 				fragmentOffsetFrequencyList.add(new FragmentOffsetFrequency(new IonType.SuffixIon(maxCharge, maxOffset), maxProb));
	 		}
	 		
	 		Collections.sort(insignificantFragmentOffsetFrequencyList);
	 		ArrayList<FragmentOffsetFrequency> noiseOffsetFrequencyList = new ArrayList<FragmentOffsetFrequency>(NUM_NOISE_IONS);
	 		
	 		int numNoise = 0;
	 		for(FragmentOffsetFrequency off : insignificantFragmentOffsetFrequencyList)
	 		{
	 			if(off.getIonType().getCharge() == 1)
	 				noiseOffsetFrequencyList.add(off);
	 			if(++numNoise >= NUM_NOISE_IONS)
	 				break;
	 		}
	 		Collections.sort(fragmentOffsetFrequencyList, Collections.reverseOrder());
	 		fragOFFTable.put(partition, fragmentOffsetFrequencyList);
	 		insignificantFragOFFTable.put(partition, noiseOffsetFrequencyList);
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