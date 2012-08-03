package edu.ucsd.msjava.msscorer;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;

import edu.ucsd.msjava.msgf.Histogram;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.IonType;
import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.WindowFilter;
import edu.ucsd.msjava.msutil.IonType.PrefixIon;

public class NewRankScorer implements NewAdditiveScorer {
	public static final int VERSION = 7061;
	public static final String DATE = "12/21/2011";
	// Optional
	protected WindowFilter filter = new WindowFilter(6, 50);
	
	// Type of the data
	protected SpecDataType dataType;
	
	// Parameters to be used for scoring
	protected int numSegments = 1;
	protected Histogram<Integer> chargeHist = null;
	protected TreeSet<Partition> partitionSet = null;
	protected TreeMap<Integer,ArrayList<PrecursorOffsetFrequency>> precursorOFFMap = null;	// charge -> precursorOffsetList
	protected Hashtable<Partition,ArrayList<FragmentOffsetFrequency>> fragOFFTable = null;	// partition -> ionTypes
	protected Hashtable<Partition,ArrayList<FragmentOffsetFrequency>> insignificantFragOFFTable = null;	// for noise error distribution
	protected Hashtable<Partition,Hashtable<IonType,Float[]>> rankDistTable = null;
	
	protected Tolerance mme = new Tolerance(0.5f);
	
	// Deconvolution
	protected boolean applyDeconvolution = false;
	protected float deconvolutionErrorTolerance = 0;
	
	protected int numPrecurOFF = 0;
	protected int maxRank = 0;
	
	// For edge scoring
	protected int errorScalingFactor = 0;	// if 0, don't user errors, 10 for low accuracy, 100 for high accuracy
	protected Hashtable<Partition,Float[]> ionErrDistTable = null;
	protected Hashtable<Partition,Float[]> noiseErrDistTable = null;
	protected Hashtable<Partition,Float[]> ionExistenceTable = null;
	
	// Ion Types
	private HashMap<Partition, IonType> mainIonTable;
	private HashMap<Partition, IonType[]> ionTypeTable;

	public NewRankScorer() 
	{
	}
	
	public NewRankScorer(String paramFileName)
	{
		readFromFile(new File(paramFileName), false);
	}
	
	public NewRankScorer(InputStream is)
	{
		readFromInputStream(is, false);
	}
	
	public<T extends Matter> NewScoredSpectrum<T>  getScoredSpectrum(Spectrum spec)
	{
		return new NewScoredSpectrum<T>(spec, this);
	}
	
	public SpecDataType getSpecDataType()
	{
		return dataType;
	}
	
	public TreeSet<Partition> getParitionSet()
	{
		return partitionSet;
	}
	
	public void filterPrecursorPeaks(Spectrum spec)
	{
		for(PrecursorOffsetFrequency off : getPrecursorOFF(spec.getCharge()))
			spec.filterPrecursorPeaks(mme, off.getReducedCharge(), off.getOffset());
	}
	
	public NewRankScorer mme(Tolerance mme)
	{
		this.mme = mme;
		return this;
	}
	
	public boolean applyDeconvolution()
	{
		return this.applyDeconvolution;
	}
	
	public float deconvolutionErrorTolerance()
	{
		return this.deconvolutionErrorTolerance;
	}
	
	public NewRankScorer doNotUseError()	{ this.errorScalingFactor = 0; return this;}
	
	public boolean supportEdgeScores()
	{
		return errorScalingFactor != 0;
	}
	
	public float getNodeScore(Partition part, IonType ionType, int rank) {
		// ion score
		Hashtable<IonType,Float[]> rankTable = rankDistTable.get(part);	// rank -> probability
		assert(rankTable != null);
		int rankIndex = rank > maxRank ? maxRank-1 : rank-1;
		float ionScore = getScoreFromTable(rankIndex, rankTable, ionType, false);

		return ionScore;
	}

	public float getMissingIonScore(Partition part, IonType ionType) {
		Hashtable<IonType,Float[]> table = rankDistTable.get(part);
		assert(table != null);
		int rankIndex = maxRank;
		return getScoreFromTable(rankIndex, table, ionType, false);
	}
	
	public float getErrorScore(Partition part, float error) {
		int errIndex = Math.round(error*errorScalingFactor);
		if(errIndex > errorScalingFactor)
			errIndex = errorScalingFactor;
		else if(errIndex < -errorScalingFactor)
			errIndex = -errorScalingFactor;
		Float[] ionErrHist = this.ionErrDistTable.get(part);
//		float noiseProb = (errorScalingFactor-Math.abs(errIndex))/(errorScalingFactor*errorScalingFactor);
//		if(noiseProb == 0)
//			noiseProb = 1f/(errorScalingFactor*errorScalingFactor);
//		return (float)Math.log(ionErrHist[errIndex+errorScalingFactor]/noiseProb);
		errIndex += errorScalingFactor;
		Float[] noiseErrHist = this.noiseErrDistTable.get(part);
		return (float)Math.log(ionErrHist[errIndex]/noiseErrHist[errIndex]);
	}	
	
	public float getIonExistenceScore(Partition part, int index, float probPeak)
	{
		Float[] ionExistenceProb = this.ionExistenceTable.get(part);
		float noiseExistenceProb;
		if(index == 0)	// nn
			noiseExistenceProb = (1-probPeak)*(1-probPeak);
		else if(index == 3) // yy
			noiseExistenceProb = probPeak*probPeak;
		else
			noiseExistenceProb = probPeak*(1-probPeak);
		return (float)Math.log(ionExistenceProb[index]/noiseExistenceProb);
	}
	
	private float getScoreFromTable(int index, Hashtable<IonType,Float[]> table, IonType ionType, boolean isError)
	{
		Float[] frequencies = table.get(ionType);
		assert(frequencies != null): ionType.getName()+" is not supported!";
		float ionFrequency = frequencies[index];
		Float[] noiseFrequencies = table.get(IonType.NOISE);
		assert(noiseFrequencies != null);
		float noiseFrequency = noiseFrequencies[index];
		if(!isError)
			noiseFrequency *= Math.min(ionType.getCharge(), numSegments);
		assert(ionFrequency > 0 && noiseFrequency > 0): "Ion frequency must be positive:" +
				index+" "+ionType.getName()+" "+ionFrequency+" "+noiseFrequency;
		return (float)Math.log(ionFrequency/noiseFrequency);
	}

	public void readFromFile(File paramFile)
	{
		readFromFile(paramFile, false);
	}
	
	protected void readFromFile(File paramFile, boolean verbose)
	{
		InputStream is = null;
		try {
			is = new BufferedInputStream(new FileInputStream(paramFile));
		} catch (IOException e) {
			e.printStackTrace();
		}
		readFromInputStream(is, verbose);
	}
	
	private void readFromInputStream(InputStream is, boolean verbose)
	{
		DataInputStream in = new DataInputStream(is);

		// Read the date
		try {
//			int year = in.readInt();	// version information
//			int month = in.readInt();
//			int date = in.readInt();
//			if(verbose)
//			System.out.println("CreationDate: " + year + "/" + (month+1) + "/" + date);
			
			int version = in.readInt();
			if(verbose)
				System.out.println("Version: " + version);

			// Read activation method
			StringBuffer bufMet = new StringBuffer();
			byte lenActMethod = in.readByte();
			for(byte i=0; i<lenActMethod; i++)
				bufMet.append(in.readChar());
			ActivationMethod activationMethod = ActivationMethod.get(bufMet.toString());
			
			// Read instrument type
			StringBuffer bufInst = new StringBuffer();
			byte lenInst = in.readByte();
			for(byte i=0; i<lenInst; i++)
				bufInst.append(in.readChar());
			InstrumentType instType = InstrumentType.get(bufInst.toString());

			// Read enzyme
			Enzyme enzyme;
			StringBuffer bufEnz = new StringBuffer();
			byte lenEnz = in.readByte();
			if(lenEnz != 0)
			{
				for(byte i=0; i<lenEnz; i++)
					bufEnz.append(in.readChar());
				enzyme = Enzyme.getEnzymeByName(bufEnz.toString());
			}
			else
				enzyme = null;
			
			// Read protocol
			Protocol protocol;
			StringBuffer bufProtocol = new StringBuffer();
			byte lenProtocol = in.readByte();
			if(lenProtocol != 0)
			{
				for(byte i=0; i<lenProtocol; i++)
					bufProtocol.append(in.readChar());
				protocol = Protocol.get(bufProtocol.toString());
			}
			else
				protocol = Protocol.NOPROTOCOL;
			
			this.dataType = new SpecDataType(activationMethod, instType, enzyme, protocol);
			
			// MME
			boolean isTolerancePPM = in.readBoolean();
			float mmeVal = in.readFloat();
			mme = new Tolerance(mmeVal, isTolerancePPM);
			
			// Apply deconvolution
			boolean applyDeconvolution = in.readBoolean();
			float deconvolutionErrorTolerance = in.readFloat();
			this.applyDeconvolution = applyDeconvolution;
			this.deconvolutionErrorTolerance = deconvolutionErrorTolerance;
			
			// Charge histogram
			if(verbose)
				System.out.println("ChargeHistogram");
			chargeHist = new Histogram<Integer>();
			int size = in.readInt();	// size
			for(int i=0; i<size; i++)
			{
				int charge = in.readInt();
				int numSpecs = in.readInt();
				if(verbose)
					System.out.println(charge+"\t"+numSpecs);
				chargeHist.put(charge, numSpecs);
			}

			// Partition info
			if(verbose)
				System.out.println("PartitionInfo");
			partitionSet = new TreeSet<Partition>();
			size = in.readInt();
			numSegments = in.readInt();
			for(int i=0; i<size; i++)
			{
				int charge = in.readInt();
				float parentMass = in.readFloat();
				int segNum = in.readInt();
				partitionSet.add(new Partition(charge, parentMass, segNum));
				if(verbose)
					System.out.println(charge+"\t"+parentMass+"\t"+segNum);
			}

			// Precursor offset frequency function
			if(verbose)
				System.out.println("PrecursorOFF");
			precursorOFFMap = new TreeMap<Integer,ArrayList<PrecursorOffsetFrequency>>();
			size = in.readInt();
			this.numPrecurOFF = size;
			for(int i=0; i<size; i++)
			{
				int charge = in.readInt();
				int reducedCharge = in.readInt();
				float offset = in.readFloat();
				boolean isTolPPM = in.readBoolean();
				float tolVal = in.readFloat();
				
				float frequency = in.readFloat();
				ArrayList<PrecursorOffsetFrequency> offList = precursorOFFMap.get(charge);
				if(offList == null)
				{
					offList = new ArrayList<PrecursorOffsetFrequency>();
					precursorOFFMap.put(charge, offList);
				}
				offList.add(new PrecursorOffsetFrequency(reducedCharge, offset, frequency).tolerance(new Tolerance(tolVal, isTolPPM)));
				if(verbose)
					System.out.println(charge+"\t"+reducedCharge+"\t"+offset+"\t"+new Tolerance(tolVal, isTolPPM).toString()+"\t"+frequency);
			}

			// Fragment ion offset frequency function
			if(verbose)
				System.out.println("FragmentOFF");
			fragOFFTable = new Hashtable<Partition,ArrayList<FragmentOffsetFrequency>>();
			for(Partition partition : partitionSet)
			{
				if(verbose)
					System.out.println(partition.getCharge()+"\t"+partition.getSegNum()+"\t"+partition.getParentMass());
				ArrayList<FragmentOffsetFrequency> fragmentOFF = new ArrayList<FragmentOffsetFrequency>();
				size = in.readInt();
				for(int i=0; i<size; i++)
				{
					boolean isPrefix = in.readBoolean();
					int charge = in.readInt();
					float offset = in.readFloat();
					IonType ionType;
					if(isPrefix)
						ionType = new IonType.PrefixIon("P_"+charge+"_"+Math.round(offset),charge, offset);
					else
						ionType = new IonType.SuffixIon("S_"+charge+"_"+Math.round(offset),charge, offset);
					float frequency = in.readFloat();
					fragmentOFF.add(new FragmentOffsetFrequency(ionType, frequency));
					if(verbose)
						System.out.println(ionType.getName()+"\t"+frequency);
				}
				fragOFFTable.put(partition, fragmentOFF);
			}

			determineIonTypes();
			// Rank distributions
			rankDistTable = new Hashtable<Partition,Hashtable<IonType,Float[]>>(); 
			maxRank = in.readInt();
			if(verbose)
				System.out.println("RankDistribution,"+maxRank);
			for(Partition partition : partitionSet)
			{
				if(verbose)
					System.out.println(partition.getCharge()+"\t"+partition.getSegNum()+"\t"+partition.getParentMass());
				Hashtable<IonType,Float[]> table = new Hashtable<IonType,Float[]>();
				ArrayList<IonType> ionTypeList = new ArrayList<IonType>();
				IonType[] ionTypes = getIonTypes(partition);
				if(ionTypes == null || ionTypes.length == 0)
					continue;
				
				for(IonType ion : ionTypes)
					ionTypeList.add(ion);
				ionTypeList.add(IonType.NOISE);
				for(IonType ion : ionTypeList)
				{
					if(verbose)
						System.out.print(ion.getName());
					Float[] frequencies = new Float[maxRank+1];
					for(int i=0; i<frequencies.length; i++)
					{
						frequencies[i] = in.readFloat();
						if(verbose)
							System.out.print("\t"+frequencies[i]);
					}
					table.put(ion, frequencies);
					if(verbose)
						System.out.println();
				}
				rankDistTable.put(partition, table);
			}
			
			// Error distribution

			errorScalingFactor = in.readInt();
			if(errorScalingFactor > 0)
			{
				if(verbose)
					System.out.println("ErrorDistribution,"+errorScalingFactor);
				
				ionErrDistTable = new Hashtable<Partition,Float[]>();
				noiseErrDistTable = new Hashtable<Partition,Float[]>();
				ionExistenceTable = new Hashtable<Partition,Float[]>();
				
				for(Partition partition : partitionSet)
				{
					if(verbose)
						System.out.println(partition.getCharge()+"\t"+partition.getSegNum()+"\t"+partition.getParentMass());
					Float[] ionErrDist = new Float[errorScalingFactor*2+1];
					for(int i=0; i<ionErrDist.length; i++)
						ionErrDist[i] = in.readFloat();
					ionErrDistTable.put(partition, ionErrDist);
					Float[] noiseErrDist = new Float[errorScalingFactor*2+1];
					for(int i=0; i<noiseErrDist.length; i++)
						noiseErrDist[i] = in.readFloat();
					noiseErrDistTable.put(partition, noiseErrDist);
					Float[] ionExTable = new Float[4];
					for(int i=0; i<ionExTable.length; i++)
						ionExTable[i] = in.readFloat();
					ionExistenceTable.put(partition, ionExTable);
				}				
			}
			
			int validation = in.readInt();
			if(validation != Integer.MAX_VALUE)
			{
				System.err.println("Parameter is wrong!");
				System.exit(-1);
			}
			in.close();		
		} catch (IOException e) {
			e.printStackTrace();
		}
	}	
	
	// Builders
	protected NewRankScorer tolerance(Tolerance mme)
	{
		this.mme = mme;
		return this;
	}

	protected NewRankScorer filter(WindowFilter filter)
	{
		this.filter = filter;
		return this;
	}
	
	// Getters and Setters
	public Tolerance getMME()	{ return mme; }
	
	protected Histogram<Integer> getChargeHist() {
		return chargeHist;
	}

	protected TreeSet<Partition> getPartitionSet() {
		return partitionSet;
	}

	protected int getNumPrecursorOFF() {	return this.numPrecurOFF;	}
	protected int getMaxRank()	{ return this.maxRank; }
	protected int getNumErrorBins()	{ return this.errorScalingFactor; }
	protected int getNumSegments() {	return this.numSegments;	}
	
	int getSegmentNum(float peakMz, float parentMass)
	{
		int segNum = (int)(peakMz/parentMass*numSegments);
		if(segNum >= numSegments)
			segNum = numSegments-1;
		return segNum;
	}
	
	protected ArrayList<PrecursorOffsetFrequency> getPrecursorOFF(int charge) 
	{
		if(precursorOFFMap == null || precursorOFFMap.size() == 0)
			return new ArrayList<PrecursorOffsetFrequency>();
		Entry<Integer, ArrayList<PrecursorOffsetFrequency>> entry = precursorOFFMap.floorEntry(charge);
		if(entry == null)
			entry = precursorOFFMap.ceilingEntry(charge);
		return entry.getValue();
	}
	
	protected Partition getPartition(int charge, float parentMass, int segNum)
	{
		if(partitionSet == null || partitionSet.size() == 0)
			return null;
		Partition partition = new Partition(charge, parentMass, segNum);
		Partition matched = partitionSet.floor(partition);
		if(matched == null)	// small charge
		{
			// use the smallest charge available
			partition = new Partition(partitionSet.first().getCharge(), parentMass, segNum);
			return partitionSet.floor(partition);
		}
		if(charge == matched.getCharge())	// scoring is available at this charge
		{
			return matched;
		}
		else	// high charge
		{
			partition = new Partition(matched.getCharge(), parentMass, segNum);
			return partitionSet.floor(partition);
		}
	}
	
	protected ArrayList<FragmentOffsetFrequency> getFragmentOFF(int charge, float parentMass, int segNum)
	{
		return getFragmentOFF(getPartition(charge, parentMass, segNum));
	}

	protected ArrayList<FragmentOffsetFrequency> getFragmentOFF(Partition partition)
	{
		return this.fragOFFTable.get(partition);
	}
	
	protected Hashtable<IonType,Float[]> getRankDistTable(int charge, float parentMass, int segNum)
	{
		return getRankDistTable(getPartition(charge, parentMass, segNum));
	}
	
	protected Hashtable<IonType,Float[]> getRankDistTable(Partition partition)
	{
		return this.rankDistTable.get(partition);
	}
	
	public IonType[] getIonTypes(int charge, float parentMass, int segNum)
	{
		return getIonTypes(getPartition(charge, parentMass, segNum));
	}
	
	protected IonType[] getIonTypes(Partition partition)
	{
		if(ionTypeTable != null)
			return ionTypeTable.get(partition);

		else
		{
			ArrayList<FragmentOffsetFrequency> offList = fragOFFTable.get(partition);
			IonType[] ionTypes = new IonType[offList.size()];
			for(int i=0; i<offList.size(); i++)
				ionTypes[i] = offList.get(i).getIonType();
			return ionTypes;
		}
	}	
	
	protected IonType getMainIonType(Partition partition)
	{
		return mainIonTable.get(partition);
	}
	
	protected void determineIonTypes()
	{
		ionTypeTable = new HashMap<Partition, IonType[]>();
		
		for(Partition partition : partitionSet)
		{
			ArrayList<FragmentOffsetFrequency> offList = fragOFFTable.get(partition);
			IonType[] ionTypes = new IonType[offList.size()];
			for(int i=0; i<offList.size(); i++)
				ionTypes[i] = offList.get(i).getIonType();
			ionTypeTable.put(partition, ionTypes);
		}
		
		mainIonTable = new HashMap<Partition, IonType>();
		for(Partition partition : partitionSet)
		{
			if(partition.getSegNum() != 0)
				continue;
			HashMap<IonType, Float> ionProb = new HashMap<IonType, Float>();
			for(int seg=0; seg<numSegments; seg++)
			{
				Partition part = new Partition(partition.getCharge(), partition.getParentMass(), seg);
				ArrayList<FragmentOffsetFrequency> offList = fragOFFTable.get(part);
				for(FragmentOffsetFrequency off : offList)
				{
					Float prob = ionProb.get(off.getIonType());
					if(prob == null)
						ionProb.put(off.getIonType(), off.getFrequency());
					else
						ionProb.put(off.getIonType(), prob+off.getFrequency());
				}
			}
			IonType mainIon = null;
			float prob = -1;
			for(IonType ion : ionProb.keySet())
			{
				if(ionProb.get(ion) > prob)
				{
					mainIon = ion;
					prob = ionProb.get(ion);
				}
			}
			assert(mainIon != null);
			for(int seg=0; seg<numSegments; seg++)
			{
				Partition part = new Partition(partition.getCharge(), partition.getParentMass(), seg);
				mainIonTable.put(part, mainIon);
			}
		}
	}

	protected HashSet<Integer> getIonOffsets(Partition partition, int charge, boolean isPrefix)
	{
		HashSet<Integer> offsets = new HashSet<Integer>();
		ArrayList<FragmentOffsetFrequency> offList = fragOFFTable.get(partition);
		for(FragmentOffsetFrequency off : offList)
		{
			if(isPrefix && (off.getIonType() instanceof IonType.PrefixIon)
				|| !isPrefix && (off.getIonType() instanceof IonType.SuffixIon))
			{
				offsets.add(Math.round(off.getIonType().getOffset()));
			}
		}
		return offsets;
	}	
	
	protected IonType[] getNoiseIonTypes(Partition partition)
	{
		ArrayList<FragmentOffsetFrequency> offList = insignificantFragOFFTable.get(partition);
		IonType[] ionTypes = new IonType[offList.size()];
		for(int i=0; i<offList.size(); i++)
			ionTypes[i] = offList.get(i).getIonType();
		return ionTypes;
	}	

	protected void writeParameters(File outputFile)
	{
		if(chargeHist == null ||
				partitionSet == null ||
				precursorOFFMap == null ||
				fragOFFTable == null ||
				rankDistTable == null)
		{
			assert(false): "Parameters are not generated!";
			System.exit(-1);
			return;
		}
		
		DataOutputStream out = null;
		try {
			out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		} catch (IOException e) {
			e.printStackTrace();
		}

		// Write the date
		try {
			out.writeInt(VERSION);
			
			// Write method
			out.writeByte(dataType.getActivationMethod().getName().length());
			out.writeChars(dataType.getActivationMethod().getName());
			
			// Write instrument type
			out.writeByte(dataType.getInstrumentType().getName().length());
			out.writeChars(dataType.getInstrumentType().getName());
			
			// Write enzyme
			Enzyme enzyme = dataType.getEnzyme();
			if(enzyme != null)
			{
				out.writeByte(enzyme.getName().length());
				out.writeChars(enzyme.getName());
			}
			else
				out.writeByte((byte)0);
	
			// Write protocol
			Protocol protocol = dataType.getProtocol();
			if(protocol != null && protocol != Protocol.NOPROTOCOL)
			{
				out.writeByte(protocol.getName().length());
				out.writeChars(protocol.getName());
			}
			else
				out.writeByte((byte)0);
			
			// Maximum mass error
			out.writeBoolean(mme.isTolerancePPM());
			out.writeFloat(mme.getValue());
			
			// Apply deconvolution
			out.writeBoolean(applyDeconvolution);
			out.writeFloat(deconvolutionErrorTolerance);
			
			// Charge histogram
			out.writeInt((chargeHist.maxKey()-chargeHist.minKey()+1));	// size
			for(int charge=chargeHist.minKey(); charge<=chargeHist.maxKey(); charge++)
			{
				out.writeInt(charge);
				out.writeInt(chargeHist.get(charge));
			}

			// Partition info
			out.writeInt(partitionSet.size());
			out.writeInt(numSegments);
			for(Partition p : partitionSet)
			{
				out.writeInt(p.getCharge());
				out.writeFloat(p.getParentMass());
				out.writeInt(p.getSegNum());
			}

			// Precursor offset frequency function
			out.writeInt(numPrecurOFF);
			for(int charge=chargeHist.minKey(); charge<=chargeHist.maxKey(); charge++)
			{
				ArrayList<PrecursorOffsetFrequency> offList = precursorOFFMap.get(charge);
				if(offList != null) 
				{
					for(PrecursorOffsetFrequency off : offList)
					{
						out.writeInt(charge);	// charge
						out.writeInt(off.getReducedCharge());	// reduced charge
						out.writeFloat(off.getOffset());	// offset
						out.writeBoolean(off.getTolerance().isTolerancePPM());
						out.writeFloat(off.getTolerance().getValue());
						out.writeFloat(off.getFrequency());	// frequency
					}
				}
			}

			// Fragment ion offset frequency function
			for(Partition partition : partitionSet)
			{
				ArrayList<FragmentOffsetFrequency> fragmentOFF = getFragmentOFF(partition);
				out.writeInt(fragmentOFF.size());	// num offsets
				Collections.sort(fragmentOFF, Collections.reverseOrder());
				for(FragmentOffsetFrequency off : fragmentOFF)
				{
					out.writeBoolean(off.getIonType() instanceof PrefixIon);
					out.writeInt(off.getIonType().getCharge());
					out.writeFloat(off.getIonType().getOffset());
					out.writeFloat(off.getFrequency());
				}
			}

			// Rank distributions
			out.writeInt(maxRank);
			for(Partition partition : partitionSet)
			{
//				if(partition.getParentMass() > 4100 && partition.getCharge() == 5 && partition.getSegNum() == 1)
//					System.out.println("Debug");
				
				Hashtable<IonType,Float[]> rankDistTable = getRankDistTable(partition);
				if(rankDistTable == null)
					continue;
				IonType[] ionTypes = getIonTypes(partition);
				if(ionTypes == null || ionTypes.length == 0)
					continue;
				ArrayList<IonType> ionTypeList = new ArrayList<IonType>();
				for(IonType ion : ionTypes)
					ionTypeList.add(ion);
				ionTypeList.add(IonType.NOISE);
				for(IonType ion : ionTypeList)
				{
					Float[] frequencies = rankDistTable.get(ion);
					assert(frequencies.length == maxRank+1);
					for(Float freq : frequencies)
						out.writeFloat(freq);
				}
			}
			
			// Error distribution
//			protected int errorScalingFactor = 0;	// if 0, don't user errors, 10 for low accuracy, 100 for high accuracy
//			protected Hashtable<Partition,Float[]> ionErrDistTable = null;
//			protected Hashtable<Partition,Float[]> noiseErrDistTable = null;
//			protected Hashtable<Partition,Float[]> ionExistenceTable = null;
			out.writeInt(errorScalingFactor);
			if(errorScalingFactor > 0)
			{
				for(Partition partition : partitionSet)
				{
					Float[] ionErrDist = ionErrDistTable.get(partition);
					assert(ionErrDist.length == 2*errorScalingFactor+1);
					for(Float f : ionErrDist)
						out.writeFloat(f);
					Float[] noiseErrDist = noiseErrDistTable.get(partition);
					assert(noiseErrDist.length == 2*errorScalingFactor+1);
					for(Float f : noiseErrDist)
						out.writeFloat(f);
					Float[] ionExTable = ionExistenceTable.get(partition);
					assert(ionExTable.length == 4);
					for(Float f : ionExTable)
						out.writeFloat(f);
				}	
			}
			
			// for validation
			out.writeInt(Integer.MAX_VALUE);
			out.flush();
			out.close();		
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	protected void writeParametersPlainText(File outputFile)
	{
		PrintStream out = null;
		if(outputFile == null)
			out = System.out;
		else
		{
			try {
				out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		// Write the version info
		out.println("#MSGFScoringParameters\tv"+	
				new SimpleDateFormat("yyyyMMdd").format(Calendar.getInstance().getTime()));

		// Write method
		if(dataType.getActivationMethod() != null)
			out.println("#Activation Method: " + dataType.getActivationMethod().getName());
		
		// Write instrument type
		if(dataType.getInstrumentType() != null)
			out.println("#Instrument type: " + dataType.getInstrumentType().getName());
		
		// Write enzyme
		if(dataType.getEnzyme() != null)
			out.println("#Enzyme: " + dataType.getEnzyme().getName());
		
		// Write protocol
		if(dataType.getProtocol() != null)
			out.println("#Protocol: " + dataType.getProtocol().getName());
		
		// Write mme
		out.println("#Maximum mass error: " + mme.toString());
		
		// Write whether to apply deconvolution
		out.println("Apply deconvolution: " + applyDeconvolution);
		out.println("Deconvolution error tolerance: " + deconvolutionErrorTolerance);
		
		// Charge histogram
		out.println("#ChargeHistogram\t"+(chargeHist.maxKey()-chargeHist.minKey()+1));
		for(int charge=chargeHist.minKey(); charge<=chargeHist.maxKey(); charge++)
			out.println(charge+"\t"+chargeHist.get(charge));
		
		// Partition info
		out.println("#Partitions\t"+partitionSet.size());
		for(Partition p : partitionSet)
			out.println(p.getCharge()+"\t"+p.getSegNum()+"\t"+p.getParentMass());
		
		// Precursor offset frequency function
		out.println("#PrecursorOffsetFrequencyFunction\t"+numPrecurOFF);
		for(int charge=chargeHist.minKey(); charge<=chargeHist.maxKey(); charge++)
		{
			ArrayList<PrecursorOffsetFrequency> offList = precursorOFFMap.get(charge);
			if(offList != null)
				for(PrecursorOffsetFrequency off : offList)
					out.println(charge+"\t"+off.getReducedCharge()+"\t"+off.getOffset()+"\t"+off.getTolerance().toString()+"\t"+off.getFrequency());
		}

		// Fragment ion offset frequency function
		out.println("#FragmentOffsetFrequencyFunction\t"+partitionSet.size());
		for(Partition partition : partitionSet)
		{
			ArrayList<FragmentOffsetFrequency> fragmentOFF = getFragmentOFF(partition);
			out.println("Partition\t"+partition.getCharge()+"\t"+partition.getSegNum()+"\t"+partition.getParentMass()+"\t"+fragmentOFF.size());
			Collections.sort(fragmentOFF, Collections.reverseOrder());
			for(FragmentOffsetFrequency off : fragmentOFF)
				out.println(off.getIonType().getName()+"\t"+off.getFrequency()+"\t"+off.getIonType().getOffset());
		}
		
		// Rank distributions
		out.println("#RankDistributions\t"+partitionSet.size());
		for(Partition partition : partitionSet)
		{
			Hashtable<IonType,Float[]> rankDistTable = getRankDistTable(partition);
			IonType[] ionTypes = getIonTypes(partition);
			if(ionTypes == null || ionTypes.length == 0)
				continue;
			ArrayList<IonType> ionTypeList = new ArrayList<IonType>();
			for(IonType ion : ionTypes)
				ionTypeList.add(ion);
			ionTypeList.add(IonType.NOISE);
			out.println("Partition\t"+partition.getCharge()+"\t"+partition.getSegNum()+"\t"+partition.getParentMass()+"\t"+ionTypeList.size()+"\t"+maxRank);
			for(IonType ion : ionTypeList)
			{
				out.print(ion.getName());
				Float[] frequencies = rankDistTable.get(ion);
				for(Float freq : frequencies)
					out.print("\t"+freq);
				out.println();
			}
		}
		
		// Error distributions
		// Error distribution
		if(errorScalingFactor > 0)
		{
			out.println("#ErrorDistributions\t"+errorScalingFactor);
			for(Partition partition : partitionSet)
			{
				out.println("Partition\t"+partition.getCharge()+"\t"+partition.getSegNum()+"\t"+partition.getParentMass()+"\t"+this.getMainIonType(partition).getName());
				Float[] ionErrDist = ionErrDistTable.get(partition);
				out.print("Signal");
				for(Float f : ionErrDist)
					out.print("\t"+f);
				out.println();
				Float[] noiseErrDist = noiseErrDistTable.get(partition);
				out.print("Noise");
				for(Float f : noiseErrDist)
					out.print("\t"+f);
				out.println();
				Float[] ionExTable = ionExistenceTable.get(partition);
				out.print("IonExistence");
				for(Float f : ionExTable)
					out.print("\t"+f);
				out.println();
			}	
		}
		
		out.flush();
		out.close();
	}	
	
	public static void main(String argv[]) throws Exception
	{
		readWriteTest();
//		paramTest();
	}
	
	public static void readWriteTest() throws Exception
	{
	}
}
