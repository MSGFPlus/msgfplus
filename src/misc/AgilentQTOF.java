package misc;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import parser.MgfSpectrumParser;
import msgf.Histogram;
import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.CompositionFactory;
import msutil.Constants;
import msutil.Enzyme;
import msutil.IonType;
import msutil.Peak;
import msutil.Peptide;
import msutil.Sequence;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import msutil.IonType.PrefixIon;

public class AgilentQTOF {
	public static void main(String argv[]) throws Exception
	{
//		correctChargeZero();
//		compositionMSGF();
//		aminoacidMSGF();
//		compositionDensity();
//		ionProb();
//		graphTest();
//		specGraphMSGF();
//		compositionClusterMSGF();
//		testFixedLengthAAGraph();
//		testFloatMassMSGF();
//		testGraphErrors();
//		testGraphTheoreticalErrors();
//		testGraphErrorDist();
//		testErrors();
//		aminoacidExtMSGF();
//		trypsinTest();
	}

	public static void trypsinTest() throws Exception
	{
		SpectraIterator itr = new SpectraIterator("/home/sangtaekim/Research/Data/HeckRevision/AnnotatedSpectra/CID_Tryp_Confident.mgf", new MgfSpectrumParser());
		int[] numSpecs = new int[100];
		int[] numSpecsPrecedingKR = new int[100];
		int[] numSpecsEndingKR = new int[100];
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int charge = spec.getCharge();
			numSpecs[charge]++;
			String annotation = spec.getAnnotationStr();
//			if(annotation.charAt(0) == 'K' || annotation.charAt(0) == 'R')
//			{
//				
//			}
			char lastChar = annotation.charAt(annotation.length()-1);
			if(lastChar == 'K' || lastChar == 'R')
				numSpecsEndingKR[charge]++;
		} 
		
		for(int c=2; c<=7; c++)
			System.out.println(c+"\t"+numSpecsEndingKR[c]/(float)numSpecs[c]);
	}
	
//	public static void testErrors() throws Exception
//	{
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
////		int numBits = 15;
//		
//		for(int numBits=17; numBits>=17; numBits--)
//		{
//			FloatingPointMassFactory factory = new FloatingPointMassFactory(aaSet, numBits, 1);
//			ArrayList<FloatingPointMass> aaNodes = factory.getAllNodes();
//			float errSum = 0;
//			for(int i=1; i<aaNodes.size(); i++)
//			{
//				FloatingPointMass m = aaNodes.get(i);
//				AminoAcid aa = aaSet.getAminoAcid(factory.getEdgeIndex(m, factory.getZero()));
//				float err = m.getMass()-aa.getMass();
//				errSum += err;
//				System.out.println(m.getMass()+"\t"+m.toBinaryString()+"\t"+err/m.getMass()*1e6f);
//			}
//			float avgErr = errSum/(aaNodes.size()-1);
//			float avgErrPPM = avgErr/(float)Math.pow(2, 20-numBits)*1e6f;
//			System.out.println(numBits+"\t"+avgErr+"\t"+avgErrPPM);
//		}
//		
//		System.out.println(new FloatingPointMass(0, 17).getMass());
//		System.out.println(new FloatingPointMass(0+57.021464f,17).getMass());
//		System.out.println(new FloatingPointMass(1000f, 17).getMass());
//		System.out.println(new FloatingPointMass(1000f+57.021464f,17).getMass());
//	}
	
//	public static void testGraphErrorDist() throws Exception
//	{
//		int maxLength = 25;
//		int numBits = 15;
//		int ppmInterested = 32;
//		float percentile = 0.99f;
//		
//		ScoreDistFactory sdFactory = new ScoreDistFactory(true, false);
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
//		FloatingPointMassFactory factory = new FloatingPointMassFactory(aaSet, numBits, maxLength+2);
//		ArrayList<FloatingPointMass> allNodes = factory.getAllNodes();
//		HashMap<FloatingPointMass,ScoreDist> errMap = new HashMap<FloatingPointMass,ScoreDist>(); 
//		ScoreDist initSD = sdFactory.getInstance(0, 1);
//		initSD.addNumber(0, 1);
//		
//		errMap.put(factory.getZero(), initSD);
//		
//		for(FloatingPointMass curNode : allNodes)
//		{
//			if(curNode == factory.getZero())
//				continue;
//			
//			int minPPMErr = Integer.MAX_VALUE;
//			int maxPPMErr = Integer.MIN_VALUE;
//			
//			ArrayList<FloatingPointMass> allPrevNodes = factory.getPreviousNodes(curNode);
//			ArrayList<FloatingPointMass> prevNodes = new ArrayList<FloatingPointMass>();
//			ArrayList<Float> errors = new ArrayList<Float>();
//			for(FloatingPointMass prevNode : allPrevNodes)
//			{
//				AminoAcid aa = aaSet.getAminoAcid(factory.getEdgeIndex(curNode, prevNode));
//				float error = (curNode.getMass()-prevNode.getMass())-aa.getMass();
//				
//				ScoreDist prevDist = errMap.get(prevNode);
//				if(prevDist == null)
//					continue;
//				prevNodes.add(prevNode);
//				errors.add(error);
//				float prevMass = prevNode.getMass();
//				float prevMinErr = 1e-6f*prevDist.getMinScore()*prevMass;
//				float prevMaxErr = 1e-6f*(prevDist.getMaxScore()-1)*prevMass;
//				
//				float curMinErr = prevMinErr + error;
//				int curMinErrPPM = Math.round(curMinErr/curNode.getMass()*1e6f);
//				
//				float curMaxErr = prevMaxErr + error;
//				int curMaxErrPPM = Math.round(curMaxErr/curNode.getMass()*1e6f);
//				
//				if(curMinErrPPM < minPPMErr)
//					minPPMErr = curMinErrPPM;
//				if(curMaxErrPPM > maxPPMErr)
//					maxPPMErr = curMaxErrPPM;
//			}			
//			assert(minPPMErr <= maxPPMErr);
//			
//			ScoreDist curDist = sdFactory.getInstance(minPPMErr, maxPPMErr+1);
//			int index = -1;
//			for(FloatingPointMass prevNode : prevNodes)
//			{
//				index++;
//				float error = errors.get(index);
//				ScoreDist prevDist = errMap.get(prevNode);
//				
//				for(int i=prevDist.getMinScore(); i<prevDist.getMaxScore(); i++)
//				{
//					float prevMass = prevNode.getMass();
//					float prevErr = 1e-6f*i*prevMass;
//					float curErr = prevErr + error;
//					int curErrPPM = Math.round(curErr/curNode.getMass()*1e6f);
//					curDist.addNumber(curErrPPM, prevDist.getNumberRecs(i));
//				}
//			}
//			double sum = 0;
//			double interested = 0;
//			for(int i=curDist.getMinScore(); i<curDist.getMaxScore(); i++)
//			{
//				sum += curDist.getNumberRecs(i);
//				if(i>=-ppmInterested && i<=ppmInterested)
//					interested += curDist.getNumberRecs(i);
//			}
//			assert(sum > 0): curNode.getMass();
//			double ratio = (interested / sum);
////			System.out.print(curDist.getMinScore()+","+curDist.getMaxScore()+"\t");
//			System.out.print(curNode.getMass()+","+ratio);//+"\t"+curDist.getMeanScore());
////			System.out.print("\t"+curDist.getMinScore()+","+curDist.getMaxScore());
//			
//			int mean = Math.round(curDist.getMeanScore());
//			double adjustedInterested = 0;
//			for(int score = mean-ppmInterested; score<=mean+ppmInterested; score++)
//			{
//				if(score >= curDist.getMinScore() && score<curDist.getMaxScore())
//					adjustedInterested += curDist.getNumberRecs(score);
//			}
//			double adjustedRatio = adjustedInterested/sum;
//			System.out.print(","+mean+","+adjustedRatio);
//			System.out.print(","+curDist.getMinScore()+","+(curDist.getMaxScore()-1));
//			System.out.println();
//			errMap.put(curNode, curDist);
//		}
//	}	
	
//	public static void testGraphTheoreticalErrors() throws Exception
//	{
//		int maxLength = 20;
//		int numBits = 15;
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		FloatingPointMassFactory factory = new FloatingPointMassFactory(aaSet, numBits, maxLength+2);
//		ArrayList<FloatingPointMass> allNodes = factory.getAllNodes();
//		HashMap<FloatingPointMass,Pair<Float,Float>> errMap = new HashMap<FloatingPointMass,Pair<Float,Float>>(); 
//		errMap.put(factory.getZero(), new Pair<Float,Float>(0f,0f));
//		
//		for(FloatingPointMass curNode : allNodes)
//		{
//			float errMin = Integer.MAX_VALUE;
//			float errMax = Integer.MIN_VALUE;
//			for(FloatingPointMass prevNode : factory.getPreviousNodes(curNode))
//			{
//				AminoAcid aa = aaSet.getAminoAcid(factory.getEdgeIndex(curNode, prevNode));
//				float curError = (curNode.getMass()-prevNode.getMass())-aa.getMass();
//				
//				Pair<Float,Float> prevErr = errMap.get(prevNode);
//				if(prevErr == null)
//					continue;
//				
//				float curErrMin = prevErr.getFirst()+curError;
//				float curErrMax = prevErr.getSecond()+curError;
//				if(curErrMin < errMin)
//					errMin = curErrMin;
//				if(curErrMax > errMax)
//					errMax = curErrMax;
//			}
//			if(curNode.getMass() > 0)
//			{
//				errMap.put(curNode, new Pair<Float,Float>(errMin, errMax));
//				float errMinPPM = errMin/curNode.getMass()*1e6f;
//				float errMaxPPM = errMax/curNode.getMass()*1e6f;
//				System.out.println(curNode.getMass()+","+errMinPPM+","+errMaxPPM);
//			}
//		}
//	}
	
//	public static void testGraphErrors() throws Exception
//	{
//		int maxLength = 20;
//		int numBits = 17;
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		FloatingPointMassFactory factory = new FloatingPointMassFactory(aaSet, numBits, maxLength+2, Enzyme.TRYPSIN);
//		String fileName = System.getProperty("user.home")+"/Research/Data/AgilentQTOF/annotatedAgilentQTOF.mgf";
//		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
//		int specNum = -1;
//		
//		HashMap<Integer,Histogram<Integer>> errSummary = new HashMap<Integer,Histogram<Integer>>();
//		
//		while(itr.hasNext())
//		{
//			Spectrum spec = itr.next();
//			specNum++;
//			Peptide pep = spec.getAnnotation();
//			if(pep.size() > maxLength)
//				continue;
//
//			for(int p=0; p<1; p++)
//			{
//				Sequence<FloatingPointMass> fpSeq = factory.toCumulativeSequence(p==0, pep);
//				Composition comp = Composition.NIL;
//				for(int i=0; i<pep.size(); i++)
//				{
//					Composition c;
//					if(p==0)
//						c = pep.get(i).getComposition();
//					else
//						c = pep.get(pep.size()-1-i).getComposition();
//					comp = comp.getAddition(c);
//					float compMass = comp.getMass();
//					float fpMass = fpSeq.get(i).getMass();
//					float errPPM = (fpMass-compMass)/compMass*1e6f;
//					int x = Math.round(compMass*Constants.INTEGER_MASS_SCALER);
//					int y = Math.round(errPPM);
//					Histogram<Integer> hist = errSummary.get(x);
//					if(hist == null)
//					{
//						hist = new Histogram<Integer>();
//						errSummary.put(x, hist);
//					}
//					hist.add(y);
//				}
//			}
//		}
//		ArrayList<Integer> massList = new ArrayList<Integer>(errSummary.keySet());
//		Collections.sort(massList);
//		for(int m : massList)
//		{
//			System.out.println("Mass: " + m);
//			Histogram<Integer> hist = errSummary.get(m);
//			hist.printSorted();
//		}
//		
//		System.out.print("mass_"+numBits+"=[");
//		for(int i=0; i<massList.size()-1; i++)
//			System.out.print(massList.get(i)+",");
//		System.out.println(massList.get(massList.size()-1)+"];");
//		
//		System.out.print("maxErr_"+numBits+"=[");
//		for(int i=0; i<massList.size()-1; i++)
//		{
//			Histogram<Integer> hist = errSummary.get(massList.get(i));
//			int maxErr = Math.max(Math.abs(hist.minKey()), Math.abs(hist.maxKey()));
//			System.out.print(maxErr+",");
//		}
//		Histogram<Integer> hist = errSummary.get(massList.get(massList.size()-1));
//		int maxErr = Math.max(Math.abs(hist.minKey()), Math.abs(hist.maxKey()));
//		System.out.println(maxErr+"];");
//		
//	}
	

	public static void testFixedLengthAAGraph() throws Exception 
	{
		int numBits = 15;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		float maxErrPPM = 0;
		float avgErrPPM = 0;
		int size = 0;
		
		for(AminoAcid aa : aaSet)
		{	
			Composition c = aa.getComposition();
			size++;
			float mass = c.getMass();
			String byteStr = floatToByteStr(c.getMass(), numBits);
			float approxMass = byteStrToFloat(byteStr);
			float errPPM = (approxMass-mass)/mass*1e6f;
			if(Math.abs(errPPM) > maxErrPPM)
				maxErrPPM = Math.abs(errPPM);
			avgErrPPM += errPPM*aa.getProbability();
			System.out.println(c+"\t"+c.getMass()+"\t"+byteStr+"\t"+approxMass+"\t"+errPPM+"\t"+round(mass,numBits));
		}
		System.out.println("MaxErr: " + maxErrPPM);
		System.out.println("AvgErr: " + avgErrPPM);
	}
	
	private static float round(float f, int numBits)
	{
		int intVal = Float.floatToIntBits(f);
		int carry = (intVal & (1 << (23-numBits-1))) == 0 ? 0 : 1;
		int intRounded = ((intVal >> (23-numBits))+ carry) << (23-numBits);
		float rounded = Float.intBitsToFloat(intRounded);
		return rounded;
	}
	
	private static float byteStrToFloat(String byteStr)
	{
		return Float.intBitsToFloat(Integer.parseInt(byteStr, 2));
	}
	
	private static String floatToByteStr(float f, int numBits)
	{
		int intValue = Float.floatToIntBits(f);
		String byteStr = Integer.toBinaryString(intValue);
		for(int i=byteStr.length(); i<32; i++)
			byteStr = "0" + byteStr;
			
		String sign = byteStr.substring(0, 1);
		String exp = byteStr.substring(1, 9);
		StringBuffer fraction = new StringBuffer(byteStr.substring(9));

		if(numBits < fraction.length())
		{
			char carry = fraction.charAt(numBits);
			for(int i=fraction.length()-1; i>=0; i--)
			{
				if(i >= numBits)
					fraction.setCharAt(i, '0');
				else if(carry == '1')
				{
					if(fraction.charAt(i) == '1')
						fraction.setCharAt(i, '0');
					else
					{
						fraction.setCharAt(i, '1');
						break;
					}
				}
			}
		}
		
		return sign+exp+fraction;
	}
	

	public static void graphTest() throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		int maxLength = 20;
		long time = System.currentTimeMillis();
		CompositionFactory allCompositions = new CompositionFactory(aaSet, null, maxLength);
		
		System.out.println("Composibion Building: " + (System.currentTimeMillis()-time));
		
		Tolerance tol = Tolerance.parseToleranceStr("30ppm");
		float maxJumpMass = 700;
		String fileName = "/home/sangtaekim/Research/Data/AgilentQTOF/annotatedAgilentQTOF.mgf";
		WindowFilter filter = new WindowFilter(6, 50);
		
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		int numSpecs = 0;
		int sumEdges = 0;
		while(itr.hasNext())
		{
			Spectrum s = itr.next();
			Spectrum spec = filter.apply(s);
			Peptide pep = spec.getAnnotation();
			
			if(pep.size() > maxLength || spec.getCharge() != 2)
				continue;
			
			numSpecs++;

			int numEdges = 0;
			// construct a graph
			for(int i=0; i<spec.size(); i++)
			{
				Peak p1 = spec.get(i);
				for(int j=i+1; j<spec.size(); j++)
				{
					Peak p2 = spec.get(j);
					float diff = p2.getMass() - p1.getMass();
					if(allCompositions.getNodes(diff, tol).size() > 0)
						numEdges++;
				}
			}
			sumEdges += numEdges;
			System.out.println(pep+"\t"+numEdges);
		}
		
		System.out.println("AvgNumEdges: " + sumEdges/(float)numSpecs);
		
	}
	
	public static void ionProb() throws Exception
	{
		IonType[] ions = {IonType.Y}; // IonType.getIonType("y2")
		Tolerance tol = Tolerance.parseToleranceStr("30ppm");
		int numJumps = 0;
		float maxJumpMass = 500;
		int maxLength = 20;
		
		String fileName = "/home/sangtaekim/Research/Data/AgilentQTOF/annotatedAgilentQTOF.mgf";
		WindowFilter filter = new WindowFilter(6, 50);
		
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		int numCleavages = 0;
		int numCleavagesWithPeaks = 0;
		int numSpecs = 0;
		int numPeaks = 0;
		int numSpecsWithSpecGraphPath = 0;
		while(itr.hasNext())
		{
			Spectrum s = itr.next();
//			Spectrum spec = filter.apply(s);
			Spectrum spec = s;
			Peptide pep = spec.getAnnotation();
			
			if(pep.size() > maxLength || spec.getCharge() != 2)
				continue;
			
			numSpecs++;
			numPeaks += spec.size();
			
			float suffixMass = 0;
			boolean specGraphHasPath = true;
			float prevSuffixMassWithPeak = 0;
			int prevI = pep.size();
			for(int i=pep.size()-1; i>0; i--)
			{
				numCleavages++;
				suffixMass += pep.get(i).getMass();
				float prefixMass = spec.getParentMass() - (float)Composition.H2O - suffixMass;
//				float prefixMass = pep.getMass(0, i);
				boolean peakExists = false;
				for(IonType ion : ions)
				{
					float mass;
					if(ion instanceof PrefixIon)
						mass = ion.getMz(prefixMass);
					else
						mass = ion.getMz(suffixMass);
					ArrayList<Peak> peakList = spec.getPeakListByMass(mass, tol);
					if(peakList.size() > 0)	// peak exists
					{
						peakExists = true;
						prevI = i;
						prevSuffixMassWithPeak = suffixMass;
						break;
					}
				}
				if(peakExists)
					numCleavagesWithPeaks++;
				else
				{
					if(prevI-i <= numJumps)
						numCleavagesWithPeaks++;
					if(suffixMass - prevSuffixMassWithPeak > maxJumpMass)
					{
						specGraphHasPath = false;
					}
				}
			}
			if(specGraphHasPath)
				numSpecsWithSpecGraphPath++;
		}
		
		System.out.print("Ion:"); 
		for(IonType ion : ions)
			System.out.print(" "+ ion.getName());
		System.out.println();
		System.out.println("NumSpec: " + numSpecs);
		System.out.println("RatioSpecWithSpecGraphPath: " + numSpecsWithSpecGraphPath/(float)numSpecs);
		System.out.println("AverageNumPeaks: " + numPeaks/(float)numSpecs);
		System.out.println("IonProb: " + numCleavagesWithPeaks/(float)numCleavages);
	}
	
	public static void correctChargeZero() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/AgilentQTOF/annotatedAgilentQTOF.mgf";
		
		String outputFileName = "/home/sangtaekim/Research/Data/AgilentQTOF/annotatedAgilentQTOF_NoC0.mgf";
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			if(spec.getCharge() == 0)
			{
				Peptide pep = spec.getAnnotation();
				int charge = Math.round(pep.getParentMass()/spec.getPrecursorPeak().getMz());
				spec.setCharge(charge);
			}
			spec.outputMgf(out);
		}
		out.close();
	}

	public static void compositionDensity() throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		int maxLength = 10;
		long time = System.currentTimeMillis();
		CompositionFactory allCompositions = new CompositionFactory(aaSet, null, maxLength);
		System.out.println("Composibion Building: " + (System.currentTimeMillis()-time));

		float clusterSizePPM = 30;
		int clusterIndex = 1;
		ArrayList<Composition> cluster = new ArrayList<Composition>();
		System.out.println("0\t0\t1");	// 0
		for(int i=1; i<allCompositions.getData().length; i++)
		{
			Composition comp = new Composition(allCompositions.getData()[i]);
			if(cluster.size() == 0)
			{
				cluster.add(comp);
			}
			else
			{
				float firstCompMass = cluster.get(0).getMass();
				float massDiffPPM = (comp.getMass() - firstCompMass)/firstCompMass*1e6f;
				if(massDiffPPM < clusterSizePPM)
					cluster.add(comp);
				else
				{
					float avgMass = (cluster.get(0).getMass()+cluster.get(cluster.size()-1).getMass())/2;
					System.out.println(clusterIndex+"\t"+avgMass+"\t"+cluster.size()+"\t"+i);
					clusterIndex++;
					cluster.clear();
				}
			}
		}
		
//		int counter = 0;
//		for(int c : allCompositions.getData())
//		{
//			counter++;
//			if(counter < 10 || counter < 100 && counter % 10 == 0 || counter < 1000 && counter % 100 == 0 || counter < 10000 && counter % 1000 == 0 || counter % 10000 == 0)
//			{
//				Composition comp = new Composition(c);
//				System.out.println(comp.getMass()+"\t"+allCompositions.getCompositions(comp.getMass(), tolerance).size());
//			}
//		}
//		Histogram<Integer> hist = new Histogram<Integer>();
//		for(int c : allCompositions.getData())
//		{
//			Composition comp = new Composition(c);
//			hist.add(comp.getNominalMass());
//		}
//		hist.printSorted();
	}
	

}
