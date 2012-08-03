package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import edu.ucsd.msjava.msgf.AminoAcidGraph;
import edu.ucsd.msjava.msgf.Histogram;
import edu.ucsd.msjava.msgf.IntMassFactory;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msgf.NominalMassFactory;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.IonProbability;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScoredSpectrum;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.CompositionFactory;
import edu.ucsd.msjava.msutil.Constants;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.IonType;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.SpectraContainer;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.WindowFilter;
import edu.ucsd.msjava.msutil.Modification.Location;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.suffixarray.SuffixArray;
import edu.ucsd.msjava.suffixarray.SuffixArraySequence;


public class Zubarev {
	public static void main(String argv[]) throws Exception
	{
//		splitHCDETD();
//		compositionMSGF();
//		addScanNum();
//		summarizeAnnotatedSpec();
//		testDBSearch();
//		countNTT();
//		isMascotTwoPassSearch2();
//		errorTest();
//		rescalingTest();
//		compareGraphSize();
//		countTryptic();
//		paramTest();
//		computeDBSize();
//		testEdgeScoresNominalMass();
//		testEdgeScoresRandom();
//		testRandomPeakProb();
//		testTagging();
//		filtrationPowerNominalMass();
//		filtrationPowerComposition(10);
		deconvolutionTest();
	}

	public static void deconvolutionTest() throws Exception
	{
		for(int charge=3; charge<=3; charge++)
		{
			System.out.println("Charge: " + charge);
			System.out.println("Before deconvolution");
			deconvolutionTest(charge, false);
			System.out.println("After deconvolution");
			deconvolutionTest(charge, true);
			System.out.println();
		}
	}
	
	public static void deconvolutionTest(int charge, boolean applyDeconvolution) throws Exception
	{
		String specFileName = System.getProperty("user.home")+"/Research/Data/Heck_DDDT/AnnotatedSpectra/HCD_HighRes_Tryp.mgf";
//		specFileName = System.getProperty("user.home")+"/Research/Data/Zubarev/AnnotatedSpectra/Zubarev_HCD_Annotated.mgf";
		SpectraContainer container = new SpectraContainer();

		SpectraIterator itr = new SpectraIterator(specFileName, new MgfSpectrumParser());
		
//		IonType[] ionTypes = {IonType.B, IonType.getIonType("b+n"), IonType.Y, IonType.getIonType("y+n"), IonType.A, IonType.getIonType("b2"), IonType.getIonType("y2"), IonType.getIonType("b3"), IonType.getIonType("y3")};
		IonType[] ionTypes = IonType.getAllKnownIonTypes(3, true, false).toArray(new IonType[0]);
		float toleranceBetweenIsotopes = 0.02f;
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			if(spec.getCharge() != charge)
				continue;
			if(applyDeconvolution)
				container.add(spec.getDeconvolutedSpectrum(toleranceBetweenIsotopes));
			else
				container.add(spec);
		}
		
		IonProbability probGen = new IonProbability(container.iterator(), ionTypes, new Tolerance(30, true));
		float[] prob = probGen.getIonProb();
		int i=0;
		for(IonType ion : ionTypes)
		{
			if(prob[i] > 0.05)
				System.out.println(ion.getName()+"\t"+prob[i]);
			i++;
		}
	}

	public static void filtrationPowerComposition(float tolerancePPM) throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		int maxMass = 700;
		float[] filtFactorNM = new float[maxMass+1];
		filtFactorNM[0] = 1;
		for(int mass=1; mass<filtFactorNM.length; mass++)
		{
			for(AminoAcid aa : aaSet.getAAList(Location.Anywhere))
			{
				int prevMass = mass - aa.getNominalMass();
				if(prevMass < 0)
					continue;
				filtFactorNM[mass] += filtFactorNM[prevMass]*0.05f;
			}
		}
		
		int maxLength = 13;
		Tolerance tolerance = new Tolerance(tolerancePPM, true);
		CompositionFactory factory = new CompositionFactory(aaSet, null, maxLength);
		HashMap<Composition,Float> filtFactor = new HashMap<Composition,Float>(); 
		filtFactor.put(Composition.NIL, 1f);
		for(int compNum : factory.getData())
		{
			Composition c = new Composition(compNum);
			if(c.equals(Composition.NIL))
				continue;
			float prob = 0;
			for(AminoAcid aa : aaSet.getAAList(Location.Anywhere))
			{
				Composition aaComp = aa.getComposition();
				Composition prevComp = c.getSubtraction(aaComp);
				if(prevComp == null)
					continue;
				Float prevProb = filtFactor.get(prevComp);
				if(prevProb == null)
					continue;
				prob += prevProb*0.05f;
			}
			filtFactor.put(c, prob);
		}
		System.out.println("Mass\tFiltrationFactor");
		for(int compNum : factory.getData())
		{
			Composition c = new Composition(compNum);
			if(c.getMass() > 700)
				break;
			float prob = 0;
			for(Composition neighbor : factory.getNodes(c.getMass(), tolerance))
			{
				prob += filtFactor.get(neighbor);
			}
			System.out.println(c.getMass()+"\t"+filtFactorNM[c.getNominalMass()]+"\t"+prob);
		}
	}
	
	public static void filtrationPowerNominalMass() throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		int maxMass = 700;
		float[] filtFactor = new float[maxMass+1];
		filtFactor[0] = 1;
		for(int mass=1; mass<filtFactor.length; mass++)
		{
			for(AminoAcid aa : aaSet.getAAList(Location.Anywhere))
			{
				int prevMass = mass - aa.getNominalMass();
				if(prevMass < 0)
					continue;
				filtFactor[mass] += filtFactor[prevMass]*0.05f;
			}
		}
		System.out.println("Mass\tFiltrationFactor");
		for(int mass=0; mass<filtFactor.length; mass++)
			System.out.println(NominalMass.getMassFromNominalMass(mass)+"\t"+filtFactor[mass]);
	}
	
	public static void testRandomPeakProb() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/Zubarev/AnnotatedSpectra/Zubarev_HCD_Annotated.mgf";
		fileName = System.getProperty("user.home")+"/Research/Data/HeckRevision/AnnotatedSpectra/CID_Tryp_Confident.mgf";

		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		float sumPeakRatio = 0;
		int numSpecs = 0;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			spec = new WindowFilter(10,50).apply(spec);
			numSpecs++;
			boolean[] booleanSpec = new boolean[spec.getAnnotation().getNominalMass()+1];
			booleanSpec[0] = true;
			booleanSpec[booleanSpec.length-1] = true;
			int numPeaks = 0;
			for(Peak p : spec)
			{
				int mz = NominalMass.toNominalMass(p.getMz());
				if(mz >= 0 && mz < booleanSpec.length)
				{
					if(booleanSpec[mz] == false)
					{
						numPeaks++;
						booleanSpec[mz] = true;
					}
				}
			}
			int numPairs = 0;
			int[] numReal = new int[4];	// yy: 3, yn: 2, ny: 1, nn:0
			for(int i=1; i<booleanSpec.length; i++)
			{
				for(AminoAcid aa : aaSet)
				{
					int aaMass = aa.getNominalMass();
					int prevMass = i-aaMass;
					if(prevMass >= 0)
					{
						numPairs++;
						if(booleanSpec[i] == true)
						{
							if(booleanSpec[prevMass] == true)
								numReal[3]++;
							else
								numReal[1]++;
						}
						else
						{
							if(booleanSpec[prevMass] == true)
								numReal[2]++;
							else
								numReal[0]++;
						}
					}
				}
			}
			
			System.out.println(spec.size()+" "+numPeaks+" "+spec.getParentMass());
			float peakRatio = numPeaks/(float)booleanSpec.length;
			System.out.print(numReal[0]/(float)numPairs+","+(1-peakRatio)*(1-peakRatio));
			System.out.print("\t"+numReal[1]/(float)numPairs+","+(1-peakRatio)*peakRatio);
			System.out.print("\t"+numReal[2]/(float)numPairs+","+peakRatio*(1-peakRatio));
			System.out.print("\t"+numReal[3]/(float)numPairs+","+peakRatio*peakRatio);
			System.out.println();
			int sum = numReal[0]+numReal[1]+numReal[2]+numReal[3]; 
			assert(sum == numPairs): sum+"!="+numPairs;
			sumPeakRatio += peakRatio;
//			System.out.println(peakRatio);
		}		
		System.out.println("Average\t"+sumPeakRatio/numSpecs);
	}
	
	public static void testEdgeScoresRandom() throws Exception
	{
		float[] nominalMass = new float[200];
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		for(AminoAcid aa : aaSet)
			nominalMass[aa.getNominalMass()] = aa.getMass();
		
		Histogram<Integer> errHist = new Histogram<Integer>();
		
		AminoAcid aaQ = aaSet.getAminoAcid('Q');
		AminoAcid aaK = aaSet.getAminoAcid('K');
		String fileName = System.getProperty("user.home")+"/Research/Data/Zubarev/AnnotatedSpectra/Zubarev_HCD_Annotated.mgf";
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		int numPairs = 0;
		HashSet<String> pepSet = new HashSet<String>();

		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
//			spec = new WindowFilter(6,50).apply(spec);
			Peptide peptide = spec.getAnnotation();
			
			if(pepSet.contains(peptide.toString()))
				continue;
			else
				pepSet.add(peptide.toString());
			
			boolean[] nominalBY = new boolean[10000];
			int nominalB=1, nominalY=19;
			for(int i=0; i<peptide.size()-1; i++)
			{
				nominalB += peptide.get(i).getNominalMass();
				nominalY += peptide.get(peptide.size()-1-i).getNominalMass();
				nominalBY[nominalB] = true;
				nominalBY[nominalY] = true;
			}
			for(int i=0; i<spec.size()-1; i++)
			{
				Peak p1 = spec.get(i);
				float p1Mass = p1.getMz();
				int nominalP1 = NominalMass.toNominalMass(p1.getMz());
				if(nominalBY[nominalP1] == true)
					continue;
				for(int j=i+1; j<spec.size(); j++)
				{
					Peak p2 = spec.get(j);
					float p2Mass = p2.getMz();
					int nominalP2 = NominalMass.toNominalMass(p2.getMz());
					if(nominalBY[nominalP2] == true)
						continue;
					int nominalDiff = nominalP2-nominalP1;
					assert(nominalDiff >= 0);
					if(nominalDiff > 186)
						break;
					if(nominalMass[nominalDiff] == 0)
						continue;
					numPairs++;
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
//					System.out.println(err);
					errHist.add(Math.round(err*100));
				}
			}
		}
		System.out.println("NumPairs: " + numPairs);
		errHist.printSortedRatio();
	}
	
	public static void testEdgeScoresNominalMass() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/Zubarev/AnnotatedSpectra/Zubarev_HCD_Annotated.mgf";
//		fileName = System.getProperty("user.home")+"/Research/Data/HeckRevision/AnnotatedSpectra/CID_Tryp_Confident.mgf";
		int scalingFactor = 100;
//		scalingFactor = 10;
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		int numYY = 0;
		int numNY = 0;
		int numYN = 0;
		int numNN = 0;
		
		int numY=0;
		int numB=0;
		int numBY = 0;
		int numNo = 0;

		int numSpecs = 0;
		int numCorrectLength3Tags = 0;
		
		float rescalingConstant = Constants.INTEGER_MASS_SCALER;
//		rescalingConstant = Constants.INTEGER_MASS_SCALER_HIGH_PRECISION;
		
		Histogram<Integer> yHist = new Histogram<Integer>();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
//			spec = new WindowFilter(6, 50).apply(spec);
//			if(spec.getCharge() !=3)
//				continue;
			Peptide peptide = spec.getAnnotation();
			if(Math.abs(spec.getParentMass()-(float)Composition.H2O-peptide.getMass()) > 0.5f)
				continue;
			
			numSpecs++;
//			Collections.reverse(peptide);
//			Collections.shuffle(peptide);
			
			int intSRM = 0;
			int intPeptideMass = 0;
			for(int i=0; i<peptide.size(); i++)
				intPeptideMass += Math.round(peptide.get(i).getMass()*rescalingConstant);
			float[] b = new float[peptide.size()+1];
			float[] y = new float[peptide.size()+1];
			
			b[0] = IonType.B.getOffset(); 
			y[0] = IonType.Y.getOffset();
			for(int i=0; i<peptide.size()-1; i++)
			{
				int byPresence = 0;
				intSRM += Math.round(peptide.get(peptide.size()-1-i).getMass()*rescalingConstant);
				float yMass = IonType.Y.getMz(intSRM/rescalingConstant);
				Peak p = spec.getPeakByMass(yMass, new Tolerance(0.5f));
				if(p != null)
				{
					y[i+1] = p.getMz();
					byPresence += 2;
//					y[i+1] = -1;
				}
				else
				{
//					float y2Mass = IonType.getIonType("y2").getMz(intSRM/rescalingConstant);
//					Peak p2 = spec.getPeakByMass(y2Mass, new Tolerance(0.5f));
//					if(p2 != null)
//						y[i+1] = p2.getMz()*2-(float)Composition.H;
//					else
						y[i+1] = -1;
				}
				
				int intPRM = intPeptideMass - intSRM;
				float bMass = IonType.B.getMz(intPRM/rescalingConstant);
				p = spec.getPeakByMass(bMass, new Tolerance(0.5f)); 
				if(p != null)
				{
					b[peptide.size()-(i+1)] = p.getMz();
					byPresence += 1;
				}
				else
				{
					b[peptide.size()-(i+1)] = -1;
				}
				
				if(byPresence == 1)	// no y and yes b
					y[i+1] = spec.getParentMass() + (float)Composition.H2 - b[peptide.size()-(i+1)];
//				else if(byPresence == 2)	// yes y and no b
//					b[peptide.size()-(i+1)] = spec.getParentMass() + (float)Composition.H2 - y[i+1];
					
				if(byPresence == 3)
					numBY++;
				else if(byPresence == 2)
					numY++;
				else if(byPresence == 1)
					numB++;
				else
					numNo++;
			}
			
			intSRM += Math.round(peptide.get(0).getMass()*Constants.INTEGER_MASS_SCALER);
			assert(intSRM == intPeptideMass): peptide+": " + intSRM + " != " + intPeptideMass + " " + (spec.getParentMass()-(float)Composition.H2O);
			b[peptide.size()] = spec.getParentMass()-(float)Composition.H2O+(float)Composition.H;
			y[peptide.size()] = spec.getParentMass()+(float)Composition.H;
			
			boolean useY = true;
			float[] mainIon;
			if(useY)
				mainIon = y;
			else
				mainIon = b;
			
			for(int i=1; i<=peptide.size(); i++)
			{
				if(mainIon[i] >= 0)
				{
					if(mainIon[i-1] >=0)
					{
						numYY++;
						AminoAcid aa;
						if(useY)
							aa = peptide.get(peptide.size()-i);
						else
							aa = peptide.get(i-1);
						float expMass = mainIon[i]-mainIon[i-1];
						float theoMass = aa.getMass();
//						if(y[i] > 500f)
//							continue;
						float diff = expMass-theoMass;
//						if(diff < 0.1f && diff > -0.1f)
//							System.out.println(diff);
						yHist.add(Math.round(diff*scalingFactor));
//						System.out.println("YY");
					}
					else
					{
						numNY++;
//						System.out.println("NY");
					}
				}
				else
				{
					if(mainIon[i-1] >=0)
					{
//						System.out.println("YN");
						numYN++;
					}
					else
					{
//						System.out.println("NN");
						numNN++;
					}
				}
			}
			
			// length3 tags
			boolean tagCorrect = false;
			for(int i=3; i<=peptide.size(); i++)
			{
				if(mainIon[i] >= 0)
				{
					for(int j=1; j<=3; j++)
					if(mainIon[i-j] >=0)
					{
						AminoAcid aa;
						if(useY)
							aa = peptide.get(peptide.size()-1-(i-j));
						else
							aa = peptide.get(i-j);
						float expMass = mainIon[i-j+1]-mainIon[i-j];
						float theoMass = aa.getMass();
						float diff = expMass-theoMass;
//						if(diff < 5.f/scalingFactor)
						{
							tagCorrect = true;
							break;
						}
					}
					else
						break;
				}
				else
				{
					break;
				}
			}			
			if(tagCorrect)
				numCorrectLength3Tags++;
		}
		int sum = numYY+numYN+numNY+numNN;
		System.out.println("NumYY: " + numYY+" "+numYY/(float)sum);
		System.out.println("NumYN: " + numYN+" "+numYN/(float)sum);
		System.out.println("NumNY: " + numNY+" "+numNY/(float)sum);
		System.out.println("NumNN: " + numNN+" "+numNN/(float)sum);
		System.out.println();
		int sum2 = numBY+numY+numB+numNo;
		System.out.println("NumBothBY: " + numBY+ " " + numBY/(float)sum2);
		System.out.println("NumYOnly: " + numY+ " " + numY/(float)sum2);
		System.out.println("NumBOnly: " + numB+ " " + numB/(float)sum2);
		System.out.println("NumNoBY: " + numNo+ " " + numNo/(float)sum2);
		yHist.printSortedRatio();
		System.out.println("CorrectTag: " + numCorrectLength3Tags/(float)numSpecs);
	}	
	public static void testEdgeScores() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/Zubarev/AnnotatedSpectra/Zubarev_HCD_Annotated.mgf";
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		int numYY = 0;
		int numNY = 0;
		int numYN = 0;
		int numNN = 0;
		
		Histogram<Integer> yHist = new Histogram<Integer>();
		Histogram<Integer> yHistPPM = new Histogram<Integer>();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
//			spec = new WindowFilter(6, 50).apply(spec);
			spec.setRanksOfPeaks();
			Peptide peptide = spec.getAnnotation();
//			Collections.reverse(peptide);
			double prm = 0;
			double srm = 0;
			float[] b = new float[peptide.size()+1];
			float[] y = new float[peptide.size()+1];
			
			b[0] = IonType.B.getOffset(); 
			y[0] = IonType.Y.getOffset();
			for(int i=0; i<peptide.size()-1; i++)
			{
				prm += peptide.get(i).getAccurateMass();
				float bMass = IonType.B.getMz((float)prm);
				Peak p = spec.getPeakByMass(bMass, new Tolerance(0.5f)); 
				if(p != null)
					b[i+1] = p.getMz();
				else
					b[i+1] = -1;
				
				srm += peptide.get(peptide.size()-1-i).getAccurateMass();
				float yMass = IonType.Y.getMz((float)srm);
				p = spec.getPeakByMass(yMass, new Tolerance(0.5f));
				if(p != null)
					y[i+1] = p.getMz();
				else
					y[i+1] = -1;
			}
			srm += peptide.get(0).getAccurateMass();
			
			b[peptide.size()] = IonType.B.getMz((float)srm);
			y[peptide.size()] = IonType.Y.getMz((float)srm);
			
			boolean useY = false;
			float[] mainIon;
			if(useY)
				mainIon = y;
			else
				mainIon = b;
			
			for(int i=1; i<peptide.size(); i++)
			{
				if(mainIon[i] >= 0)
				{
					if(mainIon[i-1] >=0)
					{
						numYY++;
						AminoAcid aa;
						if(useY)
							aa = peptide.get(peptide.size()-i);
						else
							aa = peptide.get(i-1);
						float expMass = mainIon[i]-mainIon[i-1];
						float theoMass = aa.getMass();
//						if(y[i] > 500f)
//							continue;
						float diff = expMass-theoMass;
						float diffPPM = diff*1e6f/((mainIon[i]+mainIon[i-1])/2);
						if(Math.abs(diff) > 0.4f)
							continue;
						yHist.add(Math.round(diff*100));
						yHistPPM.add(Math.round(diffPPM/10));
					}
					else
						numNY++;
				}
				else
				{
					if(mainIon[i-1] >=0)
						numYN++;
					else
						numNN++;
				}
			}
		}
		System.out.println("NumYY: " + numYY);
		System.out.println("NumYN: " + numYN);
		System.out.println("NumNY: " + numNY);
		System.out.println("NumNN: " + numNN);
		yHist.printSortedRatio();
//		yHistPPM.printSortedRatio();
	}
	
//	public static void computeDBSize() throws Exception
//	{
//		String dbFileName = System.getProperty("user.home")+"/Research/Data/IPI/IPI_human_3.79.fasta";
//		dbFileName = System.getProperty("user.home")+"/Research/Data/Arabidopsis/TAIR9_pep_20090619.fasta";
//		dbFileName = System.getProperty("user.home")+"/Research/Data/YeastDB/Yeast_60.fasta";
//		SuffixArray sa = new SuffixArray(new SuffixArraySequence(dbFileName));
//		for(int length=5; length<=50; length++)
//			System.out.println(length+"\t"+sa.getNumDistinctSeq(length));
//	}
	
	public static void paramTest() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/Zubarev/SACTest/SACTest_Decoy.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		in.readLine();
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			String annotation = token[5];
			String pep = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
			float pepMass = aaSet.getPeptide(pep).getMass();
			float specPepMass = (Float.parseFloat(token[2])-(float)Composition.H)*Integer.parseInt(token[4]) - (float)Composition.H2O;
			if(Math.abs(specPepMass-pepMass) > 0.5)
			{
				System.out.println(pepMass+" != " + specPepMass);
			}
		}
	}
	
	public static void countTryptic() throws Exception
	{
		Enzyme trypsin = Enzyme.TRYPSIN;
		String fileName = System.getProperty("user.home")+"/Research/Data/Zubarev/MSGFDB_HCD_Target_MSGFNM.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		float threshold = 8.363478e-10f;
		threshold = 6.115349e-10f;
		in.readLine();
		String s;
		int numSpec = 0;
		int numTrypN = 0;
		int numTrypC = 0;
		int numNonTryptic = 0;
		int numSemiTryptic = 0;
		int numFullTryptic = 0;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 10)
				continue;
			float specProb = Float.parseFloat(token[12]);
			if(specProb > threshold)
				continue;
			numSpec++;
			String annotation = token[5];
			char nTermChar = annotation.charAt(0);
			boolean isNTermTryp = false;
			if(nTermChar == 'K' || nTermChar == 'R' || nTermChar == '_')
			{
				numTrypN++;
				isNTermTryp = true;
			}
			boolean isCTermTryp = false;
			char cTermChar = annotation.charAt(annotation.lastIndexOf('.')-1);
			if(cTermChar == 'K' || cTermChar == 'R')
			{
				numTrypC++;
				isCTermTryp = true;
			}
			if(isNTermTryp && isCTermTryp)
				numFullTryptic++;
			else if(isNTermTryp || isCTermTryp)
				numSemiTryptic++;
			else
				numNonTryptic++;
		}
		System.out.println("NumSpec\t"+numSpec);
		System.out.println("NTerm\t"+numTrypN/(float)numSpec);
		System.out.println("CTerm\t"+numTrypC/(float)numSpec);
		System.out.println("FullTryp\t"+numFullTryptic/(float)numSpec);
		System.out.println("SemiTryp\t"+numSemiTryptic/(float)numSpec);
		System.out.println("NonTryp\t"+numNonTryptic/(float)numSpec);
	}
	
	public static void compareGraphSize() throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		float rescalingFactor = 274.335215f;
		Tolerance pmTolerance = new Tolerance(30, true);
		Enzyme enzyme = Enzyme.TRYPSIN;
		int maxLength = 20;
		float numPeptides = 1;
		
		System.out.println("Length\tIntMass\tComposition");
		for(int length=20; length<=20; length++)
		{
			NominalMassFactory intMassFactoryNominalMass = new NominalMassFactory(aaSet, enzyme, length);
			IntMassFactory intMassFactory = new IntMassFactory(aaSet, enzyme, length, rescalingFactor);
			CompositionFactory compFactory = new CompositionFactory(aaSet, enzyme, length);
			numPeptides += Math.pow(20, length);
			System.out.println(length+"\t"+intMassFactoryNominalMass.size()+"\t"+intMassFactory.size()+"\t"+compFactory.size()+"\t"+numPeptides);
		}
	}
	
	public static void rescalingTest() throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		for(AminoAcid aa : aaSet)
		{
			double rescaled =  aa.getAccurateMass()*Constants.INTEGER_MASS_SCALER;
			System.out.println(aa.getResidueStr()+"\t"+aa.getAccurateMass()+"\t"+rescaled+"\t"+(aa.getNominalMass()-rescaled));
		}
	}
	
	public static void errorTest() throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		float constant = 260f;
		float bestConstant = constant;
		float bestMaxError = Float.MAX_VALUE; 
//		while(constant < 270f)
//		{
//			float maxErr = Float.MIN_VALUE;
//			for(AminoAcid aa : aaSet)
//			{
//				int nominalMass = Math.round(aa.getMass()*constant);
//				float mass = aa.getMass();
//				float err = mass-nominalMass/constant;
//				float ppmErr = err/mass*1e6f;
//				if(err > Math.abs(ppmErr))
//					maxErr = ppmErr;
////				System.out.println(aa.getResidue()+" "+mass+" "+nominalMass+" "+err+" "+ppmErr);
//			}
//			if(maxErr < bestMaxError)
//				bestConstant = constant;
//			constant += 0.1f;
//		}
		bestConstant = 274.335215f;
		System.out.println(bestConstant);
		for(AminoAcid aa : aaSet)
		{
			int nominalMass = Math.round(aa.getMass()*bestConstant);
			float mass = aa.getMass();
			float err = mass-nominalMass/bestConstant;
			float ppmErr = err/mass*1e6f;
			System.out.println(aa.getResidueStr()+" "+mass+" "+nominalMass/bestConstant+" "+err+" "+ppmErr+" "+nominalMass);
		}
	}
	
	public static void isMascotTwoPassSearch2() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/Zubarev/Mascot_HCD_PTM/Mascot_HCD_PTM_Target.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		in.readLine();
		String s;
		HashSet<String> allProt = new HashSet<String>();
		HashSet<String> protNoPTM = new HashSet<String>();
		HashSet<String> protWithPTM = new HashSet<String>();
		int numNoPTM = 0;
		int numWithPTM = 0;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			String annotation = token[2];
			String protein = token[4].substring(0, token[4].indexOf(':'));
			int sIndex = token[4].indexOf(':', token[4].indexOf(':')+1);
			int eIndex = token[4].indexOf(':',sIndex+1);
			int matchedPosition = Integer.parseInt(token[4].substring(sIndex+1, eIndex));
			
			allProt.add(protein);
			String peptide = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
			if(peptide.equals(peptide.toUpperCase()))
			{
				numNoPTM++;
				protNoPTM.add(protein);
			}
			else
			{
				numWithPTM++;
				protWithPTM.add(protein);
			}
		}
		
		int numProtWithPTM = 0;
		for(String prot : protWithPTM)
		{
			if(!protNoPTM.contains(prot))
			{
				numProtWithPTM++;
				System.out.println(prot);
			}
		}
		
		int numProtWithNoPTM = 0;
		for(String prot : protNoPTM)
		{
			if(!protWithPTM.contains(prot))
				numProtWithNoPTM++;
		}
		
		System.out.println("NumProt: " + allProt.size());
		System.out.println("NumProtWithExclusivePTM: " + numProtWithPTM);
		System.out.println("NumProtWithExclusiveNoPTM: " + numProtWithNoPTM);
		System.out.println("NumNoPTM: " + numNoPTM);
		System.out.println("NumPTM: " + numWithPTM);
	}
	
	public static void isMascotTwoPassSearch() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/Zubarev/Mascot/Zubarev_HCD_Target.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		in.readLine();
		String s;
		HashSet<String> allProt = new HashSet<String>();
		HashSet<String> protTryptic = new HashSet<String>();
		HashSet<String> protSemiTryptic = new HashSet<String>();
		int numTryptic = 0;
		int numSemiTryptic = 0;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			String annotation = token[2];
			String protein = token[4].substring(0, token[4].indexOf(':'));
			int sIndex = token[4].indexOf(':', token[4].indexOf(':')+1);
			int eIndex = token[4].indexOf(':',sIndex+1);
			int matchedPosition = Integer.parseInt(token[4].substring(sIndex+1, eIndex));
			
			allProt.add(protein);
			boolean isNTermTryptic = false;
			char nTerm = annotation.charAt(0);
			if(nTerm == 'K' || nTerm == 'R' || nTerm =='-' || (nTerm == 'M' && matchedPosition==2))
				isNTermTryptic = true;
			char cTerm = annotation.charAt(annotation.length()-3);
			boolean isCTermTryptic = false;
			if(cTerm == 'K' || cTerm == 'R' || annotation.charAt(annotation.length()-1) == '-')
				isCTermTryptic = true;
			
			if(isNTermTryptic && isCTermTryptic)
			{
				protTryptic.add(protein);
				numTryptic++;
			}
			else
			{
				protSemiTryptic.add(protein);
				numSemiTryptic++;
			}
		}
		
		int numProtWithSemiTryptic = 0;
		for(String prot : protSemiTryptic)
		{
			if(!protTryptic.contains(prot))
			{
				numProtWithSemiTryptic++;
				System.out.println(prot);
			}
		}
		
		int numProtWithTryptic = 0;
		for(String prot : protTryptic)
		{
			if(!protSemiTryptic.contains(prot))
				numProtWithTryptic++;
		}
		
		System.out.println("NumProt: " + allProt.size());
		System.out.println("NumProtWithExclusiveSemiTryptic: " + numProtWithSemiTryptic);
		System.out.println("NumProtWithExclusiveFullTryptic: " + numProtWithTryptic);
		System.out.println("NumTryptic: " + numTryptic);
		System.out.println("NumSemiTryptic: " + numSemiTryptic);
	}
	
	public static void countNTT() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/Zubarev/Mascot/Zubarev_HCD_Target.txt";
		String s;
		BufferedLineReader in = new BufferedLineReader(fileName);
		in.readLine();
		int numSpec = 0;
		int numFT = 0;	// [KR].~~[KR]
		int numSTN = 0;	// [^KR].~~[KR]
		int numSTC = 0;	// [KR].~~[^KR]
		int numNT = 0;	// [^KR].~~[^KR]
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 5)
				continue;
			String annotation = token[2];
			boolean isNTermTryptic = false;
			char nTerm = annotation.charAt(0);
			if(nTerm == 'K' || nTerm == 'R' || nTerm =='-')
				isNTermTryptic = true;
			char cTerm = annotation.charAt(annotation.length()-3);
			boolean isCTermTryptic = false;
			if(cTerm == 'K' || cTerm == 'R')
				isCTermTryptic = true;
			
			numSpec++;
			if(isNTermTryptic && isCTermTryptic)
				numFT++;
			else if(!isNTermTryptic && isCTermTryptic)
				numSTN++;
			else if(isNTermTryptic && !isCTermTryptic)
				numSTC++;
			else
				numNT++;
		}
		System.out.println("NumFT\t"+numFT/(float)numSpec);
		System.out.println("NumSTN\t"+numSTN/(float)numSpec);
		System.out.println("NumSTC\t"+numSTC/(float)numSpec);
		System.out.println("NumNT\t"+numNT/(float)numSpec);
	}
	
	public static void testDBSearch() throws Exception
	{
		float precursorMass = 413.3487f;
		int charge = 2;
		float peptideMass = (precursorMass-(float)Composition.PROTON)*charge - (float)Composition.H2O;
		File databaseFile = new File("/home/sangtaekim/Research/Data/IPI/IPI_human_3.79.fasta");
		SuffixArray sa = new SuffixArray(new SuffixArraySequence(databaseFile.getPath()));
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		Tolerance tolerance = new Tolerance(30, true);
		System.out.println(sa.getNumCandidatePeptides(aaSet, peptideMass, tolerance));
	}
	
	public static void summarizeAnnotatedSpec() throws Exception
	{
		String specFileName = System.getProperty("user.home")+"/Research/Data/Zubarev/AnnotatedSpectra/Zubarev_HCD_Annotated.mgf";
		SpectraIterator itr = new SpectraIterator(specFileName, new MgfSpectrumParser());
		Histogram<Integer> massDiffHist = new Histogram<Integer>();
		Histogram<Integer> fmDiffHist = new Histogram<Integer>();
		int numC12 = 0;
		int numC13 = 0;
		int numC14 = 0;
		int numSpec = 0;
		int numBIons = 0;
		int numYIons = 0;
		int numObservedBIons = 0;
		int numObservedBIonsBySpecMass = 0;
		int numObservedYIons = 0;
		
		int numObservedYIonsIntMass = 0;
		int numObservedBIonsIntMass = 0;
		
		int numObserbedBIonsIntMassIntPM = 0;
		
		Tolerance fragTolerance = new Tolerance(30, true);
//		Tolerance fragTolerance = new Tolerance(0.5f);
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		Enzyme enzyme = Enzyme.TRYPSIN;
		int maxLength = 50;
		IntMassFactory factory = new IntMassFactory(aaSet, enzyme, maxLength, (float)(274.335215));
		int numFT = 0;	// fully tryptic
		int numNT = 0;	// non-tryptic
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			numSpec++;
			Peptide pep = spec.getAnnotation();
			if(enzyme.isCleaved(pep))
				numFT++;
			else
				numNT++;
			
			float theoPepMass = pep.getMass();
			float specPepMass = spec.getParentMass() - (float)Composition.H2O;
			int diff = Math.round(specPepMass - theoPepMass);
			if(diff == 1)
			{
				specPepMass -= (float)Composition.ISOTOPE;
				numC13++;
			}
			else if(diff == 2)
			{
				specPepMass -= (float)(Composition.ISOTOPE*2);
				numC14++;
			}
			else
				numC12++;
			float massDiffPPM = (specPepMass-theoPepMass)/theoPepMass*1e6f;
			int intMassDiffPPM = Math.round(massDiffPPM);
			massDiffHist.add(intMassDiffPPM);

			Spectrum filteredSpec = new WindowFilter(6,50).apply(spec);
			double prm = 0;
			double srm = 0;
			int srmIndex = 0;
			int pmIndex = factory.getMassIndex(spec.getParentMass()-(float)Composition.H2O);
			
			for(int i=pep.size()-1; i>=0; i--)
			{
				srm += pep.get(i).getAccurateMass();
				prm += pep.get(pep.size()-1-i).getAccurateMass();
				double prmByPM = specPepMass - srm;
				
				srmIndex += factory.getMassIndex(pep.get(i).getMass());
				float srm2 = factory.getMassFromIndex(srmIndex);
				float prm2 = specPepMass - srm2;
				float prm3 = factory.getMassFromIndex(pmIndex - srmIndex);
				
				numYIons++;
				if(filteredSpec.getPeakByMass((float)(srm+Composition.OFFSET_Y), fragTolerance) != null)
					numObservedYIons++;
				
				if(filteredSpec.getPeakByMass((float)(srm2+Composition.OFFSET_Y), fragTolerance) != null)
					numObservedYIonsIntMass++;
				
				numBIons++;
				if(filteredSpec.getPeakByMass((float)(prm+Composition.OFFSET_B), fragTolerance) != null)
					numObservedBIons++;
				if(filteredSpec.getPeakByMass((float)(prmByPM+Composition.OFFSET_B), fragTolerance) != null)
					numObservedBIonsBySpecMass++;
				if(filteredSpec.getPeakByMass((float)(prm2+Composition.OFFSET_B), fragTolerance) != null)
					numObservedBIonsIntMass++;
				if(filteredSpec.getPeakByMass((float)(prm3+Composition.OFFSET_B), fragTolerance) != null)
					numObserbedBIonsIntMassIntPM++;
				
			}
		}
//		massDiffHist.printSortedRatio();
//		System.out.println("#C12\t"+numC12);
//		System.out.println("#C13\t"+numC13);
//		System.out.println("#C14\t"+numC14);
		System.out.println("NumSpec: " + numSpec);
		fmDiffHist.printSortedRatio();
		System.out.println("Y: " + numObservedYIons/(float)numYIons);
		System.out.println("B: " + numObservedBIons/(float)numBIons);
		System.out.println("B_SpecPM: " + numObservedBIonsBySpecMass/(float)numBIons);
		System.out.println("Y_Int: " + numObservedYIonsIntMass/(float)numYIons);
		System.out.println("B_Int: " + numObservedBIonsIntMass/(float)numBIons);
		System.out.println("B_Int_IntPM: " + numObserbedBIonsIntMassIntPM/(float)numBIons);
		System.out.println("Tryptic: " + numFT/(float)numSpec);
	}
	
	public static void addScanNum() throws Exception
	{
		String specFileName = System.getProperty("user.home")+"/Research/Data/Zubarev/Zubarev_HCD.mgf";
		SpectraIterator itr = new SpectraIterator(specFileName, new MgfSpectrumParser());
		int specNum = -1;
		HashMap<String,Integer> titleSpecNumMap = new HashMap<String,Integer>();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			specNum++;
			titleSpecNumMap.put(spec.getTitle(), specNum);
		}
		
		String fileName = System.getProperty("user.home")+"/Research/Data/Zubarev/Mascot/Zubarev_HCD_FDR01.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String header = in.readLine();
		System.out.println(header+"\tSpecFileName\tSpecNum");
		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length != 7)
				continue;
			String title = token[0];
			specNum = titleSpecNumMap.get(title);
			System.out.println(s+"\tZubarev_HCD.mgf\t"+specNum);
		}
		in.close();
	}
	
//	public static void compositionMSGF() throws Exception
//	{
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		int maxLength = 20;
//		String fileName = System.getProperty("user.home")+"/Research/Data/Zubarev/29novpredictAD_F0_03.mgf";
//		Enzyme enzyme = Enzyme.TRYPSIN;
//
//		long time = System.currentTimeMillis();
//		CompositionFactory allCompositions = new CompositionFactory(aaSet, enzyme, maxLength+2);
//
//		SpectraMap specMap = new SpectraMap(fileName, new MgfSpectrumParser());
//		NewRankScorer scorer = new NewRankScorer(System.getProperty("user.home")+"/Research/Data/CIDETDHCD/AnnotatedSpectra/PNNL_ETD.param");
//	
//		String resFileName = System.getProperty("user.home")+"/Research/Data/Zubarev/Mascot/19Feb_ADfulldigest_HCDETDTop5_2h_f0_01_hcd_target_all.txt";
//		
//		BufferedLineReader in = new BufferedLineReader(resFileName);
//		System.out.println(in.readLine()+"\tSpecProb");
//		String s;
//		int lineNum = 0;
//		int numProcessedSpec = 0;
//		while((s=in.readLine()) != null)
//		{
//			lineNum++;
//			String[] token = s.split("\t");
//			if(token.length < 5)
//				continue;
//			String annotation = token[2];
//			String pepStr = annotation.substring(annotation.indexOf('.')+1,annotation.lastIndexOf('.'));
//			if(pepStr.length() > maxLength)
//				continue;
//			Peptide pep = aaSet.getPeptide(pepStr);
//			if(pep == null)
//				continue;
//			int scanNum = Integer.parseInt(token[0]);
//			int charge = Integer.parseInt(token[1]);
//			Spectrum spec = specMap.getSpectrumByScanNum(scanNum);
//			spec.setCharge(charge);
//			if(Math.abs(spec.getParentMass() - pep.getParentMass()) > 1.5*charge)
//			{
//				System.out.println(s+"\t"+"Parent mass mismatch: " + spec.getParentMass() + " != " + pep.getParentMass());
//				continue;
//			}
//			
//			time = System.currentTimeMillis();
//			CompositionGraph graph = new CompositionGraph(pep.getParentMass(), scorer.getMME(), allCompositions, enzyme);
//			ScoredSpectrum<Composition> scoredSpec = scorer.getScoredSpectrum(spec);
//			scoredSpec.precomputeScores(graph.getIntermediateNodeList());
//			GeneratingFunction<Composition> gf = new GeneratingFunction<Composition>(scoredSpec, graph).enzyme(enzyme).doNotBacktrack().doNotCalcNumber();
//			if(gf.getEnzyme()== null || !enzyme.isCleaved(pep))
//				graph.allowNonEnzymaticCleavage();
//			gf.computeGeneratingFunction();
//			int pepScore = gf.getScore(pep);
//			float specProb = gf.getSpectralProbability(pepScore);
//			specProb *= gf.getNeighboringAAProbability(annotation);
//
//			System.out.println(s+"\t"+specProb);
//			numProcessedSpec++;			
//		}
//	}		
	
	public static void splitHCDETD() throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/Zubarev");
		for(File f : dir.listFiles())
		{
			if(!f.getName().endsWith("mgf") || f.getName().contains("29novpredict"))
				continue;
			System.out.println(f.getName());
			int numHCD = 0;
			int numETD = 0;
			SpectraIterator itr = new SpectraIterator(f.getPath(), new MgfSpectrumParser());
			String fileName = f.getPath();
			String name = fileName.substring(0, fileName.lastIndexOf('.'));
			String hcdFileName = name+"_hcd.mgf";
			String etdFileName = name+"_etd.mgf";
			PrintStream hcdOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(hcdFileName)));
			PrintStream etdOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(etdFileName))); 
			while(itr.hasNext())
			{
				Spectrum spec = itr.next();
				String title = spec.getTitle();
				if(!title.contains("experiment: 1"))
				{
					spec.outputMgf(hcdOut);
					numHCD++;
				}
				else
				{
					spec.outputMgf(etdOut);
					numETD++;
				}
			}
			hcdOut.close();
			etdOut.close();
			System.out.println(numHCD+" HCD and " + numETD + " ETD spectra.");
		}
	}
	
	static class Tag implements Comparable<Tag> {
		int[] index;
		String seq;
		int score;
		
		public Tag(int i1, int i2, char residue)
		{
			index = new int[2];
			index[0] = i1;
			index[1] = i2;
			seq = String.valueOf(residue);
		}
		
		private Tag() {}		

		public Tag score(int score)
		{
			this.score = score;
			return this;
		}
		
		public int getScore()
		{
			return score;
		}
		
		public static Tag join(Tag tag1, Tag tag2)
		{
			Tag newTag = new Tag();
			newTag.index = Arrays.copyOf(tag1.index, tag1.index.length+1);
			newTag.index[newTag.index.length-1] = tag2.getLast();
			newTag.seq = tag1.seq + tag2.seq;
			return newTag;
		}
		
		int getFirst()	{ return index[0]; }
		int getLast()	{ return index[index.length-1];}
		String getSeq()	{ return seq; }
		String getRevSeq()
		{
			StringBuffer rev = new StringBuffer();
			for(int i=seq.length()-1; i>=0; i--)
				rev.append(seq.charAt(i));
			return rev.toString();
		}

		public int compareTo(Tag arg0) {
			return score - arg0.score;
		}
	}
	
	public static void testTagging() throws Exception
	{
		int dataset = 0;	// 0: HCD, 1: CID, 2: ETD
		long time = System.currentTimeMillis();
		String fileName = null;
		if(dataset == 0)
			fileName = System.getProperty("user.home")+"/Research/Data/Zubarev/AnnotatedSpectra/Zubarev_HCD_Annotated.mgf";
		else if(dataset == 1)
			fileName = System.getProperty("user.home")+"/Research/Data/HeckRevision/AnnotatedSpectra/CID_Tryp_Confident.mgf";
		else if(dataset == 2)
			fileName = System.getProperty("user.home")+"/Research/Data/HeckRevision/AnnotatedSpectra/ETD_Tryp_Confident.mgf";
		
		int numMaxTags = 50;
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		int numSpecs = 0;
		int numSpecsWithCorrectTags = 0;
		int sumNumTags = 0;
		
		NewRankScorer scorer = null;
		if(dataset == 0)
			scorer = NewScorerFactory.get(ActivationMethod.HCD, Enzyme.TRYPSIN);
		else if(dataset == 1)
			scorer = NewScorerFactory.get(ActivationMethod.CID, Enzyme.TRYPSIN);
		else if(dataset == 2)
			scorer = NewScorerFactory.get(ActivationMethod.ETD, Enzyme.TRYPSIN);
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		float tagErrTolerance = 0.5f;
		if(dataset == 0)
			tagErrTolerance = 0.02f;
		
		Tolerance fragTolerance = new Tolerance(0.5f);
		if(dataset == 0)
			fragTolerance = new Tolerance(30, true);
		
		float[] nominalMass = new float[200];
		char[] residue = new char[200];
		int[] intAAMass = new int[18];
		int aaMassIndex = 0;
		
		IonType[] ion = null;
		if(dataset == 0 || dataset == 1)
		{
			ion = new IonType[3];
			ion[0] = IonType.Y;
			ion[1] = IonType.getIonType("y2");
			ion[2] = IonType.B;
		}
		else if(dataset == 2)
		{
			ion = new IonType[4];
			ion[0] = IonType.Z;
			ion[1] = IonType.getIonType("z+H");
			ion[2] = IonType.getIonType("z-H");
			ion[1] = IonType.getIonType("y2");
			ion[2] = IonType.B;
		}
		
		for(AminoAcid aa : aaSet)
		{
			int normAA = aa.getNominalMass();
			if(nominalMass[normAA] == 0)
			{
				nominalMass[normAA] = aa.getMass();
				residue[normAA] = aa.getResidue();
				intAAMass[aaMassIndex++] = normAA;
			}
		}
		AminoAcid aaQ = aaSet.getAminoAcid('Q');
		AminoAcid aaK = aaSet.getAminoAcid('K');
		
		NominalMassFactory factory = new NominalMassFactory(aaSet, Enzyme.TRYPSIN, 50);
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			spec = new WindowFilter(10,50).apply(spec);
			if(spec.getCharge() != 3)
				continue;
//			scorer.filterPrecursorPeaks(spec);
			Peptide peptide = spec.getAnnotation();
//			if(Math.abs(peptide.getParentMass()-spec.getParentMass()) <= 0.2f)
//				continue;

//			scorer.filterPrecursorPeaks(spec);
			
			float expPepMass = spec.getParentMass() - (float)Composition.H2O;// - (float)Composition.ISOTOPE;
			int nomExpPepMass = NominalMass.toNominalMass(expPepMass);
			float[] srm = new float[nomExpPepMass+1];
			srm[0] = 0;
			srm[nomExpPepMass] = expPepMass;
			float[] intensity = new float[srm.length];
			
			// compute scorers
			NewScoredSpectrum<NominalMass> scoredSpec = scorer.getScoredSpectrum(spec);
			ArrayList<NominalMass> nodeList = new ArrayList<NominalMass>();
			for(int i=0; i<=nomExpPepMass; i++)
				nodeList.add(new NominalMass(i));
			AminoAcidGraph graph = new AminoAcidGraph(factory, spec.getParentMass(), scoredSpec);
//			AminoAcidGraph graph = new AminoAcidGraph(factory, nomExpPepMass);
//			graph.computeNodeScores(scoredSpec);
			
			intensity[0] = 1000;
			intensity[nomExpPepMass] = 1000;
			numSpecs++;
			
			// y
			for(Peak p : spec)
			{
				float mz = p.getMz();
				int nomSRM = NominalMass.toNominalMass(mz)-19;
				if(nomSRM >0 && nomSRM < srm.length)
				{
					if(p.getIntensity() > intensity[nomSRM])
					{
						srm[nomSRM] = mz - IonType.Y.getOffset();
						intensity[nomSRM] = p.getIntensity();
					}
				}
			}
			
			// y2
			for(Peak p : spec)
			{
				float mz = p.getMz();
				float mass = IonType.getIonType("y2").getMass(mz);
				int nomSRM = NominalMass.toNominalMass(mass);
				if(nomSRM >0 && nomSRM < srm.length)
				{
					if(intensity[nomSRM] == 0 && p.getIntensity() > intensity[nomSRM])
					{
						srm[nomSRM] = mass;
						intensity[nomSRM] = p.getIntensity();
					}
				}
			}

			// b
			for(Peak p : spec)
			{
				float mz = p.getMz();
				float mass = IonType.B.getMass(mz);
				float compMass = expPepMass - mass;
				int nomSRM = NominalMass.toNominalMass(compMass);
				if(nomSRM >0 && nomSRM < srm.length)
				{
					if(intensity[nomSRM] == 0 && p.getIntensity() > intensity[nomSRM])
					{
						srm[nomSRM] = compMass;
						intensity[nomSRM] = p.getIntensity();
					}
				}
			}
			
			ArrayList<Tag> length1Tag = new ArrayList<Tag>();
			// generate length 1 tags
			for(int i=0; i<srm.length-1; i++)
			{
				if(intensity[i] == 0)
					continue;
				float p1 = srm[i];
				for(int aaIndex=0; aaIndex<intAAMass.length; aaIndex++)
				{
					int j = i+intAAMass[aaIndex];
					if(j >= srm.length)
						break;
					if(intensity[j] == 0)
						continue;
					float p2 = srm[j];
					float massDiff = p2-p1;
					int nominalDiff = j-i;
					
					float aaMass = nominalMass[nominalDiff];
					char aaResidue = residue[nominalDiff];
					if(aaMass > 0)
					{
						if(nominalDiff == 128)	// K or Q
						{
							if(Math.abs(massDiff-aaQ.getMass()) > Math.abs(massDiff-aaK.getMass()))
							{
								aaMass = aaK.getMass();
								aaResidue = 'K';
							}
							else
							{
								aaMass = aaQ.getMass();
								aaResidue = 'Q';
							}
						}
						float error = massDiff-aaMass;
						if(error > -tagErrTolerance && error < tagErrTolerance)
						{
							length1Tag.add(new Tag(i,j, aaResidue));
						}
					}
				}
			}

			// generate length 2 tags
			ArrayList<Tag> length2Tag = new ArrayList<Tag>();
			for(int i=0; i<length1Tag.size()-1; i++)
			{
				Tag t1 = length1Tag.get(i);
				int index1 = t1.getLast();
				for(int j=i+1; j<length1Tag.size(); j++)
				{
					Tag t2 = length1Tag.get(j);
					int index2 = t2.getFirst();
					if(index2 == index1)
					{
						length2Tag.add(Tag.join(t1, t2));
					}
					else if(index2 > index1)
						break;
				}
			}

			// generate length 3 tags
			ArrayList<Tag> length3Tag = new ArrayList<Tag>();
			for(int i=0; i<length2Tag.size()-1; i++)
			{
				Tag t1 = length2Tag.get(i);
				int index1 = t1.getLast();
				for(int j=0; j<length1Tag.size(); j++)
				{
					Tag t2 = length1Tag.get(j);
					int index2 = t2.getFirst();
					if(index2 == index1)
					{
						length3Tag.add(Tag.join(t1, t2));
					}
					else if(index2 > index1)
						break;
				}
			}
			
			// compute tag scores
			for(Tag t : length3Tag)
			{
				int score = 0;
				for(int i : t.index)
				{
					if(i==0 || i==nomExpPepMass)
						score += 10;
					else
					{
						NominalMass p = new NominalMass(nomExpPepMass-i);
						NominalMass s = new NominalMass(i);
						score += scoredSpec.getNodeScore(p, s);
					}
				}
				t.score(score);
			}
			Collections.sort(length3Tag, Collections.reverseOrder());
			
			// matching tags
			float[] theoSRM;
			theoSRM = new float[peptide.size()+1];
			theoSRM[0] = 0;
			float m = 0;
			for(int i=0; i<peptide.size(); i++)
			{
				m += peptide.get(peptide.size()-1-i).getMass();
				theoSRM[i+1] = m;
			}
			
			String pepStr = peptide.toString();
			pepStr = pepStr.replaceAll("I", "L");
			int numTags = 0;
			boolean hasCorrectTag = false;
			sumNumTags += Math.min(length3Tag.size(), numMaxTags);
			for(Tag t : length3Tag)
			{
				String tagSeq = t.getRevSeq();
				int matchIndex = pepStr.indexOf(tagSeq);
				if(matchIndex >= 0)
				{
					float error = theoSRM[peptide.size()-(matchIndex+3)] - srm[t.getFirst()];
//						assert(Math.abs(error) < 2);
					if(error > 0 && fragTolerance.isTolerancePPM())
						error = error*1e6f/theoSRM[peptide.size()-(matchIndex+3)];
					if(error < fragTolerance.getValue())
					{
						hasCorrectTag = true;
						numSpecsWithCorrectTags++;
						break;
					}
				}
				if(++numTags >= numMaxTags)
					break;
			}
			if(!hasCorrectTag)
			{
				System.out.println(spec.getScanNum()+"\t"+peptide);
				int rank = 0;
				for(Tag t : length3Tag)
				{
					rank++;
					String tagSeq = t.getRevSeq();
					int matchIndex = pepStr.indexOf(tagSeq);
					boolean isCorrect = false;
					if(matchIndex >= 0)
					{
						float error = theoSRM[peptide.size()-(matchIndex+3)] - srm[t.getFirst()];
						if(error > 0 && fragTolerance.isTolerancePPM())
							error = error*1e6f/theoSRM[peptide.size()-(matchIndex+3)];
						if(error < fragTolerance.getValue())
							isCorrect = true;
					}
					System.out.println(rank+"\t"+t.getRevSeq()+"\t"+srm[t.getFirst()]+"\t"+isCorrect);
					if(isCorrect)
						break;
				}
			}
//			System.out.println(peptide+"\t"+hasCorrectTag);
		}
		System.out.println("NumSpecs\t"+numSpecs);
		System.out.println("NumSpecsWithCorrectTags\t"+numSpecsWithCorrectTags+"\t"+numSpecsWithCorrectTags/(float)numSpecs);
		System.out.println("AvgNumTags\t"+sumNumTags/(float)numSpecs);
		System.out.println("Time\t"+(System.currentTimeMillis()-time)/(float)numSpecs/1000);
	}	
}
