package unexplainedpeak;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import msgf.NominalMass;
import msgf.ScoredSpectrum;
import msgf.Tolerance;
import msscorer.NewRankScorer;
import msscorer.NewScorerFactory;
import msscorer.PrecursorOffsetFrequency;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Enzyme;
import msutil.IonType;
import msutil.Peak;
import msutil.Peptide;
import msutil.Spectrum;

public class GenerateIonDistGraph {
	public static ArrayList<IonType> gen(String filename, int rankLimit, int specCharge, Tolerance tol, int probabiltiy) throws IOException{
		return gen(filename, rankLimit, null, specCharge, tol, probabiltiy);	
	}
	
	public static ArrayList<IonType> gen(String filename, int rankLimit, ArrayList<IonType> sigIons, int specCharge, Tolerance tol,  float probabiltiy) throws IOException{
		//String filename = "/home/kwj/workspace/inputs/HeckRevision/AnnotatedSpectra/CID_Tryp_Confident.mgf";
		//String filename = "/home/kwj/workspace/inputs/shewLengthAll.mgf";
	//	String filename = "/home/kwj/workspace/inputs/HeckRevision/AnnotatedSpectra/CID_Tryp_Confident_new.mgf";
		
		//int rankLimit = 32;
		//float  = 0.03f;// 0.03
	
		Tolerance tolerance = new Tolerance(0.5f, false);
		HashMap<Integer, HashMap<IonType, Float>> iondist;
		
		OLDPeakRanker ranker = new OLDPeakRanker(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys(), tolerance);
		
		System.out.println("getting spectra");
		ArrayList<Spectrum> spectra = Util.getSpectra(filename, specCharge);
		
		System.out.println("getting spectra done");
		
		if(sigIons == null){
			OffsetFrequencyFunction offset = new OffsetFrequencyFunction(spectra, rankLimit, tolerance);
			sigIons = offset.getSignificantIons(probabiltiy);
		}
		/*
		if(filterPrecursorIons){
			for(Spectrum s : spectra){	
				for(IonType off : sigIons){
					if(off instanceof IonType.PrecursorIon){
						s.filterPrecursorPeaks(tolerance, 0, off.getOffset());
					}
				}
			}
		}
		*/
		AnnotatedPeak.setConsideredIons(sigIons);
		OLDPeakRanker.registerPrecursorIons(sigIons);
		
		iondist = getIonDist(spectra, tolerance, ranker, rankLimit, 0, false);
		
		print(iondist, sigIons, rankLimit);
		
		return sigIons;
		/*
		
		System.out.println("*** algorithm 1");
		
		iondist = getIonDist(spectra, tolerance, ranker, rankLimit, 1, false);
		
		print(iondist, sigIons, rankLimit);
		
		System.out.println("*** algorithm 2");
		
		iondist = getIonDist(spectra, tolerance, ranker, rankLimit, 2, false);
		
		print(iondist, sigIons, rankLimit);*/
		
		/*
		// should be at the end since spectra will be changed
		System.out.println("*** exclude precursor ions");
		
		iondist = getIonDist(spectra, tolerance, ranker, rankLimit, 0, true);
		
		print(iondist, sigIons, rankLimit);
		
		System.out.println("*** algorithm 1 after exclusion");
		
		iondist = getIonDist(spectra, tolerance, ranker, rankLimit, 1, false);
		
		print(iondist, sigIons, rankLimit);
		
		System.out.println("*** algorithm 2 after exclusion");
		
		iondist = getIonDist(spectra, tolerance, ranker, rankLimit, 2, false);
		
		print(iondist, sigIons, rankLimit);*/
		
	}
	
	public static HashMap<Integer, HashMap<IonType, Float>> getIonDist(ArrayList<Spectrum> spectra, Tolerance tolerance, OLDPeakRanker ranker, int rankLimit, int algorithmNumber, boolean excludePrecursorIons){
		HashMap<Integer, HashMap<IonType, Float>> iondist = new HashMap<Integer, HashMap<IonType, Float>>();
		for(Spectrum spec: spectra){
			HashMap<Integer, Peak> peaks = null;
			
			if(algorithmNumber == 0 && !excludePrecursorIons)
				peaks = ranker.getPeaksinIntensityOrder(spec, rankLimit);
			else if(algorithmNumber == 0 && excludePrecursorIons)
				peaks = ranker.getPeaksinIntensityOrderExcludingPrecursorIons(spec, rankLimit);
			else if(algorithmNumber == 1)
				peaks = ranker.getPeaksbyAlgorithm1(spec, rankLimit);
			else if(algorithmNumber == 2)
				peaks = ranker.getPeaksbyAlgorithm2(spec, rankLimit);
			
			for(int rank : peaks.keySet()){
				Peak p = peaks.get(rank);
				AnnotatedPeak annotatedPeak = new AnnotatedPeak(p.getMz(), p.getIntensity(), p.getCharge(), spec.getAnnotation(), tolerance);
				if(annotatedPeak.isExplained()){
					HashMap<IonType, Float> v;
					if(iondist.containsKey(rank))
						v = iondist.get(rank);
					else v = new HashMap<IonType, Float>();
					
					for(IonType ion: annotatedPeak.getExplainingIonTypes()){
						if(annotatedPeak.getExplainingIonTypes().size() > 1 && ion instanceof IonType.InternalIon) continue;
						if(annotatedPeak.getExplainingIonTypes().size() > 1 && ion instanceof IonType.CyclicIon) continue;
						
						float num;
						if(v.containsKey(ion)) num = v.get(ion);
						else num = 0;
						
						v.put(ion, ++num);
						//if(ion instanceof IonType.CyclicIon && annotatedPeak.getExplainingIonTypes().size() == 1){
						//	System.out.println(spec.getAnnotationStr() + "\t" + annotatedPeak.getExplainingAminoAcids().get(annotatedPeak.getExplainingIonTypes().indexOf(ion)));
						//}
					}
					
					iondist.put(rank, v);
					
				}
			}
		}
		
		for(int i : iondist.keySet()){
			HashMap<IonType, Float> v = iondist.get(i);
			for(IonType ion : v.keySet()){
				v.put(ion, v.get(ion)/spectra.size());
			}
		}
			
		
		
		return iondist;
	}
	
	static void print(HashMap<Integer, HashMap<IonType, Float>> iondist, ArrayList<IonType> sigIons, int rankLimit){
		System.out.println("dist = [");
		float[] noise = new float[rankLimit];
		for(IonType ion : sigIons){
			System.out.println("%" + Util.getIonString(ion));
			for(int rank = 1; rank <= rankLimit; rank++){
				if(iondist.containsKey(rank) && iondist.get(rank).containsKey(ion)){
						System.out.print("\t" + String.format("%.2f", iondist.get(rank).get(ion)*100));
						noise[rank-1] += iondist.get(rank).get(ion)*100;
				}
				else 
					System.out.print("\t0.0");
			}
			System.out.println();
		}
		System.out.println("%noise");
		for(int rank = 1; rank <= rankLimit; rank++){
			System.out.print("\t" + String.format("%.2f", Math.max(0,100 - noise[rank-1])));
		}
		System.out.println();
		System.out.println("];");
	}
}
