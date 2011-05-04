package unexplainedpeak;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import parser.MgfSpectrumParser;
import unexplainedpeak.OffsetFrequencyFunction.offsetPeak;
import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;

public class Util {
	
	//TODO fix after using mass no mz
	public static String getIonString(IonType io){
		 return Character.toLowerCase(io.toString().charAt(0)) + "_"+io.getCharge()+"("+(io.getOffset()) * io.getCharge()+")";
	}
	
	public static ArrayList<Spectrum> getSpectra(String filename, int specCharge) throws IOException{
		return getSpectra(filename, null, specCharge);
	}
	
	static ArrayList<Spectrum> getSpectra(String filename, ArrayList<Spectrum> spectra, int specCharge) throws IOException{
		//System.out.println("getting spectra");
		
		ArrayList<Spectrum> spectraToWrite = new ArrayList<Spectrum>();
		HashSet<Integer> scanNums = new HashSet<Integer>();
		
		if(spectra != null)
			for(Spectrum spec : spectra)
				scanNums.add(spec.getScanNum());
			
		int s=0;
		Iterator<Spectrum> iterator = new SpectraIterator(filename, new MgfSpectrumParser());
		while(iterator.hasNext()){
			
		//	if(s++>100) break;
			Spectrum spec = iterator.next();
			WindowFilter filter = new WindowFilter(6,50);
			
			if(spec.getCharge() != specCharge) continue;
			if(spec.getAnnotation().isModified()) continue;
			//if(spec.getAnnotationStr().contains("P")) continue;

			if(spectra != null && !scanNums.contains(spec.getScanNum())) continue;
		//	spec = filter.apply(spec);
			spec.setRanksOfPeaks();
			spectraToWrite.add(spec);
		}
		//System.out.println("Done");
		return spectraToWrite;
	}

	static ArrayList<Peak> pickPeaks(Spectrum spec, int rankLimit){
		return pickPeaks(spec, rankLimit, false, null, null);
	}
	
	
	//TODO erase them
	static ArrayList<IonType> sigIons;
	
	public static void sigIons(ArrayList<IonType> si){sigIons = si;}
	
	
	static ArrayList<Peak> pickPeaks(Spectrum spec, int rankLimit, boolean considerAAMassDiff, AminoAcidSet aaSet, Tolerance tolerance){
		ArrayList<Peak> peaks = new ArrayList<Peak>();
		
		if(spec.size() <= rankLimit) return spec;
		
		spec.setRanksOfPeaks();
		
		
		if(considerAAMassDiff){
			/*
			HashMap<Integer, Peak> tmp = new HashMap<Integer, Peak>();
			for(Peak p : spec){
				if(p.getRank() <= rankLimit)
					tmp.put(p.getRank(), p);
			}
			
			float maxProb = 0;
			Peak peak = null;
			for(int i=1; i<=rankLimit; i++){
				Peak p = tmp.get(i);
				HashMap<IonType, Float> dist = IonFinding.getIonDistribution().getDist(spec.getParentMass(), p.getMz());
				if(dist == null || !dist.containsKey(sigIons.get(0))) continue;
				float prob = dist.get(sigIons.get(0)) +  dist.get(sigIons.get(1));
				if(prob > maxProb){
					maxProb = prob;
					peak = p;
				}
				
			}
			
			if(peak!=null)peaks.add(peak);
			//System.out.println(peaks);
			if(peaks.isEmpty()) return pickPeaks(spec, 1);
			else return peaks;
			
			*/
			/*
			float probability = 0.1f;
			float avg=0, var = 0;
			for(Peak p : spec)
				avg += p.getIntensity();
			
			avg /= spec.size();
			
			for(Peak p : spec)
				var += Math.pow((p.getIntensity() - avg), 2);
			
			var /= spec.size()-1;
			
			float threshold = (float)(avg + Math.sqrt(var / probability));
			*/
			HashMap<Integer, Peak> tmp = new HashMap<Integer, Peak>();
			for(Peak p:spec){
				//if(p.getIntensity() >= threshold){
				if(p.getRank() < 100){
				    float ratio = p.getMz()/spec.getAnnotation().getParentMass();
					if(ratio < 0.52) continue;
					
					tmp.put(p.getRank(), p);
				}
			}
			
			int minRank = Integer.MAX_VALUE; 
			Peak p = null;
			Peak pp = null;
			
			for(int r=1; r<=spec.size(); r++){
				Peak p1 = tmp.get(r);

				if(p1 == null) continue;
				
				
				for(int s=r+1; s<=spec.size();s++){
					Peak p2 = tmp.get(s);
					if(p2 == null) continue;
					float mzDiff = Math.abs(p1.getMz() - p2.getMz());
					float toleranceAsDa = tolerance.getToleranceAsDa(mzDiff);
					for(int c=1; c<=spec.getCharge()-1;c++){
						for(AminoAcid aa : aaSet){
							if(Math.abs(mzDiff - aa.getMass()/c) <= toleranceAsDa){
								if(minRank > s){
									minRank = s;
									p = p1;
									pp = p2;
								}
								break;
								//if(!peaks.contains(p1))peaks.add(p1);
							}
						}
			//	float ratio = p.getMz()/spec.getAnnotation().getParentMass();
			//	if(ratio >= 0.6 && ratio < 0.8)
			//		tmp.add(p);
					
					}
				}
			}
			
			if(p != null && pp != null){
				peaks.add(p);
			//	peaks.add(pp);
			}
			if(peaks.isEmpty()) return pickPeaks(spec, rankLimit);
			
			/*for(int r=1; r<=spec.size(); r++){
				Peak p1 = tmp.get(r);keyfloat ratio = p.getMz()/spec.getAnnotation().getParentMass();
				if(ratio >= 0.6 && ratio < 0.7)
					tmp.add(p);
				fofloat ratio = p.getMz()/spec.getAnnotation().getParentMass();
				if(ratio >= 0.6 && ratio < 0.8)
					tmp.add(p);
				
			}
			
			int minrank = Integer.MAX_VALUE;if(peaks.isEmpty()) return pickPeaks(spec, rankLimit);
			Peak minrankpeak = null;
			for(Peak p:tmp){
				if(minrank > p.getRank()){
					minrank = p.getRank();
					minrankpeak = p;
				}
					
			}
			//System.out.println(tmp + "\t" + minrankpeak);
			if(!tmp.isEmpty())peaks.add(minrankpeak);
			else{
				for(Peak p:spec)
					if(p.getRank() <= rankLimit) peaks.add(p);
			}
				r(int s=r+1; s<=spec.size();s++){
					Peak p2 = tmp.get(s);
					if(p1 == null || p2 == null) continue;
					float mzDiff = Math.abs(p1.getMz() - p2.getMz());
					float toleranceAsDa = tolerance.getToleranceAsDa(mzDiff);
					for(int c=1; c<=spec.getCharge()-1;c++){
						for(AminoAcid aa : aaSet){
							if(Math.abs(mzDiff - aa.getMass()/c) <= toleranceAsDa)
								if(!peaks.contains(p1))peaks.add(p1);
						}
					}
				}
				if(peaks.size() >=rankLimit) return peaks;
			}*/
			
			//if(peaks.isEmpty()) return pickPeaks(spec, rankLimit);
			
		}else{
			for(Peak p:spec)
				if(p.getRank() <= rankLimit) peaks.add(p);
		
		}
		
		
		return peaks;
	}
	

}