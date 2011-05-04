package msgappednovo.train;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;


import parser.MgfSpectrumParser;

import msgappednovo.parameters.PeakParameter;
import msgappednovo.parameters.SpectrumParameter;
import msgappednovo.train.InterPeakOffsetFrequencyFunction.NewGeneralizedOffsetPeak;
import msgf.Tolerance;
import msutil.IonType;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;

public class IonSelector {
	private boolean outputIOFF = false;
	private int maxRank;
	private float sigprob = 0.15f;
	//private float sigprobForPrecursorIon = 0.2f;
	private int maxIonNum;
	private WindowFilter filter = null;
	private HashMap<IonType, Float> sigIonIntensityMap = null;
	private String specfilename;
	private Tolerance tol;


	public IonSelector(String specfilename, Tolerance tol, WindowFilter filter, int maxRank, int maxIonNum){
		this.specfilename = specfilename;
		this.tol = tol;
		this.filter = filter;
		this.maxRank = maxRank;
		this.maxIonNum = maxIonNum;
	}
	

	private void sortIons(){
		ArrayList<Float> intensities = new ArrayList<Float>();
		HashMap<IonType, Float> tmpMap = new HashMap<IonType, Float>();
		
		for(IonType ion : sigIonIntensityMap.keySet()){
			intensities.add(sigIonIntensityMap.get(ion));
		}
		
		Collections.sort(intensities, Collections.reverseOrder());
		
		for(float intensity : intensities){
		//	if(intensity < 0.15f) continue;
			for(IonType ion : sigIonIntensityMap.keySet()){
				if(intensity == sigIonIntensityMap.get(ion)){
					tmpMap.put(ion, intensity);
					if(tmpMap.size() > maxIonNum-1) break;
				}
			}
			if(tmpMap.size() > maxIonNum-1) break;
		}
		sigIonIntensityMap = tmpMap;
	}
	
	private void normalizeSigIonIntensityMap(){
		
		float maxRatio = 0;
		for(IonType ion : sigIonIntensityMap.keySet()){
			maxRatio = Math.max(maxRatio, sigIonIntensityMap.get(ion));
		}
		
		for(IonType ion : sigIonIntensityMap.keySet()){
			float ratio = sigIonIntensityMap.get(ion);
			ratio /= maxRatio;
			sigIonIntensityMap.put(ion, ratio);
		}
		sortIons();
	}
	
	private void trainIonPeakIntensityRatio(int specCharge){
		Iterator<Spectrum> iterator;
		HashMap<IonType, Float> numMap = new  HashMap<IonType, Float>();
		
		if(sigIonIntensityMap.isEmpty()) return;
		
		int sn = 0;
		ArrayList<IonType> sigIons = new ArrayList<IonType>();
		
		for(IonType ion : sigIonIntensityMap.keySet()){
			if(! (ion instanceof IonType.PrecursorIon)) sigIons.add(ion);
			sigIonIntensityMap.put(ion, 0f);
		}
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 

			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
		
				sn++;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training ion intensities : " + sn);
				}
				
				spec.setRanksOfPeaks();
				PeakGenerator pgen = new PeakGenerator(spec);
				
				
				float highestintensity = 0;
				
				for(Peak p : spec){
					if(p.getRank() == 1){
						highestintensity = p.getIntensity();
						break;
					}
				}
				
				if(filter != null)
					spec = filter.apply(spec);
				
				
				for(Peak p : spec){
					if(p.getIntensity() == 0) continue;
					if(p.getRank() > maxRank) continue;
					
					boolean isExplained = false;
					for(IonType ion : sigIons){	
						if(!ion.equals(IonType.NOISE)){
							if(pgen.isExplainedBy(p, ion, tol, tol)){
								Float ratio = sigIonIntensityMap.get(ion);
								Float num = numMap.get(ion);
								
								if(ratio == null) ratio = 0f;
								if(num == null) num = 0f;
								
								ratio += p.getIntensity()/highestintensity;
								num ++;
								
								numMap.put(ion, num);
								sigIonIntensityMap.put(ion, ratio);
								isExplained = true;
							}
						}
					}
					/*if(!isExplained){
						Float ratio = sigIonIntensityMap.get(IonType.NOISE);
						Float num = numMap.get(IonType.NOISE);
						
						if(ratio == null) ratio = 0f;
						if(num == null) num = 0f;
						
						ratio += p.getIntensity()/highestintensity;
						num ++;
						
						numMap.put(IonType.NOISE, num);
						sigIonIntensityMap.put(IonType.NOISE, ratio);
					}*/
				}
			}
			
			for(IonType ion : sigIons){
				float ratio = sigIonIntensityMap.get(ion)/numMap.get(ion);
				sigIonIntensityMap.put(ion, ratio);
			}
			
			normalizeSigIonIntensityMap();
			
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 				
	}
	
	
	private int findSigIons(int specCharge){
		sigIonIntensityMap = new  HashMap<IonType, Float>();
		
		Iterator<Spectrum> iterator;
		HashMap<IonType, Float> tmpSigIonMap = new  HashMap<IonType, Float>();
		HashMap<IonType, Float> sigIonMap = new  HashMap<IonType, Float>();
		
		int sn = 0;
		HashMap<InterPeakOffset, Integer> poffsetsMap = new HashMap<InterPeakOffset, Integer>();
		HashMap<InterPeakOffset, Integer> soffsetsMap = new HashMap<InterPeakOffset, Integer>();
		HashMap<InterPeakOffset, Integer> roffsetsMap = new HashMap<InterPeakOffset, Integer>();
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			float normalizer = 0;
			float normalizerForPrecursor = 0;
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				
				Spectrum filteredSpec = spec;
				if(filter != null)filteredSpec = filter.apply(spec);
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
		
				sn++;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("finding significant ions : " + sn);
				}
				
				spec.setRanksOfPeaks();
				PeakGenerator pgen = new PeakGenerator(spec);
				
				
				for(int charge = 1; charge <= spec.getCharge(); charge++){
					
					ArrayList<Peak> comparedPeaks = new ArrayList<Peak>();
					float maxMz = PeakParameter.maxMzWith(charge, spec);
					for(Peak p : filteredSpec){
						if(p.getIntensity() > 0 && p.getMz() < maxMz&& p.getRank() <= maxRank){
							comparedPeaks.add(p);
						}
					}
					
					
					ArrayList<Peak> pbps = pgen.getTheoreticalPrefixBasePeaks(charge);
					ArrayList<Peak> sbps = pgen.getTheoreticalSuffixBasePeaks(charge);
					Peak rbp = pgen.getTheoreticalPrecursorBasePeak(charge);
					
					if(charge == 1){
						normalizer += pbps.size();
						normalizerForPrecursor ++;
					}
					
					ArrayList<InterPeakOffset> poffsets = new ArrayList<InterPeakOffset>();
					ArrayList<InterPeakOffset> soffsets = new ArrayList<InterPeakOffset>();
					ArrayList<InterPeakOffset> roffsets = new ArrayList<InterPeakOffset>();
					
					poffsets.addAll(InterPeakOffsetFrequencyFunction.getGeneralizedOffsets(pbps, comparedPeaks, false, 0, tol));
					soffsets.addAll(InterPeakOffsetFrequencyFunction.getGeneralizedOffsets(sbps, comparedPeaks, false, 0, tol));
					roffsets.addAll(InterPeakOffsetFrequencyFunction.getGeneralizedOffsets(rbp, comparedPeaks, false, 0, tol));
						
					updateOffsetsMap(poffsetsMap, poffsets);
					updateOffsetsMap(soffsetsMap, soffsets);
					updateOffsetsMap(roffsetsMap, roffsets);//TODO something is wrong with precursor ion selection
				}
			}
		
			String iofffilename = null;
			if(outputIOFF) iofffilename = specfilename + "_"  + "Prefix";
			
			for(NewGeneralizedOffsetPeak gofPeak : InterPeakOffsetFrequencyFunction.getOffSetFrequencyFunction(poffsetsMap, normalizer, 0 , iofffilename)){
				InterPeakOffset gof = gofPeak.getGeneralizedOffset();
				tmpSigIonMap.put(new IonType.PrefixIon(gof.getBaseCharge(), gof.getOffset()/gof.getBaseCharge()), gofPeak.getProbability());
			}	
		
			if(outputIOFF) iofffilename = specfilename + "_"  + "Suffix";
			for(NewGeneralizedOffsetPeak gofPeak : InterPeakOffsetFrequencyFunction.getOffSetFrequencyFunction(soffsetsMap, normalizer, 0, iofffilename)){
				InterPeakOffset gof = gofPeak.getGeneralizedOffset();
				tmpSigIonMap.put(new IonType.SuffixIon(gof.getBaseCharge(), gof.getOffset()/gof.getBaseCharge()), gofPeak.getProbability());
			}	
			
			if(outputIOFF) iofffilename = specfilename + "_"  + "Precursor";
			for(NewGeneralizedOffsetPeak gofPeak : InterPeakOffsetFrequencyFunction.getOffSetFrequencyFunction(roffsetsMap, normalizerForPrecursor, 0, iofffilename)){
				InterPeakOffset gof = gofPeak.getGeneralizedOffset();
				tmpSigIonMap.put(new IonType.PrecursorIon(gof.getBaseCharge(), gof.getOffset()/gof.getBaseCharge()), gofPeak.getProbability());
			}	
			
			
			ArrayList<IonType> knownIonTypes = IonType.getAllKnownIonTypes(4, true);
			
			for(IonType ion : tmpSigIonMap.keySet()){
				if(ion instanceof IonType.PrecursorIon) continue;
				
				boolean isKnown = false;
				
				for(IonType kion : knownIonTypes){
					if(ion.getCharge() != kion.getCharge()) continue;
					if(kion instanceof IonType.PrecursorIon) continue;
					
					if((ion.isPrefixIon() && kion.isPrefixIon()) || (!ion.isPrefixIon() && !kion.isPrefixIon())){
						if(Math.abs(ion.getOffset() - kion.getOffset()) < 0.4/ion.getCharge()){
							sigIonMap.put(kion, tmpSigIonMap.get(ion));
							isKnown = true;
							break;
						}
					}
				}
				if(!isKnown) sigIonMap.put(ion, tmpSigIonMap.get(ion));
				
			}
			
			if(!sigIonMap.isEmpty()){
				ArrayList<Float> probs = new ArrayList<Float>();
				for(IonType ion : sigIonMap.keySet()) probs.add(sigIonMap.get(ion));
				
				Collections.sort(probs);
					
				for(IonType ion : sigIonMap.keySet()){
					if(ion instanceof IonType.PrecursorIon){
					//	if(sigIonMap.get(ion) >= sigprobForPrecursorIon)
					//		sigIonIntensityMap.put(ion, 0f);
					}else	if(sigIonMap.get(ion) >= sigprob)
						sigIonIntensityMap.put(ion,sigIonMap.get(ion));
				}
			//	sigIonIntensityMap.put(IonType.NOISE, 0f);
				
				normalizeSigIonIntensityMap();
				System.out.println(sigIonIntensityMap);
			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
		return sn;
	}
	
	static private void updateOffsetsMap(HashMap<InterPeakOffset, Integer> offsetsMap, ArrayList<InterPeakOffset> offsets){
		if(offsets==null) return;

		for(InterPeakOffset gof : offsets){
			Integer n = offsetsMap.get(gof);
			if(n==null) n = 0;
			
			n++;
			
			offsetsMap.put(gof, n);
		}
	}
	
	public int train(int specCharge, boolean considerIonIntensity){
		int sn = findSigIons(specCharge);
		if(considerIonIntensity) trainIonPeakIntensityRatio(specCharge);
		
		int n = 0;
		for(IonType ion : sigIonIntensityMap.keySet()){
			if(ion instanceof IonType.PrecursorIon) continue;
			n++;
		}
		
		System.out.println(sigIonIntensityMap);
		//System.out.println(n);
		
		return sn;
	}
	
	public HashMap<IonType, Float> getSigIons(){
		return sigIonIntensityMap;
	}

	static public void main(String[] args){
		System.out.println(IonType.getAllKnownIonTypes(4, true));
	//	String inputmgf = "/home/kwj/workspace/inputs/Training/CID_Tryp_Confident.mgf";
		//int charge = 3;
		//IonSelector ionfinder = new IonSelector(inputmgf, new Tolerance(0.5f, false), new WindowFilter(6,50), 8);
		//ionfinder.train(charge, true);
	}
	
}
