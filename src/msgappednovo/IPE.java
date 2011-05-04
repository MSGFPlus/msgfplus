package msgappednovo;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;


import msgappednovo.features.FeatureComparator;
import msgappednovo.features.LinkingFeature;
import msgappednovo.features.Feature;
import msgappednovo.features.DensityFeature;
import msgappednovo.features.IonDependencyFeature;
import msgappednovo.features.NullFeature;
import msgappednovo.parameters.PeakParameter;
import msgappednovo.parameters.SpectrumParameter;
import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.BufferedLineReader;
import parser.MgfSpectrumParser;

public class IPE {
	static  private int maxIterationNum = 0;
	static private HashMap<Integer, HashMap<IonType, Float>> ionChargeMap = null; // charge key 
	static private HashMap<Integer, ArrayList<IonType>> sortedIonChargeMap = null;
	static  private HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<IonType, HashMap<Integer, ArrayList<Feature>>>>>> independentConditionMap = null;

	private String parafile;
	private Tolerance tol;
	private Tolerance pmtol;
	private AminoAcidSet aaSet;
	private float minKLDivergence = 0;//1e-2f;
	private WindowFilter filter = null;
	
	
	public IPE(String parafile, Tolerance tol, Tolerance pmtol, AminoAcidSet aaSet){
		this.parafile = parafile;
		this.tol = tol;
		this.pmtol = pmtol;
		this.aaSet = aaSet;
	}
	
	public IPE filter(WindowFilter b) {filter = b; return this;}
	
	public AminoAcidSet getAASet() {return aaSet;}
	
	public ArrayList<IonType> getSigIonsOrderedByIntensityWithOutNoiseIon(int charge) {
		if(ionChargeMap == null) readFromFile();
		if(ionChargeMap.get(charge) == null) return new ArrayList<IonType>();
		
		if(sortedIonChargeMap == null) sortedIonChargeMap = new HashMap<Integer, ArrayList<IonType>>();
		if(sortedIonChargeMap.containsKey(charge)) return sortedIonChargeMap.get(charge);
		else sortedIonChargeMap.put(charge, new ArrayList<IonType>());
		
		ArrayList<IonType> ret = sortedIonChargeMap.get(charge);

		ArrayList<Float> intensities = new ArrayList<Float>();
		
		for(IonType ion : ionChargeMap.get(charge).keySet()){
			intensities.add(ionChargeMap.get(charge).get(ion));
		}
		
		Collections.sort(intensities);
		
		for(int i = intensities.size()-1; i>=0; i--){
			for(IonType ion : ionChargeMap.get(charge).keySet()){
				if(intensities.get(i) == ionChargeMap.get(charge).get(ion)){
					ret.add(ion);
				}
			}
		//	if(ret.size() > num) break;
		}
		ret.remove(IonType.NOISE);
		return ret;
	}
	
	
	private void updateSpectrum(Spectrum spec, HashMap<Peak, HashMap<IonType, Float>> profile, boolean isLastIteration){
		 float high = 0, low = 0;
		
		 spec.setRanksOfPeaks();
		 
		ArrayList<IonType> sigIons = new ArrayList<IonType>();
		
		sigIons.addAll(getSigIonsOrderedByIntensityWithOutNoiseIon(spec.getCharge()));
		
		if(!sigIons.contains(IonType.NOISE))sigIons.add(IonType.NOISE);
		HashMap<IonType, Float> ionfactorMap = ionChargeMap.get(spec.getCharge());
			
		 for(Peak p : spec){
			 if(p.getRank() == 1){
				 high = p.getIntensity();
				 break;
			 }
		}
		 
		 for(Peak p : profile.keySet()){
			HashMap<IonType, Float> probs = profile.get(p);
			
			float intensity = 0; 
				
			float factor = 1;
			for(int i=0;i<sigIons.size();i++){
			
				IonType ion = sigIons.get(i);

				if(ion instanceof IonType.PrecursorIon) continue;
				if(ion.equals(IonType.NOISE)) continue;
				//if(isLastIteration && i > 1) continue;
				
				//if(isLastIteration)
					factor = ionfactorMap.get(ion);
				//else{
				//	factor = (float) (1 -  Math.pow((float)i/sigIons.size(),2));
				//}
				//	if(isLastIteration) factor = 10000000;
				intensity += probs.get(ion) *  factor;

			}
			intensity = intensity * (high - low) + low;
			p.setIntensity(intensity);
		 }
	}
	
	private HashMap<Peak, HashMap<IonType, Float>> getProfileForEachIteration(Spectrum spec, int maxRank, int iterationNum){
		HashMap<Peak, HashMap<IonType, Float>> profile = new HashMap<Peak, HashMap<IonType, Float>>();	
		Spectrum filteredspec = spec;
		if(filter != null) filteredspec = filter.apply(spec);
		
		spec.setRanksOfPeaks();
		
		SpectrumParameter spar = new SpectrumParameter(spec);
		HashMap<PeakParameter, HashMap<IonType, HashMap<Integer, ArrayList<Feature>>>> sigconMap = independentConditionMap.get(iterationNum).get(spar);
		//System.out.println(spar + " " + iterationNum);
		if(sigconMap == null) return profile;
		
		PeakParameter.refresh();
		
	//	if(iterationNum > 0 ) maxRank = (int)(maxRank * 0.5);
		for(Peak bp : spec){
			if(bp.getRank() > maxRank) continue;
			PeakParameter ppar = new PeakParameter(bp, spec,iterationNum);
			
			if(ppar.getBasePeakGroupNum() > PeakParameter.getMaxGroupNum()) continue;

			HashMap<IonType, HashMap<Integer, ArrayList<Feature>>> consMap = sigconMap.get(ppar);
			
			HashMap<IonType, ArrayList<Feature>> matchedConMap = new HashMap<IonType, ArrayList<Feature>>();
			
			for(IonType ion : consMap.keySet()){
				HashMap<Integer, ArrayList<Feature>> cons = consMap.get(ion);
				ArrayList<Feature> matchedCons = new ArrayList<Feature>();
				for(int setnum : cons.keySet()){	
					for(Feature con : cons.get(setnum)){
						if(con instanceof DensityFeature) continue;
						if(con.holdsFor(bp, filteredspec, tol, pmtol, iterationNum)){
							matchedCons.add(con);
							break;
						}
					}
				}
				matchedConMap.put(ion, matchedCons);
			}
			
			profile.put(bp, Feature.getIonProbabilitiesNew(bp, matchedConMap, spec));
	
		}
		return profile;
	}
	
	private void updateIndependentConditionMap(Feature con, int setNum){
		int iterationNum = con.getIterationNum();
		SpectrumParameter spar = con.getSpectrumParameter();
		PeakParameter ppar = con.getBasePeakParameter();
		
		if(!independentConditionMap.containsKey(iterationNum))
			independentConditionMap.put(iterationNum, new HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<IonType, HashMap<Integer, ArrayList<Feature>>>>>());
		HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<IonType, HashMap<Integer, ArrayList<Feature>>>>>
			s1 = independentConditionMap.get(iterationNum);
		
		if(!s1.containsKey(spar))
			s1.put(spar, new HashMap<PeakParameter, HashMap<IonType, HashMap<Integer, ArrayList<Feature>>>>());
		HashMap<PeakParameter, HashMap<IonType, HashMap<Integer, ArrayList<Feature>>>>
			s2 = s1.get(spar);
	
		if(!s2.containsKey(ppar))
			s2.put(ppar, new HashMap<IonType, HashMap<Integer, ArrayList<Feature>>>());
		HashMap<IonType, HashMap<Integer, ArrayList<Feature>>> s3 = s2.get(ppar);
	
		ArrayList<IonType> sigIons = new ArrayList<IonType> ();
		for(IonType ion : ionChargeMap.get(spar.getSpecCharge()).keySet()) sigIons.add(ion);

		sigIons.add(IonType.NOISE);
		//System.out.println(sigIons);
		
		for(IonType ion : sigIons){
			if(!s3.containsKey(ion)) s3.put(ion, new HashMap<Integer, ArrayList<Feature>>());
		//	if(NewCondition.getKLDivergenceFromNullCondition(ion, con) < minKLDivergence) continue;
			HashMap<Integer, ArrayList<Feature>> s4 = s3.get(ion);
			
			if(!s4.containsKey(setNum))
				s4.put(setNum, new ArrayList<Feature>());
			ArrayList<Feature> cons = s4.get(setNum);
			//
			int index = Collections.binarySearch(cons, con, Collections.reverseOrder(new FeatureComparator(ion)));
			if(index < 0) index = -index-1;
			if(Feature.getKLDivergenceFromNullCondition(ion, con) > minKLDivergence) cons.add(index, con);
		//	Collections.sort(cons, new ConditionComparator(ion));
		}
	}
	
	public int readFromFile(){
		if(getMaxIterationNum() > 0) return getMaxIterationNum();
		BufferedLineReader in;
		int maxIterationNum = 0;

		try {
			in = new BufferedLineReader(parafile);
			String s;
			int mode = -1;
			int setNum = -1;
			int chargeForIon =  -1;
			HashMap<IonType, Float> ionProbMap = null;
			Feature con = null;
			Feature nullCon = null;
			
			HashMap<Integer, ArrayList<IonType>> sigIonMap = new HashMap<Integer, ArrayList<IonType>>();
			
			while((s = in.readLine()) != null){
				if(s.startsWith("#ION\t")){
					chargeForIon = Integer.parseInt(s.split("\t")[1]);
					if(ionChargeMap == null)
						ionChargeMap = new HashMap<Integer, HashMap<IonType, Float>>();
					
					mode = 0; continue;
				}
				if(s.equals("#END")){
					if(independentConditionMap == null)
						independentConditionMap = new HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<IonType, HashMap<Integer, ArrayList<Feature>>>>>>();
					mode = 1; continue;
				}
				if(s.startsWith("#IONWEIGHT\t")){
					mode = 2; continue;
				}

				if(mode == 0){
					String[] token = s.split("\t");
					if(!ionChargeMap.containsKey(chargeForIon))
						ionChargeMap.put(chargeForIon, new HashMap<IonType, Float>());
					
					HashMap<IonType, Float> ions = ionChargeMap.get(chargeForIon);
					
					if(!token[0].startsWith("n"))
						ions.put(IonType.getIonType(token[0]), Float.parseFloat(token[1]));
					else ions.put(IonType.NOISE, Float.parseFloat(token[1]));
								
					if(!sigIonMap.containsKey(chargeForIon))
						sigIonMap.put(chargeForIon, new ArrayList<IonType>());
					ArrayList<IonType> is = sigIonMap.get(chargeForIon);
					is.add(IonType.getIonType(token[0]));

					
				}else if(mode == 1){
					if(s.startsWith("#")){
						setNum = Integer.parseInt(s.substring(1));
					}else if(s.startsWith("N")){
						con = NullFeature.parseFileString(s);
						nullCon = con;
					}else if(s.startsWith("G")){
						con = IonDependencyFeature.parseFileString(s);
					}else if(s.startsWith("B")){
						con = LinkingFeature.parseFileString(s, aaSet);
					}else if(s.startsWith("D")){
						con = DensityFeature.parseFileString(s);
					}else{
						ionProbMap = new HashMap<IonType, Float>();
						String[] token = s.split("\t");
						
						for(String t : token){
							String[] k = t.split(" ");
							ionProbMap.put(sigIonMap.get(chargeForIon).get(Integer.parseInt(k[0])), Float.parseFloat(k[1]));
						}
						con.setIonProbMap(ionProbMap);
						con.registerNullCondition(nullCon);
						updateIndependentConditionMap(con, setNum);
						maxIterationNum = Math.max(con.getIterationNum() + 1, maxIterationNum);
						
					}
				}
			}
			/*for(int k : ionChargeMap.keySet()){
				System.out.println(ionChargeMap.get(k));
				System.out.println(sigIonMap.get(k));
			}*/
			in.close();
		}catch (IOException e) {
			System.exit(1);
			e.printStackTrace();
		}
		IPE.setMaxIterationNum(maxIterationNum);
		return maxIterationNum;
	}
	
	// Spectrum spec is NOT updated
	public HashMap<Peak, HashMap<IonType, Float>> getProfile(Spectrum spec, int maxBasePeakIntensityRank, int maxIterationNum){ 
		HashMap<Peak, HashMap<IonType, Float>> profile = null;
		 
		maxIterationNum = Math.min(maxIterationNum, readFromFile());

		Spectrum s = spec.getCloneWithoutPeakList();
		for(Peak p : spec){ 
			s.add(p.clone());//TODO erase in near future...??
		}
		
		for(int iterationNum = 0; iterationNum < maxIterationNum; iterationNum++){
			 profile = getProfileForEachIteration(s, maxBasePeakIntensityRank, iterationNum);
			 updateSpectrum(s, profile, iterationNum == maxIterationNum-1);
		}
		
		return profile;	
	}
	
	public void run(String specfilename, String outfilename, int specCharge, int maxRank,  int maxIterationNum, int iterationStartNum, boolean train){
		
		long time = System.currentTimeMillis();
		
		IPE.setMaxIterationNum(0);
		maxIterationNum = Math.min(maxIterationNum, readFromFile());
		
		int sn = 0;
		try {
			PrintStream out = null;
			if(outfilename != null) out = new PrintStream(outfilename);

			Iterator<Spectrum> iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0;
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;

				
				sn++;
				//if(sn >  1000)break;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("running " +": " + sn);
				}
				
	
				for(int iterationNum = iterationStartNum; iterationNum < maxIterationNum; iterationNum++){
					 HashMap<Peak, HashMap<IonType, Float>> profile = getProfileForEachIteration(spec, maxRank, iterationNum);
					 updateSpectrum(spec, profile, (!train) && iterationNum == maxIterationNum-1);
				}
				
				if(outfilename != null){
					spec.outputMgf(out);
				}
			}
			out.close();
		}catch (FileNotFoundException e) {
			System.exit(1);
			e.printStackTrace();
		}
		
		System.out.println("Running Done : " + (float)(System.currentTimeMillis() - time)/sn/1000 + " sec/spec");
		
	}

	public static void setMaxIterationNum(int maxIterationNum) {
		IPE.maxIterationNum = maxIterationNum;
	}

	public static int getMaxIterationNum() {
		return maxIterationNum;
	}
}
