package msgappednovo.train;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import parser.MgfSpectrumParser;

import msgappednovo.features.LinkingFeature;
import msgappednovo.features.Feature;
import msgappednovo.features.DensityFeature;
import msgappednovo.features.IonDependencyFeature;
import msgappednovo.features.NullFeature;
import msgappednovo.parameters.PeakParameter;
import msgappednovo.parameters.SpectrumParameter;
import msgappednovo.train.InterPeakOffsetFrequencyFunction.NewGeneralizedOffsetPeak;
import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;

class FeatureSelector {
	private boolean outputIOFF = false;
	
	private float sigprob = 0.10f;
	
	private HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>> conditionMap = null;
	private String specfilename;
	private Tolerance tol;
	private Tolerance pmtol;
	private ArrayList<IonType> sigIons;
	private WindowFilter filter = null;
	
	protected FeatureSelector(String specfilename, ArrayList<IonType> sigIons, Tolerance tol, Tolerance pmtol){
		this.specfilename = specfilename;
		this.sigIons = sigIons;
		this.tol = tol;
		this.pmtol = pmtol;
	}
	
	protected FeatureSelector filter(WindowFilter b) {filter = b; return this;}
	
	private void updateGOFMap(HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<InterPeakOffset, Integer>>>> 
	 gofMap, SpectrumParameter spar, PeakParameter ppar, int peakIntensityRatio, InterPeakOffset gof){
		
		if(gof == null) return;
		
		if(!gofMap.containsKey(spar)) gofMap.put(spar, new HashMap<PeakParameter, HashMap<Integer, HashMap<InterPeakOffset, Integer>>>());
		HashMap<PeakParameter, HashMap<Integer, HashMap<InterPeakOffset, Integer>>> subGofMap1 = gofMap.get(spar);
	
		if(!subGofMap1.containsKey(ppar)) subGofMap1.put(ppar,  new HashMap<Integer, HashMap<InterPeakOffset, Integer>>());
		HashMap<Integer, HashMap<InterPeakOffset, Integer>> subGofMap2 = subGofMap1.get(ppar);
	
		if(!subGofMap2.containsKey(peakIntensityRatio)) subGofMap2.put(peakIntensityRatio, new HashMap<InterPeakOffset, Integer>());
		HashMap<InterPeakOffset, Integer> subGofMap3 = subGofMap2.get(peakIntensityRatio);
		
		Integer num = subGofMap3.get(gof);
		if(num == null) num = 0;
		num ++;
		
		subGofMap3.put(gof, num);
	}
	
	private void updateSigGOFCons(
		HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<IonDependencyFeature>>>
		sigGofCons,
		HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<InterPeakOffset, Integer>>>> 
		gofMap, 
		HashMap<SpectrumParameter, HashMap<PeakParameter, Integer>> 
		numMap,
		int iterationNum,
		IonType ion)
	{
		for(SpectrumParameter spar : gofMap.keySet()){
			HashMap<PeakParameter, HashMap<Integer, HashMap<InterPeakOffset, Integer>>> subGofMap1 = gofMap.get(spar);
			if(!sigGofCons.containsKey(spar)) sigGofCons.put(spar, new HashMap<PeakParameter, ArrayList<IonDependencyFeature>>());
			HashMap<PeakParameter, ArrayList<IonDependencyFeature>> subSigGofCons = sigGofCons.get(spar);
		
			for(PeakParameter ppar : subGofMap1.keySet()){
				HashMap<Integer, HashMap<InterPeakOffset, Integer>> subGofMap2 = subGofMap1.get(ppar);
				if(!subSigGofCons.containsKey(ppar)) subSigGofCons.put(ppar, new ArrayList<IonDependencyFeature>());
				ArrayList<IonDependencyFeature> gofcs = subSigGofCons.get(ppar);
			
				for(int peakIntensityRatio : subGofMap2.keySet()){
					String iofffilename = null;
					if(outputIOFF) {
						new File(specfilename+"_IOFFS").mkdir();
						iofffilename = specfilename+"_IOFFS" + File.separatorChar +  ion  + " "  + spar  + " " + ppar + " ratio: " + peakIntensityRatio + " iteration: " + iterationNum;
					}
					for(NewGeneralizedOffsetPeak gofPeak : InterPeakOffsetFrequencyFunction.getOffSetFrequencyFunction(subGofMap2.get(peakIntensityRatio), numMap.get(spar).get(ppar), sigprob, iofffilename)){
						InterPeakOffset gof = gofPeak.getGeneralizedOffset();
						IonDependencyFeature gofConp = new IonDependencyFeature(spar, ppar, gof.getBaseCharge(), true, iterationNum, peakIntensityRatio, gof);
						IonDependencyFeature gofCona = new IonDependencyFeature(spar, ppar, gof.getBaseCharge(), false, iterationNum, peakIntensityRatio, gof);
						
						if(!gofcs.contains(gofConp)) gofcs.add(gofConp);
						if(!gofcs.contains(gofCona)) gofcs.add(gofCona);
					}
				}
			}
		}
	}

	
	private void  
		findSigGOFCons(HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<IonDependencyFeature>>> 
		sigGofCons,
		IonType ion,
		int specCharge, int maxRank, int iterationNum){
		
		Iterator<Spectrum> iterator;
		
		//HashMap<NewSpectrumParameter, HashMap<NewBasePeakParameter, ArrayList<NewGOFCondition>>>
		//	sigGofCons = new HashMap<NewSpectrumParameter, HashMap<NewBasePeakParameter, ArrayList<NewGOFCondition>>>();
		HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<InterPeakOffset, Integer>>>> 
			gofMap = new HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<InterPeakOffset, Integer>>>>();
		HashMap<SpectrumParameter, HashMap<PeakParameter, Integer>> 
			numMap = new HashMap<SpectrumParameter, HashMap<PeakParameter, Integer>>();
	
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; int sn = 0;
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				Spectrum filteredspec = spec;
				
				if(spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				if(filter != null){
					filteredspec = filter.apply(spec);
				}
				
				sn++;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("Iteration: " +iterationNum+ " finding significant GOFs for " + ion + ": " + sn);
				}
				
				spec.setRanksOfPeaks();
				
				/*for(int i=0; i<spec.size(); i++){
					Peak p = spec.get(i);
					if(p.getRank() > maxRank){
						spec.remove(i--);
					}
				}
				*/
				
				SpectrumParameter spar = new SpectrumParameter(spec);
				PeakGenerator pgen = new PeakGenerator(spec);
				
				if(!numMap.containsKey(spar)) numMap.put(spar,  new HashMap<PeakParameter, Integer>());
				HashMap<PeakParameter, Integer> nums = numMap.get(spar);
				
				for(Peak bp : spec){
					if(bp.getRank() > maxRank) continue;
						
					if(ion.equals(IonType.NOISE)){
						boolean expd = false;
						for(IonType i : sigIons) if(pgen.isExplainedBy(bp, i, tol, pmtol)) {expd = true; break;}
						if(expd) continue;
					}else	if(!pgen.isExplainedBy(bp, ion, tol, pmtol)) continue;
					
					PeakParameter ppar = new PeakParameter(bp, spec, iterationNum);
					Integer n = nums.get(ppar);
					if(n == null) n = 0;
					nums.put(ppar, n+1);
					
					int charge = ion.getCharge();
					
					for(int chargeOffset=1 - charge; chargeOffset<=spec.getCharge() - charge; chargeOffset++){
						for(int i=0;i<2;i++){
							Peak nbp = PeakGenerator.getChargeChangedBasePeak(bp, charge, chargeOffset);
							if(i==1) nbp = PeakGenerator.getComplementaryBasePeak(nbp, nbp.getCharge(), spec);
							
							float minMass = nbp.getMz() + InterPeakOffsetFrequencyFunction.MIN/nbp.getCharge();
							float maxMass = nbp.getMz() + InterPeakOffsetFrequencyFunction.MAX/nbp.getCharge();
							
							ArrayList<Peak> compPeaks = filteredspec.getPeakListByMassRange(minMass, maxMass);
							if(nbp.getCharge() == filteredspec.getCharge())
								compPeaks.add(filteredspec.getPrecursorPeak());
							
							for(Peak cp : compPeaks){ 
								if(cp.equals(bp)) continue;			 
							//	if(cp.getRank() > maxRank) continue;
								int peakIntensityRatio = PeakParameter.getPeakIntensityRatioNum(bp, cp, spec);

								InterPeakOffset off = InterPeakOffsetFrequencyFunction.getGeneralizedOffset(nbp, cp, i==1, chargeOffset, tol);
								updateGOFMap(gofMap, spar, ppar, peakIntensityRatio, off);
							}
						}
					}
				}	
			}
			
			updateSigGOFCons(sigGofCons, gofMap, numMap, iterationNum, ion);
			
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
	}
	
	protected void train(int specCharge, int maxRank, int iterationNum, AminoAcidSet aaSet){
		conditionMap = new HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>>();
		HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<IonDependencyFeature>>>
		sigGofCons = new HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<IonDependencyFeature>>>();
			
		for(IonType ion : sigIons){
			if(ion instanceof IonType.PrecursorIon) continue;
			findSigGOFCons(sigGofCons, ion, specCharge, maxRank, iterationNum);
		}
		
		//findSigGOFCons(sigGofCons, IonType.NOISE, specCharge, 30, iterationNum);
		
		HashSet<Integer> charges = new HashSet<Integer>();
		
		for(IonType ion : sigIons){
			charges.add(ion.getCharge());
		}
		
		for(SpectrumParameter spar : sigGofCons.keySet()){
			if(!conditionMap.containsKey(spar)) conditionMap.put(spar, new HashMap<PeakParameter, ArrayList<Feature>>());
			HashMap<PeakParameter, ArrayList<Feature>> subConditionMap = conditionMap.get(spar);
		
			for(PeakParameter ppar : sigGofCons.get(spar).keySet()){
				if(!subConditionMap.containsKey(ppar)) subConditionMap.put(ppar, new ArrayList<Feature>());
				ArrayList<Feature> cons = subConditionMap.get(ppar);
				
				cons.addAll(sigGofCons.get(spar).get(ppar));
			}
		}
		
		float[] n = new float[2];
		for(SpectrumParameter spar : SpectrumParameter.getAllSpectrumParameters(specCharge)){
			if(!conditionMap.containsKey(spar)) conditionMap.put(spar, new HashMap<PeakParameter, ArrayList<Feature>>());
			HashMap<PeakParameter, ArrayList<Feature>> subConditionMap = conditionMap.get(spar);
			
			for(PeakParameter ppar : PeakParameter.getAllBasePeakParameters(spar.getSpecCharge())){
				if(!subConditionMap.containsKey(ppar)) subConditionMap.put(ppar, new ArrayList<Feature>());
				ArrayList<Feature> cons = subConditionMap.get(ppar);
				
				/*
				for(int ratio : NewPeakParameter.getAllPeakIntensityRatioNums()){
					Peptide pep = new Peptide("PEPTIDE");
					HashSet<InterPeakOffset> gofs = new HashSet<InterPeakOffset>();
					
					for(IonType ion1 : sigIons){
						if(ion1.equals(IonType.NOISE)) continue;
						boolean ip1 = ion1.isPrefixIon();
						
						Peak bp = new Peak(ion1.getMz(pep.getMass(0, 2)), 1, ion1.getCharge()); 
						for(IonType ion2 : sigIons){
							if(ion2.equals(IonType.NOISE) || ion2.equals(ion1)) continue;
							boolean ip2 = ion2.isPrefixIon();
							Peak nbp = NewPeakGenerator.getChargeChangedBasePeak(bp, ion1.getCharge(), ion2.getCharge()-ion1.getCharge());
							if(ip1 != ip2) nbp = NewPeakGenerator.getComplementaryBasePeak(nbp, nbp.getCharge(), pep.getParentMass());						
						
							Peak cp;
							if(ip1 != ip2)
								cp = new Peak(ion2.getMz(pep.getMass(2, pep.size())), 1, ion2.getCharge()); 
							else
								cp = new Peak(ion2.getMz(pep.getMass(0, 2)), 1, ion2.getCharge()); 
							
							InterPeakOffset gof  = InterPeakOffsetFrequencyFunction.getGeneralizedOffset(nbp, cp, ip1 != ip2, ion2.getCharge()-ion1.getCharge(), tol);
							if(gof == null){
						//		System.out.println(ion1 + " " + ion2);
								continue;
							}
							
							if(gofs.contains(gof)) continue;
							
							gofs.add(gof);
						//	System.out.println(ion1 + " " + ion2 + " " + gof) ;
						//	System.out.println(gof);
							IonDependencyFeature gofConp = new IonDependencyFeature(spar, ppar, gof.getBaseCharge(), true, iterationNum, ratio, gof);
							IonDependencyFeature gofCona = new IonDependencyFeature(spar, ppar, gof.getBaseCharge(), false, iterationNum, ratio, gof);
							cons.add(gofConp);
							cons.add(gofCona);
						}	
					}
				}
				*/
				for(int charge : charges){
					for(int ratio : PeakParameter.getAllPeakIntensityRatioNums()){
						cons.add(new LinkingFeature(spar, ppar, charge, true, iterationNum, ratio, aaSet));
						cons.add(new LinkingFeature(spar, ppar, charge, false, iterationNum, ratio, aaSet));
						//cons.add(new DensityFeature(spar, ppar, charge, true, iterationNum, ratio));
					//	cons.add(new DensityFeature(spar, ppar, charge, false, iterationNum, ratio));
					}
					cons.add(new NullFeature(spar, ppar, charge, iterationNum));
				}
				n[0] ++; n[1] += cons.size();
			}
		}
		System.out.println("Average # cons per one peak: " + n[1] / n[0]);
	}
	
	protected HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>> getSignificantConditions(){
		return conditionMap;
	}
	
	/*
	static public void main(String[] args){
		String inputmgf = "/home/kwj/workspace/inputs/Training/CID_Tryp_Confident.mgf";
		int charge = 2;
		
		IonFinder ionfinder = new IonFinder(inputmgf, new Tolerance(0.5f, false));
		ionfinder.train(charge, 50);
		System.out.println(ionfinder.getSigIons());
		
		ConditionFinder confinder = new ConditionFinder(inputmgf, ionfinder.getSigIons(), new Tolerance(0.5f, false));
		confinder.train(charge, 150, 0, AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
	}
	*/
}
