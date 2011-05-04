package msgappednovo.train;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import parser.MgfSpectrumParser;

import msgappednovo.features.LinkingFeature;
import msgappednovo.features.Feature;
import msgappednovo.features.DensityFeature;
import msgappednovo.features.IonDependencyFeature;
import msgappednovo.features.NullFeature;
import msgappednovo.parameters.PeakParameter;
import msgappednovo.parameters.SpectrumParameter;
import msgf.Tolerance;
import msutil.IonType;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;

class FeatureProbTrainer {
	
	private String specfilename;
	private Tolerance tol;
	private Tolerance pmtol;
	private ArrayList<IonType> sigIons;
	private HashMap<IonType, Float> sigIonMap;
	private HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>> conditionMap;
	private HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>> independentConditionMap;
	private WindowFilter filter = null;
	
	protected FeatureProbTrainer(String specfilename,  HashMap<IonType, Float> sigIonMap, HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>> conditionMap,
			Tolerance tol, Tolerance pmtol){
		this.specfilename = specfilename;
		this.sigIonMap = sigIonMap;
		sigIons = new ArrayList<IonType>();
		
		ArrayList<Float> intensities = new ArrayList<Float>();
		
		for(IonType ion : sigIonMap.keySet()){
			intensities.add(sigIonMap.get(ion));
		}
		
		Collections.sort(intensities);
		
		for(int i = intensities.size()-1; i>=0; i--){
			for(IonType ion : sigIonMap.keySet()){
				if(intensities.get(i) == sigIonMap.get(ion) && !ion.equals(IonType.NOISE)){
					sigIons.add(ion);
				}
			}
		}
		
		this.conditionMap = conditionMap;
		this.tol = tol;
		this.pmtol = pmtol;
	}
	
	protected FeatureProbTrainer filter(WindowFilter b) {filter = b; return this;}
	
	protected void train(int specCharge, int maxRank,  int iterationNum){
		
		int sn = 0;
		try {
			Iterator<Spectrum> iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0;
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				Spectrum filteredspec = spec;
				
				if(spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				if(filter!=null){
					filteredspec = filter.apply(spec);
				}
				
				sn++;
				
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("Iteration: " +iterationNum+ " condition probability training " +": " + sn);
				}
				
				spec.setRanksOfPeaks();
				
				/*
				for(int i=0; i<spec.size(); i++){
					Peak p = spec.get(i);
					if(p.getRank() > maxRank){
						spec.remove(i--);
					}
				}
				*/
				
				PeakGenerator pgen = new PeakGenerator(spec);
				SpectrumParameter spar = new SpectrumParameter(spec);	
				HashMap<PeakParameter, ArrayList<Feature>> sigconMap = conditionMap.get(spar);
				
				for(Peak bp : spec){
					if(bp.getRank() > maxRank) continue;
					
					PeakParameter ppar = new PeakParameter(bp, spec,iterationNum);
					ArrayList<Feature> sigcons = sigconMap.get(ppar);
					ArrayList<Feature> matchedcons = new ArrayList<Feature>();
					NullFeature nullCondition = null;
					IonType expIon = null;
					
					for(IonType ion : sigIons){
						if(pgen.isExplainedBy(bp, ion, tol, pmtol)){
							expIon = ion;
							break;
						}
					}
					
					if(expIon == null) expIon = IonType.NOISE;
					
					for(Feature con : sigcons){
						if(con.holdsFor(bp, filteredspec, tol, pmtol, iterationNum)){
							matchedcons.add(con);
						}
						if(con instanceof NullFeature){
							nullCondition = (NullFeature) con;			
						}	
					}
					
					for(Feature con : matchedcons){
						con.registerNullCondition(nullCondition);
						con.addIonCount(expIon);
					}
				}
			}
		}catch (FileNotFoundException e) {
			System.exit(1);
			e.printStackTrace();
		}
	}
	
	protected void discardConditionsWithLessThan(int minCaseNum, int maxConditionNum){
		discardConditionsWithLessThan(minCaseNum, maxConditionNum, true);
	}
	
	protected void discardConditionsWithLessThan(int minCaseNum, int maxConditionNum, boolean discard){
		ArrayList<Feature> toErase = new ArrayList<Feature>();
		
		int[] num = new int[4];
		int[] erasednum = new int [3];
		ArrayList<Float> KLDs = new ArrayList<Float>();
		
		for(SpectrumParameter spar : conditionMap.keySet()){
			HashMap<PeakParameter, ArrayList<Feature>> cm = conditionMap.get(spar);
			for(PeakParameter ppar : cm.keySet()){
				for(Feature con : cm.get(ppar)){
					if(con instanceof NullFeature){
						num[0]++;
						con.calculateIonProbs();
					}
				}
			}
		}
	
		for(SpectrumParameter spar : conditionMap.keySet()){
			HashMap<PeakParameter, ArrayList<Feature>> cm = conditionMap.get(spar);
			for(PeakParameter ppar : cm.keySet()){
				ArrayList<Feature> cs = cm.get(ppar);
				for(int i=0; i<cs.size();i++){
					Feature con = cs.get(i);
					if(con instanceof NullFeature) continue;
					else if(con instanceof IonDependencyFeature) num[1]++;
					else if(con instanceof LinkingFeature) num[2]++;
					else if(con instanceof DensityFeature) num[3]++;
					
					float sum = con.calculateIonProbs();
					if(discard && sum < minCaseNum){
						cs.remove(i--);
						toErase.add(con);
					}else{
						KLDs.add(con.getKLDivergenceFromNullCondition());
					}
				}
			}
		}
		
		Collections.sort(KLDs, Collections.reverseOrder());
		float KLDthreshold = 0;

		if(!KLDs.isEmpty())
			KLDthreshold = KLDs.get(Math.min(KLDs.size()-1, maxConditionNum));
		
		for(SpectrumParameter spar : conditionMap.keySet()){
			HashMap<PeakParameter, ArrayList<Feature>> cm = conditionMap.get(spar);
			for(PeakParameter ppar : cm.keySet()){
				ArrayList<Feature> cs = cm.get(ppar);
				for(int i=0; i<cs.size();i++){
					Feature con = cs.get(i);
					if(con instanceof NullFeature) continue;
				
					if(discard && con.getKLDivergenceFromNullCondition() < KLDthreshold){
						toErase.add(con);
						cs.remove(i--);
					}
			
				}
			}
		}
		
		
		for(Feature con : toErase){
			if(con instanceof IonDependencyFeature) erasednum[0]++;
			else if(con instanceof LinkingFeature) erasednum[1]++;
			else if(con instanceof DensityFeature) erasednum[2]++;
		}
		
		System.out.println("# Null Condition : " + num[0]);
		System.out.println("# GOF Condition : " + num[1] + " -> " + (num[1] - erasednum[0]));
		System.out.println("# Bridging Condition : " + num[2] + " -> " + (num[2] - erasednum[1]));
		System.out.println("# Density Condition : " + num[3] + " -> " + (num[3] - erasednum[2]));
		
	}
	
	private void updateIndependentConditionMap(){
		/* 0 = null
		 * 1 = bridging
		 * 2 = density
		 * 3 = gof normal
		 * 4 ~ = gof with charge off, comp, diff base peak charges*/
		
		independentConditionMap = new HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>();

		for(SpectrumParameter spar : conditionMap.keySet()){
			int sc = spar.getSpecCharge();
			
			HashMap<PeakParameter, ArrayList<Feature>> cm = conditionMap.get(spar);
			
			if(!independentConditionMap.containsKey(spar)) independentConditionMap.put(spar,  new HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>());
			HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>> si = independentConditionMap.get(spar);
				
			for(PeakParameter ppar : cm.keySet()){
				if(!si.containsKey(ppar)) si.put(ppar, new HashMap<Integer, ArrayList<Feature>>());
				HashMap<Integer, ArrayList<Feature>> ssi = si.get(ppar);
				for(Feature con : cm.get(ppar)){
					int id;
					if(con instanceof NullFeature) id = 0;
					else if(con instanceof LinkingFeature) id = 1;
					else if(con instanceof DensityFeature) id = 2;
					else{
						IonDependencyFeature gcon = (IonDependencyFeature) con;
						if(!gcon.isComplementary() && gcon.getChargeOffset() == 0) id = 3;
						else{
							id = 4 + (gcon.isComplementary()? sc * sc : 0) + sc * (gcon.getBasePeakCharge() - 1) + (gcon.getBasePeakCharge() + gcon.getChargeOffset() - 1);
						}
					}
					if(!ssi.containsKey(id)) ssi.put(id, new ArrayList<Feature>());
					ArrayList<Feature> cons = ssi.get(id);
					
					cons.add(con);
				}
			}	
		}
	}
	
	protected void writeInFile(String file, int charge, int iterationNum){
		updateIndependentConditionMap();
		
		PrintWriter out;
		try {
			out = new PrintWriter(new FileWriter(file, true));
			
			if(iterationNum == 0){
				out.println("#ION\t" + charge);
				for(IonType ion : sigIons){
					if(ion instanceof IonType.SuffixIon){
						out.println("s/"+ion.getCharge()+"/"+ion.getOffset() + "\t" + sigIonMap.get(ion));
					}else if(ion instanceof IonType.PrefixIon){
						out.println("p/"+ion.getCharge()+"/"+ion.getOffset() + "\t" + sigIonMap.get(ion));
					}else if(ion instanceof IonType.PrecursorIon){
						out.println("r/"+ion.getCharge()+"/"+ion.getOffset() + "\t" + sigIonMap.get(ion));
					}
				}
				
				//if(sigIonMap.containsKey(IonType.NOISE));
				//	out.println("n/"+0+"/"+0 + "\t" + sigIonMap.get(IonType.NOISE));
				out.println("#END");
			}
			
			for(SpectrumParameter spar : independentConditionMap.keySet()){
				HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>> si = independentConditionMap.get(spar);
				for(PeakParameter ppar : si.keySet()){
					HashMap<Integer, ArrayList<Feature>> ssi = si.get(ppar);
					for(int i=0; i<100; i++){
						ArrayList<Feature> cons = ssi.get(i);
						if(cons != null && !cons.isEmpty()){
							out.println("#"+i);
							Collections.sort(cons);
							for(int n = cons.size()-1;n>=0;n--){
								Feature con = cons.get(n);
								out.println(con.toFileString());// + "$" + con.getKLDivergenceFromNullCondition() + "$$"+con.getNullCondition());
								for(int j=0; j<sigIons.size(); j++){
									out.print(j+" "+con.getProbability(sigIons.get(j))+"\t");
								}
								out.println();
							}
						}
					}
				}
			}
			out.close();
		}catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/*
	static public void main(String[] args){
		int charge = 2;
		int maxRank = 150;
		int iterationNum = 0;
		
		String inputmgf = "/home/kwj/workspace/inputs/Training/CID_Tryp_Confident.mgf";
		String para = inputmgf + charge +".parNew";
		
		Tolerance tol =  new Tolerance(0.5f, false); 
		Tolerance pmtol =  new Tolerance(2.0f, false);
		
		IonFinder ionfinder = new IonFinder(inputmgf, tol);
		ionfinder.train(charge, 50);
		
		ConditionFinder confinder = new ConditionFinder(inputmgf, ionfinder.getSigIons(), tol);
		confinder.train(charge, maxRank, iterationNum, AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
		
		ConditionProbTrainer contrainer = new ConditionProbTrainer(inputmgf, ionfinder.getSigIons(), confinder.getSignificantConditions(), tol, pmtol);
		contrainer.train(charge, maxRank, iterationNum);
		contrainer.discardConditionsWithLessThan(100, 1e-2f);
		contrainer.writeInFile(para, charge, iterationNum);
	}
	*/
}
