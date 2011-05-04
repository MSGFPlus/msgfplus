package msgappednovo.features;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Set;

import msgappednovo.parameters.PeakParameter;
import msgappednovo.parameters.SpectrumParameter;
import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.Composition;
import msutil.IonType;
import msutil.Peak;
import msutil.Spectrum;

public abstract class Feature implements Comparable<Feature>{
	private SpectrumParameter spar;
	private PeakParameter ppar;
	private int basePeakCharge;
	private boolean isPresent;
	private int iterationNum;
	private HashMap<IonType, Float> ionProbMap;
	private Feature nullCondition = null;
	
	protected Feature(SpectrumParameter spar, PeakParameter ppar, int basePeakCharge, boolean isPresent, int iterationnum){
		this.spar = spar;
		this.ppar = ppar;
		this.basePeakCharge = basePeakCharge;
		this.isPresent = isPresent;
		this.iterationNum = iterationnum;
		this.ionProbMap = new HashMap<IonType, Float>();
	}

	protected boolean isPresent() {return isPresent;}
	public Set<IonType> getIonTypeSet() {return ionProbMap.keySet();}
	public int getIterationNum() {return iterationNum;}
	public int getBasePeakCharge() {return basePeakCharge;}
	
	public float getKLDivergenceFromNullCondition(){
		if(nullCondition == null) return 0;
		
		assert(!ionProbMap.keySet().contains(IonType.NOISE));
		
		int n = 1;
		for(IonType ion : ionProbMap.keySet()) if(!(ion instanceof IonType.PrecursorIon)) n++;
		
		float kl = 0;
		float[] p1 = new float[n];
		float[] p2 = new float[n];
		int i = 0;
		float[] sum = new float[2];
		
		for(IonType ion:this.ionProbMap.keySet()){
			if(ion instanceof IonType.PrecursorIon) continue;
			p1[i] = nullCondition.ionProbMap.get(ion);
			p2[i] = ionProbMap.get(ion);	
			sum[0] += p1[i];
			sum[1] += p2[i];
			i++;
		}
		
		p1[i] = Math.max(0, 1- sum[0]);
		p2[i] = Math.max(0, 1- sum[1]);
		sum[0] += p1[i];
		sum[1] += p2[i];
		
		for(int j=0; j<p1.length;j++){
			if(p1[j] == 0)continue;
			if(p2[j] == 0) p2[j] = 1e-10f;
			
			p1[j] /= sum[0]; p2[j] /= sum[1]; 
			
			kl+= (float)(p1[j] * Math.log(p1[j]/p2[j]) / Math.log(2) );
		}
		
		return kl;
	}
	
	public SpectrumParameter getSpectrumParameter() {return spar;}
	public PeakParameter getBasePeakParameter() {return ppar;}
	
	public void registerNullCondition(Feature con){nullCondition = con;}
	public Feature getNullCondition() {return nullCondition;}
	public float getProbability(IonType ion){
		assert(!ion.equals(IonType.NOISE));
		assert(!ionProbMap.containsKey(IonType.NOISE));
		if(ionProbMap.containsKey(ion))
			return ionProbMap.get(ion);
		else return 0f;
	}
	
	public int compareTo(Feature o) {
		return new Float(this.getKLDivergenceFromNullCondition()).compareTo(new Float(o.getKLDivergenceFromNullCondition()));
	}
	
	public void addIonCount(IonType ion){
		Float num = ionProbMap.get(ion);
		if(num == null) num = 0f;
		num++;
		ionProbMap.put(ion, num);
	}
	
	//returns sum
	public float calculateIonProbs(){	
		float sum = 0;
		for(IonType ion : ionProbMap.keySet()){
			sum += ionProbMap.get(ion);
		}
		
		for(IonType ion : ionProbMap.keySet()){	
			float p = ionProbMap.get(ion)/sum;
			ionProbMap.put(ion, p);
		}
		
		ionProbMap.remove(IonType.NOISE);
		
		return sum;
		
	}
	
	public void setIonProbMap(HashMap<IonType, Float> ionProbMap){
		this.ionProbMap = ionProbMap;
	}
	
	protected boolean equals(Feature con) {
		if(this == con) return true;
    	
    	return con.isPresent == this.isPresent &&
    		con.iterationNum == this.iterationNum &&
	    	con.basePeakCharge == this.basePeakCharge &&
	    	con.spar.equals(this.spar) &&
	    	con.ppar.equals(this.ppar);
	}
	
	protected boolean holdsFor(Peak bp, Spectrum spec, int iterationNum){
		
		bp.setCharge(basePeakCharge);

		return iterationNum == this.iterationNum &&
		new SpectrumParameter(spec).equals(this.spar) &&
		new PeakParameter(bp, spec, iterationNum).equals(this.ppar);
	}
	
	public abstract boolean equals(Object o);
	
	public abstract int hashCode();
	
	public abstract String toString();
	
	public abstract String toFileString();
	
	public abstract boolean holdsFor(Peak bp, Spectrum spec, Tolerance tol, Tolerance pmtol, int iterationNum);
	
	
	
	static public float maxMzWithIon(IonType ion, Spectrum spec){
		return ion.getMz((float) (spec.getParentMass() - Composition.H2O - AminoAcid.getStandardAminoAcid('G').getMass()));
	}
	
	static public float minMzWithIon(IonType ion){
		return ion.getMz(AminoAcid.getStandardAminoAcid('G').getMass());
	}
	
	// all conditions are assumed to be independent.
	static public HashMap<IonType, Float> getIonProbabilities(Peak bp, ArrayList<Feature> cons, ArrayList<IonType> ions, Spectrum spec){
		assert !ions.contains(IonType.NOISE);
		HashMap<IonType, Float> ret = new HashMap<IonType, Float>();
		
		float[] resultingProbs = new float[ions.size() + 1];
		float[] nullProbs = new float[ions.size() + 1];
		int numcons = 0;
		
		for(int i=0; i<resultingProbs.length; i++) resultingProbs[i] = 1;
		
		for(Feature con : cons){
			float[] probs = new float[ions.size() + 1];
			float sum = 0;
			
			for(int i=0; i<ions.size();i++){
				IonType ion = ions.get(i);
				if(bp.getMz() < minMzWithIon(ion) || bp.getMz() > maxMzWithIon(ion, spec)) probs[i] = 0;
				else 
					probs[i] = con.getProbability(ion);
				sum += probs[i];
			}
			
			probs[probs.length-1] = Math.max(0, 1-sum);
			sum += probs[probs.length-1];
			
			if(con instanceof NullFeature){
				for(int i=0; i<nullProbs.length; i++) nullProbs[i] = probs[i]/sum;
			}else{	
				for(int i=0; i<resultingProbs.length; i++) resultingProbs[i] *= probs[i]/sum;
				numcons++;
			}
		}
		
		for(int i=0; i<resultingProbs.length; i++)
			if(resultingProbs[i] > 0) resultingProbs[i] /= Math.pow(nullProbs[i], numcons-1);
		
		float sum = 0;
		for(float prob : resultingProbs) sum+=prob;

		//assert(sum>0) : (cons);
		if(sum == 0) sum = 1;
	//	System.out.println(sum);
		for(int i=0; i<ions.size();i++){
			IonType ion = ions.get(i);
			ret.put(ion, resultingProbs[i]/sum);
		}
		ret.put(IonType.NOISE, resultingProbs[resultingProbs.length-1]/sum);
		
		return ret;
	}
	
	
	static private float getIonProbability(IonType ion, ArrayList<Feature> cons, Peak bp, Spectrum spec){
		float p = 1;
		float py = 0;
		int n = 0;
		
		for(Feature con : cons){
			if(con instanceof NullFeature) continue;
			
			float mult = 1;
			if(ion.equals(IonType.NOISE)){
				for(IonType i : con.ionProbMap.keySet()){
					if(bp.getMz() < minMzWithIon(i) || bp.getMz() > maxMzWithIon(i, spec)) continue;
					mult -= con.getProbability(i);
				}
				mult = Math.max(mult, 0);
			}else{
				if(bp.getMz() < minMzWithIon(ion) || bp.getMz() > maxMzWithIon(ion, spec)) mult = 0;
				else mult = con.getProbability(ion);
			}
			
			if(py == 0){
				if(ion.equals(IonType.NOISE)){
					py = 1;
					for(IonType i : con.nullCondition.ionProbMap.keySet()){
						if(bp.getMz() < minMzWithIon(i) || bp.getMz() > maxMzWithIon(i, spec)) continue;
						py -= con.nullCondition.getProbability(i);
					}
					py = Math.max(py, 0);
				}else{
					 py = con.nullCondition.getProbability(ion);
				}
			}
				
			p*=mult;
			n ++;
		}
		if(py == 0) return 0;
		return (float) (p/Math.pow(py, n-1));
	}
	
	static public HashMap<IonType, Float> getIonProbabilitiesNew(Peak bp, HashMap<IonType, ArrayList<Feature>> matchedConMap, Spectrum spec){
		HashMap<IonType, Float> ret = new HashMap<IonType, Float>();
		
		for(IonType ion : matchedConMap.keySet()){ // + IOn noise
			ret.put(ion, getIonProbability(ion, matchedConMap.get(ion), bp , spec));
		}
		
		float sum = 0;
		for(IonType i : ret.keySet()){
			sum+=ret.get(i);
		//	System.out.print(": "+ret.get(i)+"\t");
		}
	//	System.out.println();
	//	System.out.println(sum);
		if(sum == 0){
			sum = 1;
			if(ret.containsKey(IonType.NOISE))ret.put(IonType.NOISE, 1f);
		}
	//
		for(IonType i : ret.keySet()){
			ret.put(i, ret.get(i)/sum);
		}

		return ret;
	}
	

	static public float getKLDivergenceFromNullCondition(IonType ion, Feature con){
		if(con.nullCondition == null) return 0;
		float kl = 0;
		float[] p1 = new float[2];
		float[] p2 = new float[2];
		float[] sum = new float[2];
		
		if(ion.equals(IonType.NOISE)){
			float[] t = new float[2];
			for(IonType i : con.ionProbMap.keySet()){
				t[0] += con.nullCondition.ionProbMap.get(i);
				t[1] += con.ionProbMap.get(i);
			}
			p1[0] = Math.max(1-t[0], 0);
			p2[0] = Math.max(1-t[1], 0);
		//	System.out.println(p1[0] + "\t" + p2[0]);
		}else{
			p1[0] = con.nullCondition.ionProbMap.get(ion);
			p2[0] = con.ionProbMap.get(ion);	
		}
	
		sum[0] += p1[0];
		sum[1] += p2[0];
		
		p1[1] = Math.max(0, 1- sum[0]);
		p2[1] = Math.max(0, 1- sum[1]);
		sum[0] += p1[1];
		sum[1] += p2[1];
		
		for(int j=0; j<p1.length;j++){
			if(p1[j] == 0)continue;
			if(p2[j] == 0) p2[j] = Float.MIN_VALUE;
			
			p1[j] /= sum[0]; p2[j] /= sum[1]; 
			
			kl+= (float)(p1[j] * Math.log(p1[j]/p2[j]) / Math.log(2) );
		}
		
		for(int j=0; j<p2.length;j++){
			if(p2[j] == Float.MIN_VALUE)continue;
			if(p1[j] == 0) p1[j] = Float.MIN_VALUE;
			
			//p1[j] /= sum[0]; p2[j] /= sum[1]; 
			
			kl+= (float)(p2[j] * Math.log(p2[j]/p1[j]) / Math.log(2) );
		}
		
		
		return kl;
		
	}
	
}


