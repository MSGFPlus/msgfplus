package unexplainedpeak;

import java.util.ArrayList;
import java.util.HashSet;

import msgf.Tolerance;
import msutil.Composition;
import msutil.IonType;
import msutil.Peak;
import msutil.Peptide;

public class AnnotatedPeak extends Peak{
	static ArrayList<IonType> consideredIons = new ArrayList<IonType> ();
	
	private boolean isExplained = false;
	private ArrayList<IonType> explainingIonTypes;
	private ArrayList<String> explainingAminoAcids;
	private Peptide annotation;
	
	public static void setConsideredIons(ArrayList<IonType> ions){
		consideredIons = ions;
	}
	

	static ArrayList<IonType> getConsideredIons(){
		return consideredIons;
	}
	
	
	public ArrayList<IonType> getExplainingIonTypes(){return explainingIonTypes;}
	public ArrayList<String> getExplainingAminoAcids() {return explainingAminoAcids;}
	
	public AnnotatedPeak(float mz, float intensity, int charge, Peptide annotation, Tolerance tolerance){
		super(mz, intensity, charge);
		setAnnotation(annotation);
		annotate(tolerance);
	}
	
	public void setAnnotation(Peptide annotation){
		this.annotation = annotation;
	}
	
	public boolean isExplained(){
		return isExplained;
	}
	
	public void annotate(Tolerance tolerance){
		explainingIonTypes = new ArrayList<IonType>();
		explainingAminoAcids = new ArrayList<String>();
		
		// prefix masses
		float[] prefixMasses = annotation.getPRMMasses(true, 0);
		for(int i=0; i<prefixMasses.length; i++){
			float prefixMass = prefixMasses[i];		
			float toleranceInDa = tolerance.getToleranceAsDa(prefixMass);
			
			for(IonType ion : consideredIons){
				if(ion instanceof IonType.PrefixIon){
					if(Math.abs(this.getMz() - ion.getMz(prefixMass)) * ion.getCharge() < toleranceInDa){
						explainingIonTypes.add(ion);
						explainingAminoAcids.add(annotation.toString().substring(0, i+1));
					}
					
				}
			}
		}
		
		// suffix masses
		float[] suffixMasses = annotation.getPRMMasses(false, 0);
		for(int i=0; i<suffixMasses.length; i++){
			float suffixMass = suffixMasses[i];		
			float toleranceInDa = tolerance.getToleranceAsDa(suffixMass);
			
			for(IonType ion : consideredIons){
				if(ion instanceof IonType.SuffixIon){
					if(Math.abs(this.getMz() - ion.getMz(suffixMass)) * ion.getCharge()  <= toleranceInDa){
						explainingIonTypes.add(ion);
						explainingAminoAcids.add(annotation.toString().substring(suffixMasses.length - i, suffixMasses.length+1));
					}
					
				}
			}
		}
		
		// internal masses
		for(int i=1; i<annotation.size()-1; i++){
			for(int j=i+1; j<annotation.size()-1;j++){
				float internalMass = annotation.getMass(i, j);
				float toleranceInDa = tolerance.getToleranceAsDa(internalMass);
				for(IonType ion : consideredIons){
					if(ion instanceof IonType.InternalIon){
						if(Math.abs(this.getMz() - ion.getMz(internalMass)) * ion.getCharge()  <= toleranceInDa){
							explainingIonTypes.add(ion);
							explainingAminoAcids.add(annotation.toString().substring(i, j));
						}
					}
				}
			}
		}
		
		// cyclic masses
		for(int i=2; i<annotation.size();i++){
			for(int j=annotation.size()+1;;j++){
				if(j-i >= annotation.size()) break;
				float cyclicMass = annotation.getMass(i, annotation.size());
				cyclicMass += annotation.getMass(0, j-annotation.size());
				float toleranceInDa = tolerance.getToleranceAsDa(cyclicMass);
				for(IonType ion : consideredIons){
					if(ion instanceof IonType.CyclicIon){
						if(Math.abs(this.getMz() - ion.getMz(cyclicMass))  * ion.getCharge() <= toleranceInDa){
							explainingIonTypes.add(ion);
							explainingAminoAcids.add(annotation.toString().substring(i, annotation.size()) + "->" + annotation.toString().substring(0, j-annotation.size()));
						}
					}
				}
			}
		}
		
		// precursor ion
		for(IonType ion:consideredIons){
			if(ion instanceof IonType.PrecursorIon){
				float precrusorMass = annotation.getParentMass();
				float toleranceInDa = tolerance.getToleranceAsDa(precrusorMass);
				if(Math.abs(this.getMz() - ion.getMz(precrusorMass)) * ion.getCharge()  <= toleranceInDa){
					explainingIonTypes.add(ion);
				}
			}
		}
		
		isExplained = !explainingIonTypes.isEmpty();
	}
}
