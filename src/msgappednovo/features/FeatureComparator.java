package msgappednovo.features;

import java.util.Comparator;

import msutil.IonType;

public class FeatureComparator implements Comparator<Feature>{ // 
	private IonType ion;
	public FeatureComparator(IonType ion){
		this.ion = ion;
	}
	@Override
	public int compare(Feature arg0, Feature arg1) {
		return Float.compare(Feature.getKLDivergenceFromNullCondition(ion, arg0), Feature.getKLDivergenceFromNullCondition(ion, arg1));
	}

}
