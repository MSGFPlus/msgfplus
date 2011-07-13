package msscorer;

import java.io.BufferedInputStream;
import java.io.InputStream;
import java.util.Hashtable;

import msutil.ActivationMethod;
import msutil.Enzyme;
import msutil.InstrumentType;

public class NewScorerFactory {
	private NewScorerFactory() {}
	
	private static class Condition {
		public Condition(ActivationMethod method, InstrumentType instType, Enzyme enzyme) {
			this.method = method;
			this.instType = instType;
			this.enzyme = enzyme;
		}
		@Override
		public boolean equals(Object obj) {
			if(obj instanceof Condition)
			{
				Condition other = (Condition)obj;
				if(this.method == other.method &&
					this.instType == other.instType &&
					this.enzyme == other.enzyme)
					return true;
			}
			return false;
		}
		@Override
		public int hashCode() {
			return method.hashCode()*enzyme.hashCode();
		}
		@Override
		public String toString() {
			return method.getName()+"_"+instType.getName()+"_"+enzyme.getName();			
		}
		ActivationMethod method;
		InstrumentType instType;
		Enzyme enzyme;
	}
	
	private static Hashtable<Condition, NewRankScorer> scorerTable = new Hashtable<Condition, NewRankScorer>();
	
	/**
	 * @deprecated Use get(ActivationMethod method, InstrumentType instType, Enzyme enzyme) instead
	 * @param method
	 * @param enzyme
	 * @return
	 */
	public static NewRankScorer get(ActivationMethod method, Enzyme enzyme)
	{
		if(method != ActivationMethod.HCD)
			return get(method, InstrumentType.LOW_RESOLUTION_LTQ, enzyme);
		else
			return get(method, InstrumentType.HIGH_RESOLUTION_LTQ, enzyme);
	}
	
	public static NewRankScorer get(ActivationMethod method, InstrumentType instType, Enzyme enzyme)
	{
		if(method == null)
			method = ActivationMethod.CID;
		if(enzyme == null)
			enzyme = Enzyme.TRYPSIN;
		if(method == ActivationMethod.HCD)
			instType = InstrumentType.HIGH_RESOLUTION_LTQ;
		Condition condition = new Condition(method, instType, enzyme);
		NewRankScorer scorer = scorerTable.get(condition);
		if(scorer == null)
		{
			InputStream is = ClassLoader.getSystemResourceAsStream("resources/ionstat/"+condition+".param");
			if(is == null)	// param file does not exist. Change enzyme.
			{
				// change enzyme
				Enzyme alternativeEnzyme;
				if(enzyme.isCTerm())
					alternativeEnzyme = Enzyme.TRYPSIN;
				else
					alternativeEnzyme = Enzyme.LysN;
				Condition newCond = new Condition(method, instType, alternativeEnzyme);
				is = ClassLoader.getSystemResourceAsStream("resources/ionstat/"+newCond+".param");
				
				if(is == null)	// param file still does not exist. Change method.
				{
					ActivationMethod alternativeMethod;
					if(method.isElectronBased())
						alternativeMethod = ActivationMethod.ETD;
					else
						alternativeMethod = ActivationMethod.CID;
					newCond = new Condition(alternativeMethod, instType, enzyme);
					is = ClassLoader.getSystemResourceAsStream("resources/ionstat/"+newCond+".param");
					
					if(is == null)
					{
						newCond = new Condition(alternativeMethod, instType, alternativeEnzyme);
						is = ClassLoader.getSystemResourceAsStream("resources/ionstat/"+newCond+".param");						
					}
				}
			}
			assert(is != null): "param file is missing!: " + method.getName()+" "+enzyme.getName();
			scorer = new NewRankScorer(new BufferedInputStream(is));
			assert(scorer != null): "scorer is null:" + method.getName()+" "+enzyme.getName();
			scorerTable.put(condition, scorer);
		}
		return scorer;
	}
}
