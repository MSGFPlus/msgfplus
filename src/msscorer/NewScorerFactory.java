package msscorer;

import java.io.BufferedInputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Hashtable;

import msutil.ActivationMethod;
import msutil.Enzyme;

public class NewScorerFactory {
	private NewScorerFactory() {}
	
	private static class Condition {
		public Condition(ActivationMethod method, Enzyme enzyme) {
			this.method = method;
			this.enzyme = enzyme;
		}
		@Override
		public boolean equals(Object obj) {
			if(obj instanceof Condition)
			{
				Condition other = (Condition)obj;
				if(this.method == other.method &&
					this.enzyme == other.enzyme)
					return true;
			}
			return false;
		}
		@Override
		public int hashCode() {
			return method.hashCode()*enzyme.hashCode();
		}
		ActivationMethod method;
		Enzyme enzyme;
		int charge;
	}
	
	private static Hashtable<Condition, NewRankScorer> scorerTable = new Hashtable<Condition, NewRankScorer>();
	
	public static NewRankScorer get(ActivationMethod method, Enzyme enzyme)
	{
		if(method == null)
			method = ActivationMethod.CID;
		if(enzyme == null)
			enzyme = Enzyme.TRYPSIN;
		Condition condition = new Condition(method, enzyme);
		NewRankScorer scorer = scorerTable.get(condition);
		if(scorer == null)
		{
			InputStream is = ClassLoader.getSystemResourceAsStream("resources/ionstat/"+method.getName()+"_"+enzyme.getName()+".param");
			if(is == null)	// param file does not exist. Change enzyme.
			{
				// change enzyme
				String alternativeEnzyme;
				if(enzyme.isCTerm())
					alternativeEnzyme = Enzyme.TRYPSIN.getName();
				else
					alternativeEnzyme = Enzyme.LysN.getName();
				is = ClassLoader.getSystemResourceAsStream("resources/ionstat/"+method.getName()+"_"+alternativeEnzyme+".param");
				
				if(is == null)	// param file still does not exist. Change method.
				{
					String alternativeMethod;
					if(method.isElectronBased())
						alternativeMethod = ActivationMethod.ETD.getName();
					else
						alternativeMethod = ActivationMethod.CID.getName();
					is = ClassLoader.getSystemResourceAsStream("resources/ionstat/"+alternativeMethod+"_"+enzyme.getName()+".param");
					
					if(is == null)
					{
						is = ClassLoader.getSystemResourceAsStream("resources/ionstat/"+alternativeMethod+"_"+alternativeEnzyme+".param");						
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
