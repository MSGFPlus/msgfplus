package edu.ucsd.msjava.msscorer;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.InputStream;
import java.util.Hashtable;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Protocol;

public class NewScorerFactory {
	private static final String IONSTAT_RESOURCE_DIR = "ionstat/";
	private NewScorerFactory() {}
	
	public static class SpecDataType {
		public SpecDataType(ActivationMethod method, InstrumentType instType, Enzyme enzyme) {
			this(method, instType, enzyme, Protocol.NOPROTOCOL);
		}

		public SpecDataType(ActivationMethod method, InstrumentType instType, Enzyme enzyme, Protocol protocol) {
			this.method = method;
			this.instType = instType;
			this.enzyme = enzyme;
			this.protocol = protocol;
		}
		
		@Override 
		public boolean equals(Object obj) {
			if(obj instanceof SpecDataType)
			{
				SpecDataType other = (SpecDataType)obj;
				if(this.method == other.method &&
					this.instType == other.instType &&
					this.enzyme == other.enzyme &&
					this.protocol == other.protocol
					)
					return true;
			}
			return false;
		}
		
		@Override
		public int hashCode() {
			return method.hashCode()*(enzyme == null ? 1 : enzyme.hashCode())*instType.hashCode()*(protocol == null ? 1 : protocol.hashCode());
		}
		
		@Override
		public String toString() {
			if(protocol == Protocol.NOPROTOCOL)
				return method.getName()+"_"+instType.getName()+"_"+(enzyme == null ? "null" : enzyme.getName());
			else
				return method.getName()+"_"+instType.getName()+"_"+(enzyme == null ? "null" : enzyme.getName())+"_"+protocol.getName();
		}
		
		public ActivationMethod getActivationMethod()	{ return method; }
		public InstrumentType getInstrumentType()		{ return instType; }
		public Enzyme getEnzyme()						{ return enzyme; }
		public Protocol getProtocol()			{ return protocol; }
		
		private ActivationMethod method;
		private InstrumentType instType;
		private Enzyme enzyme;
		private Protocol protocol;
	}
	
	private static Hashtable<SpecDataType, NewRankScorer> scorerTable = new Hashtable<SpecDataType, NewRankScorer>();
	
	/**
	 * @deprecated Use get(ActivationMethod method, InstrumentType instType, Enzyme enzyme) instead
	 * @param method
	 * @param enzyme
	 * @return
	 */
	public static NewRankScorer get(ActivationMethod method, Enzyme enzyme)
	{
		if(method != ActivationMethod.HCD)
			return get(method, InstrumentType.LOW_RESOLUTION_LTQ, enzyme, Protocol.NOPROTOCOL);
		else
			return get(method, InstrumentType.HIGH_RESOLUTION_LTQ, enzyme, Protocol.NOPROTOCOL);
	}
	
	public static NewRankScorer get(ActivationMethod method, InstrumentType instType, Enzyme enzyme, Protocol protocol)
	{
		if(method == null || method == ActivationMethod.PQD)
			method = ActivationMethod.CID;
		if(enzyme == null)
			enzyme = Enzyme.TRYPSIN;
		if(instType == null)
			instType = InstrumentType.LOW_RESOLUTION_LTQ;
		if(method == ActivationMethod.HCD)
			instType = InstrumentType.HIGH_RESOLUTION_LTQ;
		
		SpecDataType condition = new SpecDataType(method, instType, enzyme, protocol);
		NewRankScorer scorer = scorerTable.get(condition);
		if(scorer != null)
			return scorer;
		
		File userParamFile = new File("params/"+condition+".param");
		if(userParamFile.exists())
		{
			System.out.println("Loading user param file: " + userParamFile.getName());
			scorer = new NewRankScorer(userParamFile.getPath());
			scorerTable.put(condition, scorer);
			return scorer;
		}
		InputStream is = ClassLoader.getSystemResourceAsStream(IONSTAT_RESOURCE_DIR+condition+".param");
		if(is != null)
		{
			scorer = new NewRankScorer(new BufferedInputStream(is));
			scorerTable.put(condition, scorer);
			return scorer;
		}
		return get(method, instType, enzyme);
	}

	private static NewRankScorer get(ActivationMethod method, InstrumentType instType, Enzyme enzyme)
	{
		if(method != null && method == ActivationMethod.FUSION)
			return null;
		
		SpecDataType condition = new SpecDataType(method, instType, enzyme);
		NewRankScorer scorer = scorerTable.get(condition);
		if(scorer == null)
		{
			InputStream is = ClassLoader.getSystemResourceAsStream(IONSTAT_RESOURCE_DIR+condition+".param");
			if(is == null)	// param file does not exist. Change enzyme.
			{
				// change enzyme
				Enzyme alternativeEnzyme;
				if(enzyme.isCTerm())
					alternativeEnzyme = Enzyme.TRYPSIN;
				else
					alternativeEnzyme = Enzyme.LysN;
				SpecDataType newCond = new SpecDataType(method, instType, alternativeEnzyme);
				is = ClassLoader.getSystemResourceAsStream(IONSTAT_RESOURCE_DIR+newCond+".param");
				
				if(is == null)	// if all the above failed, try to use CIDorETD-LowRes-Tryp, CIDorETD-LowRes-LysN, or CID-TOF-Tryp
				{
					if((method == ActivationMethod.HCD)
							&& (instType == InstrumentType.TOF || instType == InstrumentType.HIGH_RESOLUTION_LTQ) 
							&& enzyme.isCTerm())
						newCond = new SpecDataType(ActivationMethod.CID, InstrumentType.TOF, Enzyme.TRYPSIN);
					else if(method.isElectronBased() && enzyme.isCTerm())
						newCond = new SpecDataType(ActivationMethod.ETD, InstrumentType.LOW_RESOLUTION_LTQ, Enzyme.TRYPSIN);
					else if(method.isElectronBased() && enzyme.isNTerm())
						newCond = new SpecDataType(ActivationMethod.ETD, InstrumentType.LOW_RESOLUTION_LTQ, Enzyme.LysN);
					else if(!method.isElectronBased() && enzyme.isNTerm())
						newCond = new SpecDataType(ActivationMethod.CID, InstrumentType.LOW_RESOLUTION_LTQ, Enzyme.LysN);
					else
						newCond = new SpecDataType(ActivationMethod.CID, InstrumentType.LOW_RESOLUTION_LTQ, Enzyme.TRYPSIN);
					is = ClassLoader.getSystemResourceAsStream(IONSTAT_RESOURCE_DIR+newCond+".param");						
				}
			}
			assert(is != null): "param file is missing!: " + method.getName()+" "+enzyme.getName();
			scorer = new NewRankScorer(new BufferedInputStream(is));
			assert(scorer != null): "scorer is null:" + method.getName()+" "+enzyme.getName();
			scorerTable.put(condition, scorer);
		}
		return scorer;
	}	
	
	public static void main(String argv[])
	{
		for(ActivationMethod method : ActivationMethod.getAllRegisteredActivationMethods())
		{
			if(method == ActivationMethod.FUSION)
				continue;
			for(InstrumentType inst : InstrumentType.getAllRegisteredInstrumentTypes())
			{
				for(Enzyme enzyme : Enzyme.getAllRegisteredEnzymes())
				{
					for(Protocol protocol : Protocol.getAllRegisteredProtocols())
					{
						NewRankScorer scorer = NewScorerFactory.get(method, inst, enzyme, protocol);
						System.out.print(method.getName()+"_"+inst.getName()+"_"+enzyme.getName()+"_"+protocol.getName()+" -> ");
						if(scorer != null)
						{
							System.out.println(scorer.getSpecDataType());
						}
						else
						{
							System.err.println("Null!");
							System.exit(-1);
						}
					}
				}
			}
		}
	}
}
