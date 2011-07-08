package msutil;

import java.util.HashMap;

public class InstrumentType {
	private String name;
	boolean isHighResolution;
	
	private InstrumentType(String name, boolean isHighResolution) 
	{
		this.name = name;
		this.isHighResolution = isHighResolution;
	}

	public String getName()		{ return name; }
	public boolean isHighResolution()	{ return isHighResolution; }
	
	@Override
	public String toString() {
		return name;
	}

	@Override
	public boolean equals(Object obj) {
		if(obj instanceof InstrumentType)
			return this.name.equalsIgnoreCase(((InstrumentType)obj).name);
		return false;
	}

	@Override
	public int hashCode() {
		return this.name.hashCode();
	}

	public static InstrumentType get(String name)
	{
		return table.get(name);
	}
	
	public static HashMap<String,InstrumentType> table = new HashMap<String,InstrumentType>();
	public static final InstrumentType LOW_RESOLUTION_LTQ;
	public static final InstrumentType TOF;
	public static final InstrumentType HIGH_RESOLUTION_LTQ;
	
	static {
		LOW_RESOLUTION_LTQ = new InstrumentType("LowRes", false);
		HIGH_RESOLUTION_LTQ = new InstrumentType("HighRes", false);
		TOF = new InstrumentType("TOF", true);
		
		table.put(LOW_RESOLUTION_LTQ.getName(), LOW_RESOLUTION_LTQ);
		table.put(HIGH_RESOLUTION_LTQ.getName(), HIGH_RESOLUTION_LTQ);
		table.put(TOF.getName(), TOF);
	}

}
