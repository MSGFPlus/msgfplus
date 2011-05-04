package msutil;

import java.util.Hashtable;

public class ActivationMethod {
	private String name;
	private boolean pair;
	private boolean electronBased = false;
	
	private ActivationMethod(String name, boolean isPaired) 
	{
		this.name = name;
		this.pair = isPaired;
	}

	private ActivationMethod(String name) 
	{
		this(name, false);
	}
	
	private ActivationMethod electronBased()
	{
		this.electronBased = true;
		return this;
	}
	
	public String getName()		{ return name; }
	public boolean isPaired()	{ return pair; }
	public boolean isElectronBased() { return electronBased; }
	public static final ActivationMethod CID;
	public static final ActivationMethod ETD;
	public static final ActivationMethod HCD;
	public static final ActivationMethod CID_ETD;
	public static final ActivationMethod CID_HCD;
	
	public static ActivationMethod get(String name)
	{
		return table.get(name);
	}

	public static boolean register(String name)
	{
		return register(name, false);
	}
	
	public static boolean register(String name, boolean isPaired)
	{
		ActivationMethod m = table.get(name);
		if(m != null)
			return false;	// registration was not successful
		else
		{
			ActivationMethod newMethod = new ActivationMethod(name, isPaired);
			table.put(name, newMethod);
			return true;
		}
	}

	@Override
	public String toString() {
		return name;
	}

	@Override
	public boolean equals(Object obj) {
		if(obj instanceof ActivationMethod)
			return this.name.equalsIgnoreCase(((ActivationMethod)obj).name);
		return false;
	}

	@Override
	public int hashCode() {
		return this.name.hashCode();
	}
	

	//// static /////////////
	private static Hashtable<String, ActivationMethod> table;
	static {
		CID = new ActivationMethod("CID");
		ETD = new ActivationMethod("ETD").electronBased();
		HCD = new ActivationMethod("HCD");
		CID_ETD = new ActivationMethod("CID/ETD", true);
		CID_HCD = new ActivationMethod("CID/HCD", true);
		
		table = new Hashtable<String, ActivationMethod>();
		table.put(CID.name, CID);
		table.put(ETD.name, ETD);
		table.put(HCD.name, HCD);
		table.put(CID_ETD.name, CID_ETD);
		table.put(CID_HCD.name, CID_HCD);
	}
}
