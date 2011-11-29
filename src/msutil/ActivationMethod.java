package msutil;

import java.util.Hashtable;

public class ActivationMethod {
	private String name;
	private boolean electronBased = false;
	
	private ActivationMethod(String name) 
	{
		this.name = name;
	}

	private ActivationMethod electronBased()
	{
		this.electronBased = true;
		return this;
	}
	
	public String getName()		{ return name; }
	public boolean isElectronBased() { return electronBased; }
	
	public static final ActivationMethod CID;
	public static final ActivationMethod ETD;
	public static final ActivationMethod HCD;
	public static final ActivationMethod PQD;
	public static final ActivationMethod FUSION;
	
	public static ActivationMethod get(String name)
	{
		return table.get(name);
	}

	public static boolean register(String name)
	{
		ActivationMethod m = table.get(name);
		if(m != null)
			return false;	// registration was not successful
		else
		{
			ActivationMethod newMethod = new ActivationMethod(name);
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
	public static ActivationMethod[] getAllRegisteredActivationMethods()
	{
		return table.values().toArray(new ActivationMethod[0]);
	}
	
	private static Hashtable<String, ActivationMethod> table;
	
	static {
		CID = new ActivationMethod("CID");
		ETD = new ActivationMethod("ETD").electronBased();
		HCD = new ActivationMethod("HCD");
		PQD = new ActivationMethod("PQD");
		FUSION = new ActivationMethod("FUSION");
		
		table = new Hashtable<String, ActivationMethod>();
		table.put(CID.name, CID);
		table.put(ETD.name, ETD);
		table.put(HCD.name, HCD);
		table.put(PQD.name, PQD);
		table.put("ETD+SA", ETD);
		table.put(FUSION.name, FUSION);
	}
}
