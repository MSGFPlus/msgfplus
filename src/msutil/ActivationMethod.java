package msutil;

import java.util.Hashtable;

public class ActivationMethod {
	private String name;
	private String fullName;
	private boolean electronBased = false;
	
	private ActivationMethod(String name, String fullName) 
	{
		this.name = name;
		this.fullName = fullName;
	}

	private ActivationMethod electronBased()
	{
		this.electronBased = true;
		return this;
	}
	
	public String getName()		{ return name; }
	public String getFullName()	{ return fullName; }
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

	public static boolean register(String name, String fullName)
	{
		ActivationMethod m = table.get(name);
		if(m != null)
			return false;	// registration was not successful
		else
		{
			ActivationMethod newMethod = new ActivationMethod(name, fullName);
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
		CID = new ActivationMethod("CID", "collision-induced dissociation");
		ETD = new ActivationMethod("ETD", "electron transfer dissociation").electronBased();
		HCD = new ActivationMethod("HCD", "high-energy collision-induced dissociation");
		PQD = new ActivationMethod("PQD", "pulsed q dissociation");
		FUSION = new ActivationMethod("FUSION", "merge spectra from the same precursor");
		
		table = new Hashtable<String, ActivationMethod>();
		table.put(CID.name, CID);
		table.put(ETD.name, ETD);
		table.put(HCD.name, HCD);
		table.put(PQD.name, PQD);
		table.put("ETD+SA", ETD);
		table.put(FUSION.name, FUSION);
	}
}
