package msutil;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import params.ParamObject;

public class ActivationMethod implements ParamObject {
	private final String name;
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
	@Override
	public String getDescription()	{ return name; }
	public boolean isElectronBased() { return electronBased; }
	
	public static final ActivationMethod ASWRITTEN;
	public static final ActivationMethod CID;
	public static final ActivationMethod ETD;
	public static final ActivationMethod HCD;
	public static final ActivationMethod PQD;
	public static final ActivationMethod FUSION;
	
	public static ActivationMethod get(String name)
	{
		return table.get(name);
	}

	public static ActivationMethod addOrChange(String name, String fullName)
	{
		ActivationMethod m = table.get(name);
		if(m != null)	// change
		{
			m.fullName = fullName;
			return m;
		}
		else
		{
			ActivationMethod newMethod = new ActivationMethod(name, fullName);
			table.put(name, newMethod);
			return newMethod;
		}
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
		return table.keySet().toArray(new ActivationMethod[0]);
	}
	
	private static LinkedHashMap<String, ActivationMethod> table;
	
	static {
		ASWRITTEN = new ActivationMethod("As written in the spectrum or CID if no info", "as written in the spectrum or CID if no info");
		CID = new ActivationMethod("CID", "collision-induced dissociation");
		ETD = new ActivationMethod("ETD", "electron transfer dissociation").electronBased();
		HCD = new ActivationMethod("HCD", "high-energy collision-induced dissociation");
		FUSION = new ActivationMethod("Merge spectra from the same precursor", "Merge spectra from the same precursor");
		PQD = new ActivationMethod("PQD", "pulsed q dissociation");

		table = new LinkedHashMap<String, ActivationMethod>();
		table.put(ASWRITTEN.name, ASWRITTEN);
		table.put(CID.name, CID);
		table.put(ETD.name, ETD);
		table.put(HCD.name, HCD);
		table.put(FUSION.name, FUSION);
		table.put("ETD+SA", ETD);
	}
}
