package edu.ucsd.msjava.msutil;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import edu.ucsd.msjava.params.ParamObject;
import edu.ucsd.msjava.params.UserParam;


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
	public String getParamDescription()	{ return name; }
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

	public static ActivationMethod getByCV(String cvAccession)
	{
		return cvTable.get(cvAccession);
	}
	
	public static ActivationMethod register(String name, String fullName)
	{
		ActivationMethod m = table.get(name);
		if(m != null)
			return m;	// registration was not successful
		else
		{
			ActivationMethod newMethod = new ActivationMethod(name, fullName);
			table.put(name, newMethod);
			return newMethod;
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
		return registeredActMethods.toArray(new ActivationMethod[0]);
	}
	
	private static HashMap<String, ActivationMethod> table;
	private static HashMap<String, ActivationMethod> cvTable;
	private static ArrayList<ActivationMethod> registeredActMethods;
	
	private static void add(ActivationMethod actMethod)
	{
		if(table.put(actMethod.name, actMethod) == null)
			registeredActMethods.add(actMethod);
	}
	
	// add to the HashMap only
	private static void addAlias(String name, ActivationMethod actMethod)
	{
		table.put(name, actMethod);
	}
	
	// add to the list only
	private static void addToList(ActivationMethod actMethod)
	{
		registeredActMethods.add(actMethod);
	}	
	
	static {
		ASWRITTEN = new ActivationMethod("As written in the spectrum or CID if no info", "as written in the spectrum or CID if no info");
		CID = new ActivationMethod("CID", "collision-induced dissociation");
		ETD = new ActivationMethod("ETD", "electron transfer dissociation").electronBased();
		HCD = new ActivationMethod("HCD", "high-energy collision-induced dissociation");
		FUSION = new ActivationMethod("Merge spectra from the same precursor", "Merge spectra from the same precursor");
		PQD = new ActivationMethod("PQD", "pulsed q dissociation");

		table = new HashMap<String, ActivationMethod>();
		
		registeredActMethods = new ArrayList<ActivationMethod>();
		
		addToList(ASWRITTEN);
		add(CID);
		add(ETD);
		add(HCD);
		addToList(FUSION);
		addAlias("ETD+SA", ETD);
		
		// Parse activation methods defined by a user
		File actMethodFile = new File("params/activationMethods.txt");
		if(actMethodFile.exists())
		{
//			System.out.println("Loading " + actMethodFile.getAbsolutePath());
			ArrayList<String> paramStrs = UserParam.parseFromFile(actMethodFile.getPath(), 2);
			for(String paramStr : paramStrs)
			{
				String[] token = paramStr.split(",");
				String shortName = token[0];
				String fullName = token[1];
				ActivationMethod newMethod = new ActivationMethod(shortName, fullName);
				add(newMethod);
			}
		}
		
		cvTable = new HashMap<String, ActivationMethod>();
		cvTable.put("MS:1000133", CID);
		cvTable.put("MS:1000598", ETD);
		cvTable.put("MS:1000422", HCD);
	}
}
