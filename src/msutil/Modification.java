package msutil;

import java.util.Hashtable;


/**
 * A class representing a modification.
 * @author sangtaekim
 *
 */
public class Modification {
	private final String name;
	private Composition composition;
	private double mass;
	
	private Modification(String name, Composition composition)
	{
		this.name = name;
		this.composition = composition;
	}
	
	private Modification(String name, double mass)
	{
		this.name = name;
		this.composition = null;
		this.mass = mass;
	}
	
	public String getName()	{ return name; }
	public Composition	getComposition()	{ return composition; }
	public float getMass() { return (float)mass; }
	public double getAccurateMass() { return mass; }
	
	public static Modification register(String name, double mass)
	{
		Modification mod = new Modification(name, mass);
		modTable.put(name, mod);
		return mod;
	}

	public static Modification register(String name, Composition composition)
	{
		Modification mod = new Modification(name, composition);
		modTable.put(name, mod);
		return mod;
	}
	
	public static Modification get(String name) { return modTable.get(name); }

	// static member
	private static Modification[] modList = 
	{
		new Modification("Carbamidomethylation", new Composition(2,3,1,1,0)),
		new Modification("Carboxymethylation", new Composition(2,2,2,0,0)),
		new Modification("Oxidation", new Composition(0,0,0,1,0)),
		new Modification("Phosphorylation", new Composition(0,1,0,3,0)),
	};
	
	private static Hashtable<String,Modification> modTable; 
	static {
		modTable = new Hashtable<String, Modification>();
		for(Modification mod : modList)
			modTable.put(mod.getName(), mod);
	}
	
	public static enum Location {
		N_Term,
		C_Term,
		Protein_N_Term,
		Protein_C_Term,
		Anywhere
	}
	
	/**
	 * A class representing the modification instance.
	 * @author sangtaekim
	 *
	 */
	public static class Instance {
		private final Modification mod;
		private final char residue;	// if null, no amino acid specificity
		private Location location;	// N_Term, C_Term, Anywhere
		private boolean isFixedModification = false;
		public Instance(Modification mod, char residue, Location location)
		{
			this.mod = mod;
			this.residue = residue;
			this.location = location;
		}
		public Instance(Modification mod, char residue)
		{
			this(mod, residue, Location.Anywhere);
		}
		public Instance fixedModification()	{ isFixedModification = true; return this; }
		
		public Modification getModification()	{ return mod; }
		public char getResidue()	{ return residue; }
		public Location getLocation()	{ return location; }
		public boolean isFixedModification() { return isFixedModification; }
	}
}
