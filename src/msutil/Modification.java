package msutil;

import java.util.Hashtable;

import msgf.NominalMass;


/**
 * A class representing a modification.
 * @author sangtaekim
 *
 */
public class Modification {
	private final String name;
	private final double mass;
	private final int nominalMass;
	private Composition composition;
	
	private Modification(String name, Composition composition)
	{
		this.name = name;
		this.mass = composition.getAccurateMass();
		this.nominalMass = composition.getNominalMass();
		this.composition = composition;
	}
	
	private Modification(String name, double mass)
	{
		this.name = name;
		this.mass = mass;
		this.nominalMass = NominalMass.toNominalMass((float)mass);
	}
	
	public String getName()	{ return name; }
	public float getMass() { return (float)mass; }
	public double getAccurateMass() { return mass; }
	public int getNominalMass() { return nominalMass; }
	public Composition getComposition() { return composition; } 
	
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
		new Modification("NIPCAM", new Composition(5,9,1,1,0)),
		new Modification("Oxidation", new Composition(0,0,0,1,0)),
		new Modification("Phosphorylation", Composition.getMass("HO3P")),
		new Modification("Methylation", new Composition(1,2,0,0,0)),
		new Modification("PyrogluQ", Composition.getMass("H-3N-1")),
		new Modification("PyrogluE", Composition.getMass("H-2O-1")),
		new Modification("Carbamylation", new Composition(1,1,1,1,0)),
		new Modification("Acetylation", new Composition(2,2,0,1,0)),
		new Modification("PyroCarbamidomethyl", Composition.getMass("H-3N-1")),
	};
	
	private static Hashtable<String,Modification> modTable; 
	static {
		modTable = new Hashtable<String, Modification>();
		for(Modification mod : modList)
			modTable.put(mod.getName(), mod);
	}
	
	public static enum Location {
		Anywhere,
		N_Term,
		C_Term,
		Protein_N_Term,
		Protein_C_Term,
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
		public String toString()
		{
			return mod.getName()+" "+residue+" "+location+" "+(isFixedModification ? "Fixed" : "Variable");		
		}
		
	}
}
