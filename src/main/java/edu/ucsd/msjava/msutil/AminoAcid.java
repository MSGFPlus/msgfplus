package edu.ucsd.msjava.msutil;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;


/**
 * 
 * @author Sangtae Kim
 *
 */
public class AminoAcid extends Matter {

	// this is recommended for Serializable objects
	static final private long serialVersionUID = 1L;

	private double mass;
	private int nominalMass;
	private char residue;	// 1 letter code for standard amino acid
	private String name;
	private float probability = 0.05f;
	private Composition composition;

	/**
	 * Constructor.
	 * @param residue single letter identifier.
	 * @param name full name of the amino acid.
	 * @param id kept for historical reasons.
	 * @param composition CHNOS composition object.
	 */
	protected AminoAcid(char residue, String name, Composition composition) {
		this.mass = composition.getAccurateMass();
		this.nominalMass = composition.getNominalMass();
		this.residue  = residue;
		this.name = name;
		this.composition = composition;
	}

	/**
	 * Constructor. Generates a custom amino acid.
	 * @param residue
	 * @param mass
	 */
	protected AminoAcid(char residue, String name, double mass) {
		this.mass = mass;
		this.nominalMass = Math.round(Constants.INTEGER_MASS_SCALER*(float)mass);
		this.residue = residue;
		this.name = name;
	}
	
	/**
	 * Builder. Set probability and returns this object.
	 * @return this object.
	 */
	public AminoAcid setProbability(float probability)		{ this.probability = probability; return this; }

	/**
	 * Standard string representation of this object. Output the single letter
	 * representation.
	 * @return the single letter code for this amino acid.
	 */
	public String toString() {
		return String.valueOf(residue) + ": " + String.format("%.2f", mass); 
	}

	/**
	 * Quick way to tell whether this object is modified.
	 * @return false if this is not modified.
	 */
	public boolean isModified()          { return false; }
	
	
	/**
	 * Quick way to tell the number of variable modifications applied to this amino acid.
	 * @return the number of variable modifications applied to this amino acid.
	 */
	public int getNumVariableMods()		{ return 0; }
	
	/**
	 * Tell whether this object is associated with a terminal-specific modification
	 * @return false if this is not associated with terminal-specific modification
	 */
	public boolean hasTerminalVariableMod()
	{
		return false;
	}
	
	/**
	 * Tell whether this object is associated with a residue-specific modification
	 * @return false if this is not associated with residue-specific modification
	 */
	public boolean hasAnywhereVariableMod()
	{
		return false;
	}

	//  accessor methods
	/**
	 * Gets the mass of this amino acid. This is the mono isotopic mass.
	 * @return the mass of this amino acid
	 */
	@Override
	public float getMass()               { return (float)mass; }

	/**
	 * Gets the mass of this amino acid as double precision. This is the mono isotopic mass.
	 * @return the mass of this amino acid (double precision)
	 */
	@Override
	public double getAccurateMass()	{ return mass; }

	/**
	 * Gets the nominal mass of this object. 
	 * @return nominal mass of this object.
	 */
	@Override
	public int getNominalMass()
	{
		return nominalMass;
	}

	/**
	 * Gets the probability of this amino acid. Currently set as 1/20, uniformly.
	 * @return the probability of this amino acid.
	 */
	public float getProbability()
	{
		return probability;
	}

	//  // prohibited
	//  @Override
	//  public void add(AminoAcid other) {
	//	  assert(false);
	//  }

	@Override
	public boolean equals(Object obj) {
		if(!(obj instanceof AminoAcid))
			return false;
		AminoAcid aa = (AminoAcid)obj;
		return this == aa;
	}

	/**
	 * Gets the representation of the residue as string.
	 * @return the string representing this amino acid.    
	 */
	public String getResidueStr()             { return String.valueOf(residue); }

	/**
	 * Gets the single letter amino acid representation.
	 * @return the single letter amino acid character.    
	 */
	public char getResidue() 		{ return residue; }

	/**
	 * Gets the single letter amino acid representation of the unmodified version of this amino acid.
	 * @return the single letter amino acid character.    
	 */
	public char getUnmodResidue() 	{ return residue; }
	
	/**
	 * Gets the full string.
	 * @return the full name/description of the amino acid.
	 */
	public String getName()              { return name; }

	/**
	 * Gets the composition object for this amino acid.
	 * @return the composition object for this amino acid.
	 */
	public Composition getComposition()  { return composition; }

	// static members
	public static AminoAcid getStandardAminoAcid(char residue) { return residueMap.get(residue); }
	public static AminoAcid[] getStandardAminoAcids()	{ return standardAATable; }

	/**
	 * Returns a modified version of this amino acid (fixed modification).
	 * @param newResidue new residue character representing the modified amino acid
	 * @param mod a modification.
	 * @return a modified amino acid object.
	 */
	public AminoAcid getAAWithFixedModification(Modification mod) {
		String name = mod.getName() + " " + this.getName();
		AminoAcid modAA;
		if(mod.getComposition() == null)
			modAA = getCustomAminoAcid(residue, name, mass+mod.getAccurateMass());
		else
			modAA = getAminoAcid(residue, name, composition.getAddition(mod.getComposition()));
		return modAA;
	}

	public static AminoAcid getCustomAminoAcid(char residue, String name, double mass) {
		AminoAcid standardAA = AminoAcid.getStandardAminoAcid(residue);
		if(standardAA != null && Math.abs(mass-standardAA.getMass()) < 0.001f)
			return standardAA;
		else
			return new AminoAcid(residue, name, mass);
	}
	
	public static AminoAcid getCustomAminoAcid(char residue, float mass) {
		return new AminoAcid(residue, "Custom amino acid", mass);
	}
	
	public static AminoAcid getAminoAcid(char residue, String name, Composition composition) {
		AminoAcid standardAA = AminoAcid.getStandardAminoAcid(residue); 
		if(standardAA != null && composition.getAccurateMass() == standardAA.getAccurateMass())
			return standardAA;
		else
			return new AminoAcid(residue, name, composition);
	}
	
	@Override
	public int hashCode()
	{
		return (int)residue;
	}
	
//	@Override
//	public boolean equals(Object obj)
//	{
//		if(!(obj instanceof AminoAcid))
//			return false;
//		else
//		{
//			AminoAcid otherAA = (AminoAcid)obj;
//			return this.getResidue() == otherAA.getResidue();
//		}
//	}

	private static Hashtable<Character, AminoAcid> residueMap;
	// Static table containing Predefined Amino Acids 
	private static final AminoAcid [] standardAATable = 
	{
		new AminoAcid('G',  "Glycine", new Composition(2,3,1,1,0)),
		new AminoAcid('A',  "Alanine", new Composition(3,5,1,1,0)),
		new AminoAcid('S',  "Serine", new Composition(3,5,1,2,0)),
		new AminoAcid('P',  "Proline", new Composition(5,7,1,1,0)),
		new AminoAcid('V',  "Valine", new Composition(5,9,1,1,0)),
		new AminoAcid('T',  "Threonine", new Composition(4,7,1,2,0)),
		new AminoAcid('C',  "Cystine", new Composition(3,5,1,1,1)),
		new AminoAcid('L',  "Leucine", new Composition(6,11,1,1,0)),
		new AminoAcid('I',  "Isoleucine", new Composition(6,11,1,1,0)),
		new AminoAcid('N',  "Asparagine", new Composition(4,6,2,2,0)),
		new AminoAcid('D',  "Aspartate", new Composition(4,5,1,3,0)),
		new AminoAcid('Q',  "Glutamine", new Composition(5,8,2,2,0)),
		new AminoAcid('K',  "Lysine", new Composition(6,12,2,1,0)),
		new AminoAcid('E',  "Glutamate", new Composition(5,7,1,3,0)),
		new AminoAcid('M',  "Methionine", new Composition(5,9,1,1,1)),
		new AminoAcid('H',  "Histidine", new Composition(6,7,3,1,0)),
		new AminoAcid('F',  "Phenylalanine", new Composition(9,9,1,1,0)),
		new AminoAcid('R',  "Arginine", new Composition(6,12,4,1,0)),
		new AminoAcid('Y',  "Tyrosine", new Composition(9,9,1,2,0)),
		new AminoAcid('W',  "Tryptophan", new Composition(11,10,2,1,0)),
	};

//	public static final AminoAcid N_TERN = new AminoAcid('[',  "N-terminus", new Composition(0,0,0,0,0));
//	public static final AminoAcid C_TERM = new AminoAcid(']',  "C-terminus", new Composition(0,0,0,0,0));
//	public static final AminoAcid PROTEIN_N_TERN = new AminoAcid('{',  "Protein N-terminus", new Composition(0,0,0,0,0));
//	public static final AminoAcid PROTEIN_C_TERM = new AminoAcid('}',  "Protein C-terminus", new Composition(0,0,0,0,0));
//	public static final AminoAcid ANY = new AminoAcid('*',  "C-terminus", new Composition(0,0,0,0,0));
	
	static {
		residueMap = new Hashtable<Character, AminoAcid>();
		for(AminoAcid aa : standardAATable)
			residueMap.put(aa.getResidue(), aa);
	}
	/*
  public static Color getColor(AminoAcid aa)
  {
    int index = aa.getIndex();
    switch(index)
    {
    case 0: return new Color(200,200,200);  
    case 1: return new Color(140,255,140);  
    case 2: return new Color(255,112,66); 
    case 3: return new Color(82,82,82);
    case 4: return new Color(255,140,255);
    case 5: return new Color(184,76,0);
    case 6: return new Color(69,94,69);
    case 7: return new Color(0,76,0);
    case 8: return new Color(255,124,112);
    case 9: return new Color(160,0,66);
    case 10: return new Color(102,0,0);
    case 11: return new Color(71,71,184);
    case 12: return new Color(102,0,100);
    case 13: return new Color(184,160,66);
    case 14: return new Color(112,112,255);
    case 15: return new Color(83,76,82);
    case 16: return new Color(100,100,224);
    case 17: return Color.orange;
    case 18: return new Color(140,112,76);
    case 19: return new Color(79,70,0);
    default: return null;
    }
  }  
	 */
	/**
	 * Get the amino acid of the given integer mass.
	 * @param mass the integer mass
	 * @return the list of amino acids with the mass.
	 */
	public static ArrayList<AminoAcid> getAminoAcids(int mass) {
		if (mass2aa.containsKey(mass)) return mass2aa.get(mass);
		return new ArrayList<AminoAcid>();
	}

	/**
	 * Checks whether the character is an standard amino acid
	 * @param c the character input
	 * @return true if it is part of the standard amino acid set, false otherwise
	 */
	public static boolean isStdAminoAcid(char c) {
		return residueMap.containsKey(c);
	}

	private static HashMap<Integer,ArrayList<AminoAcid>> mass2aa;
	static {
		mass2aa = new HashMap<Integer,ArrayList<AminoAcid>>();
		for (AminoAcid aa : getStandardAminoAcids()) {
			if (!mass2aa.containsKey(aa.getNominalMass())) {
				mass2aa.put(aa.getNominalMass(), new ArrayList<AminoAcid>());
			}
			mass2aa.get(aa.getNominalMass()).add(aa);
		}
	}
}

