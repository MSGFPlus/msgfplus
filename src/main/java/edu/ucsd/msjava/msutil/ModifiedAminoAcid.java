package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.msutil.Modification.Location;

// for variable modification
public class ModifiedAminoAcid extends AminoAcid {
	private Modification mod;
	private AminoAcid targetAA;
	private boolean isNTermVariableMod = false;
	private boolean isCTermVariableMod = false;
	private boolean hasTerminalVariableMod = false;
	private boolean hasAnywhereVariableMod = false;
	private boolean isFixedModification = false;
	private final int numMods;
	
	public ModifiedAminoAcid(AminoAcid targetAA, Modification.Instance mod, char residue)
	{
		super(residue, mod.getModification().getName()+" "+targetAA.getName(), targetAA.getAccurateMass()+mod.getModification().getAccurateMass());
		this.mod = mod.getModification();
		this.targetAA = targetAA;
		this.hasTerminalVariableMod = targetAA.hasTerminalVariableMod();
		this.hasAnywhereVariableMod = targetAA.hasAnywhereVariableMod();
		super.setProbability(targetAA.getProbability());
		if(mod.isFixedModification())
			this.isFixedModification = mod.isFixedModification();
		else if(mod.getLocation() == Location.Anywhere)
			this.hasAnywhereVariableMod = true;
		else
		{
			this.hasTerminalVariableMod = true;
			if(mod.getLocation() == Location.N_Term || mod.getLocation() == Location.Protein_N_Term)
				isNTermVariableMod = true;
			if(mod.getLocation() == Location.C_Term || mod.getLocation() == Location.Protein_C_Term)
				isCTermVariableMod = true;
		}
		if(this.hasAnywhereVariableMod)
		{
			if(this.hasTerminalVariableMod)
				numMods = 2;
			else
				numMods = 1;
		}
		else
		{
			if(this.hasTerminalVariableMod)
				numMods = 1;
			else
				numMods = 0;
		}
	}
	
	@Override
	public char getUnmodResidue() 	{ return targetAA.getUnmodResidue(); }
	
	public Modification getModification()	{ return mod; }
	
	@Override
	public String getResidueStr()
	{
		if(isFixedModification)
			return String.valueOf(getUnmodResidue());
		StringBuffer buf = new StringBuffer();
		String massStr;
		float modMass = mod.getMass();
		if(modMass >= 0)
			massStr = "+" + String.format("%.3f", modMass);
		else
			massStr = String.format("%.3f", modMass);
		if(isNTermVariableMod)
		{
			buf.append(massStr+targetAA.getResidueStr());
		}
		else
		{
			buf.append(targetAA.getResidueStr()+massStr);
		}
		return buf.toString();
	}
	
	@Override
	public boolean isModified()		
	{ 
		return !isFixedModification; 
	}
	
	@Override
	public boolean hasTerminalVariableMod()
	{
		return this.hasTerminalVariableMod;
	}
	
	@Override
	public boolean hasAnywhereVariableMod()
	{
		return this.hasAnywhereVariableMod;
	}
	
	public boolean isNTermVariableMod()
	{
		return isNTermVariableMod;
	}

	public boolean isCTermVariableMod()
	{
		return isCTermVariableMod;
	}
	
	/**
	 * Quick way to tell the number of variable modifications applied to this amino acid.
	 * @return the number of variable modifications applied to this amino acid.
	 */
	public int getNumVariableMods()		
	{ 
		return numMods;
	}
	
}
