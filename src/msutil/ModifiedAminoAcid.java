package msutil;

// for variable modification
public class ModifiedAminoAcid extends AminoAcid {
	private Modification mod;
	private AminoAcid targetAA;
	private boolean isNTermNonSpecificMod = false;
	
	public ModifiedAminoAcid(AminoAcid targetAA, Modification mod, char residue)
	{
		super(residue, mod.getName()+" "+targetAA.getName(), targetAA.getAccurateMass()+mod.getAccurateMass());
		this.mod = mod;
		this.targetAA = targetAA;
		super.setProbability(targetAA.getProbability());
	}
	
	public ModifiedAminoAcid setNTermNonSpecificMod()
	{
		this.isNTermNonSpecificMod = true;
		return this;
	}
	
	@Override
	public char getUnmodResidue() 	{ return targetAA.getUnmodResidue(); }
	public Modification getModification()	{ return mod; }
	
	@Override
	public String getResidueStr()
	{
		StringBuffer buf = new StringBuffer();
		String massStr;
		float modMass = mod.getMass();
		if(modMass >= 0)
			massStr = "+" + String.format("%.3f", modMass);
		else
			massStr = String.format("%.3f", modMass);
		if(isNTermNonSpecificMod)
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
	public boolean isModified()		{ return true; }
}
