package msutil;

// for variable modification
public class ModifiedAminoAcid extends AminoAcid {
	private Modification mod;
	private char unmodResidue;
	private boolean isResidueSpecific;
	
	public ModifiedAminoAcid(AminoAcid aa, Modification mod, char residue)
	{
		super(residue, mod.getName()+" "+aa.getName(), aa.getAccurateMass()+mod.getAccurateMass());
		this.mod = mod;
		this.unmodResidue = aa.getResidue();
		super.setProbability(aa.getProbability());
		if(Character.isUpperCase(unmodResidue))
			isResidueSpecific = true;
	}
	
	public char getUnmodResidue() 	{ return unmodResidue; }
	public Modification getModification()	{ return mod; }
	public boolean isResidueSpecific()	{ return isResidueSpecific; }
	
	@Override
	public String getResidueStr()
	{
		StringBuffer buf = new StringBuffer();
		buf.append(unmodResidue);
		float modMass = mod.getMass();
		if(modMass >= 0)
			buf.append('+');
		buf.append(String.format("%.3f", modMass));
		return buf.toString();
	}
	
	@Override
	public boolean isModified()		{ return true; }
}
