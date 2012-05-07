package msutil;

import java.util.*;

import msgf.IntMassFactory;
import msgf.MassListComparator;
import msgf.NominalMass;
import msgf.Tolerance;
import msgf.IntMassFactory.IntMass;
import msutil.Modification.Location;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

/**
 * 
 * @author Sangtae Kim
 *
 */
public class Peptide extends Sequence<AminoAcid> implements Comparable<Peptide> {

	//this is recommended for Serializable objects
	static final private long serialVersionUID = 1L;
	// maximum length of a peptide
	static final int MAX_LENGTH = 30;

	// fields
	private boolean isModified; // Indicates the peptid has a modified aminoacid

	static final boolean FAIL_WHEN_PEPTIDE_IS_MODIFIED=false; // Fail loudly

	// true if this peptide contains invalid amino acid
	private boolean isInvalid = false;
	
	// true if there's n-term modification
//	private boolean hasNTermMod = false;
//	private float nTermModMass = 0f;
	
	/**
	 * Constructor. Parses sequence string and check for modifications. Not fully implemented!!
	 * Examples: QWSYL   -17QSVL   QSV+2.12QLK-3
	 * @params sequence the sequence in string representation.
	 */
	public Peptide(String sequence, AminoAcidSet aaSet) {
		isModified = false;
		int seqLen = sequence.length();
		int index = 0;
		
		float nTermModMass = 0;
		
		// sequence has an N-term fixed mod
		while(index < seqLen)
		{
			char c = sequence.charAt(index);
			if(c == '-' || c == '+')	// sequence has an N-term mod (e.g. +42ACDEFGR)
			{
				int startIndex = index;
				while(++index < seqLen)
				{
					c = sequence.charAt(index);
					if(!Character.isDigit(c) && c != '.')
						break;
				}
				nTermModMass += Float.parseFloat(sequence.substring(startIndex, index));
			}
			else
				break;
		}
		
		boolean isNTerm = true;
		for(;index < seqLen; index++) {
			char c = sequence.charAt(index);
			assert(Character.isLetter(c)):"Error in string at index "+index;
			float mod=0f;
			if (index+1<seqLen) { // Check for modification (e.g. +17, -12.5)
				char sign=sequence.charAt(index+1);
				while (sign=='-' || sign=='+') { // Modification found
					assert(index+2<seqLen):"Missing value after \""+sign+"\"";
					assert(c>='A' && c<='Z' || c>='a' && c<='z'):"Error in string at index "+index+2;
					int startModIdx=index+2;
					int endModIdx=startModIdx+1;
					// Extends substring to find modification value
					while(endModIdx<seqLen &&
							(sequence.charAt(endModIdx)=='.' ||
									sequence.charAt(endModIdx)>='0' &&  sequence.charAt(endModIdx)<='9')) {
						endModIdx++; // A+76
					}
					float modMass = Float.parseFloat(sequence.substring(startModIdx, endModIdx)); 
					if (sign=='-') modMass*=-1f;
					mod += modMass;
					index=endModIdx-1;
					if(endModIdx < sequence.length())
						sign = sequence.charAt(endModIdx);
					else
						break;
				}
				if(index+4<seqLen && sign == 'p' && sequence.charAt(index+2) == 'h')	// phos
				{
					assert(sequence.charAt(index+3) == 'o');
					assert(sequence.charAt(index+4) == 's');
					mod = 79.966331f;
					index += 4;
				}
				else if(index+4 < seqLen && sign>='a' && sign<='z' && (Character.toUpperCase(sign) == c) && (sequence.charAt(index+2) == '-'))	// mutation or phosphorylation
				{
					assert(sequence.charAt(index+3) == '>');
					char mutatedResidue = sequence.charAt(index+4);
					assert(mutatedResidue>='a' && mutatedResidue<='z');
					c = Character.toUpperCase(mutatedResidue);
					index += 4;
				}
			}
			
			AminoAcid aa;
			if(isNTerm)
			{
				aa = aaSet.getAminoAcid(Location.N_Term, c);
				isNTerm = false;
			}
			else
				aa = aaSet.getAminoAcid(c);
			
			// TODO: how to deal C-term fixed mods
			if(!Character.isUpperCase(c) || aa == null)	// not a valid amino acid
			{
				this.isInvalid = true;
				return;
			}
			if(this.size() == 0)
				mod += nTermModMass;
			
			if (mod==0f) this.add(aa);
			else { // modified
				isModified=true; // Now peptide is modified
				float mass = aa.getMass() + mod;
				AminoAcid modAA = VolatileAminoAcid.getVolatileAminoAcid(mass);
				this.add(modAA);
			}
		}
	}  

	/**
	 * Constructor. Parses sequence string and check for modifications. Not fully implemented!!
	 * Examples: QWSYL   Q-17SVL   QSV+2.12QLK-3
	 * @params sequence the sequence in string representation.
	 */
	public Peptide(String sequence) {
		this(sequence, AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
	}

	/**
	 * Constructor from an ArrayList of AminoAcids.
	 * @params aaArray the array of amino acids.
	 */
	public Peptide(ArrayList<AminoAcid> aaArray) {
		for (AminoAcid aa : aaArray) {
			assert(aa!=null):"Null aminoacid";
			this.add(aa);
		}
	}


	/**
	 * Constructor from an ArrayList of AminoAcids.
	 * @params aaArray the array of amino acids. 
	 * added by Kyowon
	 */
	public Peptide(List<AminoAcid> aaArray) {
		for (AminoAcid aa : aaArray) {
			assert(aa!=null):"Null aminoacid";
			this.add(aa);
		}
	}
	
	/**
	 * Constructor from an array of AminoAcids.
	 * @param aaArray the array of aminoacids.
	 */
	public Peptide(AminoAcid[] aaArray) {
		for(AminoAcid aa : aaArray)   this.add(aa);
	}


	/**
	 * Returns a subpeptide of the specified range (half open intervals).
	 * @param fromIndex the index of the starting subsequence (inclusive).
	 * @param toIndex the end index of the subsequence (exclusive)
	 * @return a subsequence of specified range
	 */
	public Peptide subPeptide(int fromIndex, int toIndex) {
		return (Peptide)super.subSequence(fromIndex, toIndex);
	}  

	/** Set isModified as true and returns this object
	 */
	public Peptide setModified() {
		isModified = true;
		return this;
	}

	/**
	 * Set isModified and returns this object
	 * @param isModified the flag indicating whether this peptide is modified
	 * @return
	 */
	public Peptide setModified(boolean isModified) {
		this.isModified = isModified;
		return this;
	}

	/**
	 * Construct an boolean array representing the spectrum. All positions
	 * with a theoretical peak will be true.
	 * @return
	 */
	public boolean[] getBooleanPeptide() {
		boolean[] boolPeptide = new boolean[this.getNominalMass()+1];
		int mass = 0;
		for(AminoAcid aa : this) {
			mass += aa.getNominalMass();
			boolPeptide[mass] = true;
		}
		return boolPeptide;
	}


	public boolean isGappedPeptideTrue(ArrayList<Integer> gp)
	{
		boolean[] boolPeptide = getBooleanPeptide();
		boolean isTrue = true;
		for(int m : gp)
			if(boolPeptide[m] == false)
				isTrue = boolPeptide[m];
		return isTrue;
	}

//	/**
//	 * Returns whether this has N-term mod.
//	 * @return true if this peptide has N-term modification, false otherwise.;
//	 */
//	public boolean hasNTermMod()
//	{
//		return this.hasNTermMod;
//	}
	
//	/**
//	 * Returns n-terminal modification mass;
//	 * @return n-terminal modification mass;
//	 */
//	public float getNTermModMass() 
//	{
//		return this.nTermModMass;
//	}

//	/**
//	 * Returns n-terminal nominal modification mass;
//	 * @return n-terminal nominal modification mass;
//	 */
//	public int getNominalNTermModMass() 
//	{
//		return NominalMass.toNominalMass(this.nTermModMass);
//	}
	
	/**
	 * Returns whether this peptide contains invalid amino acids.
	 * @return true if this peptide is invalid.
	 */
	public boolean isInvalid()	{ return this.isInvalid; }

	/**
	 * Checks whether the last amino acid is modified.
	 * @return false if not modifies, true otherwise.
	 */
	public boolean isCTermModified() {
		return get(this.size()-1).isModified();
	}


	/**
	 * Checks whether the C terminus is triptic.
	 * @return true if the last amino acid is triptic and not modified, false
	 *         otherwise.
	 */
	public boolean hasTrypticCTerm() {
		AminoAcid cTerm = this.get(this.size()-1);
		if(!isCTermModified() && 
				(cTerm==AminoAcid.getStandardAminoAcid('K')||cTerm==AminoAcid.getStandardAminoAcid('R')))
			return true;
		return false;
	}

	/**
	 * Checks whether this peptide has a cleavage site.
	 * @param enzyme the enzyme used 
	 * @return true if this has a cleavage site and false otherwise.
	 */
	public boolean hasCleavageSite(Enzyme enzyme) {
		AminoAcid target;
		if(enzyme.isCTerm())
			target = this.get(this.size()-1);
		else
			target = this.get(0);
		return enzyme.isCleavable(target);
	}

	/**
	 * Gets the amino acid at the given index. 
	 * @param i the index of the amino acid in the array to retrieve.
	 * @return the amino acid at the given index.
	 */
	public AminoAcid get(int i)
	{
		if(i <= -1) // N-terminal
			return null;
		else if(i >= this.size()) // C-terminal
			return null;  
		return super.get(i);
	}


	/**
	 * Defines the order of different Sequence objects. The order is determined
	 * by the first mass which is greater than the corresponding mass
	 * in other (same index) in the other sequence. If there is a tie, the 
	 * longer Sequence is greater. This is like lexographical ordering.
	 * @param other the other Sequence to compare.
	 * @return 1 if this Sequence is greater than the other Sequence, -1 if this
	 *         is smaller, 0 otherwise.
	 */
	public int compareTo(Peptide other) {
		// funky ordering 
		int minSize = java.lang.Math.min(this.size(), other.size());

		for (int i=0; i<minSize; i++) {
			int r = get(i).compareTo(other.get(i));
			if(r!=0) {
				return r;
			}
		}

		int r = size()-other.size();
		if(r>0) {
			return 1;
		}
		else if(r<0) {
			return -1;
		}
		return 0;
	}

	/**
	 * Matches this peptide against other peptide. If differ only by "I" and "L", returns true, returns false otherwise.
	 * @param pep	peptide matched to
	 * @return	true if equals ignoring I/L difference. false otherwise.
	 */
	public boolean equalsIgnoreIL(Peptide pep)
	{
		if(this.size() != pep.size())
			return false;
		for(int i=0; i<this.size(); i++)
		{
			Composition c1 = this.get(i).getComposition();
			Composition c2 = pep.get(i).getComposition();
			if(!c1.equals(c2))
				return false;
		}
		return true;
	}

	/**
	 * Returns a string represention of this peptide.
	 * @return peptide string.
	 */
	public String toString()
	{
		StringBuffer output = new StringBuffer();
		for (AminoAcid aa : this)
		{
			output.append(aa.getResidueStr());
		}
		return output.toString();
	}

	/**
	 * Converts this into a Sequence of compositions.
	 * @param isPrefix whether to convert into prefixes or suffixes
	 * @return sequence of compositions
	 */

	public Sequence<Composition> toCumulativeCompositionSequence(boolean isPrefix, Composition offset)
	{
		Sequence<Composition> seq = new Sequence<Composition>();
		Composition c = offset;
		for(int i=0; i<this.size(); i++)
		{
			if(isPrefix)
			{
				c = c.getAddition(this.get(i).getComposition());
				seq.add(c);
			}
			else
			{
				c = c.getAddition(this.get(this.size()-1-i).getComposition());
				seq.add(c);
			}
		}
		return seq;
	}

	/**
	 * Convert this peptide into a sequence of compositions. Each amino acid maps to a corresponding composition.
	 * @param pep
	 * @return sequence of compositions.
	 */
	public Sequence<Composition> toCompositionSequence()
	{
		Sequence<Composition> seq = new Sequence<Composition>();
		for(AminoAcid aa : this)
			seq.add(aa.getComposition());
		return seq;
	}

	/**
	 * Convert this peptide into a sequence of compositions in the reversed order. Each amino acid maps to a corresponding composition.
	 * @param pep
	 * @return sequence of compositions.
	 */
	public Sequence<Composition> toReverseCompositionSequence()
	{
		Sequence<Composition> seq = new Sequence<Composition>();
		for(int i=this.size()-1; i>=0; i--)
			seq.add(this.get(i).getComposition());
		return seq;
	}

	/**
	 * Convert this peptide into a sequence of integer masses. 
	 * @param pep
	 * @return sequence of integer masses.
	 */
	public Sequence<IntMass> toPrefixIntMassSequence(IntMassFactory factory)
	{
		Sequence<IntMass> seq = new Sequence<IntMass>();
		for(int i=0; i<this.size(); i++)
			seq.add(factory.getInstance(this.get(i).getMass()));
		return seq;
	}
	
	public Sequence<IntMass> toCumulativeIntMassSequence(boolean isPrefix, IntMassFactory factory)
	{
		Sequence<IntMass> seq = new Sequence<IntMass>();
		float mass = 0;
		for(int i=0; i<this.size(); i++)
		{
			if(isPrefix)
			{
				mass += this.get(i).getMass();
				seq.add(factory.getInstance(mass));
			}
			else
			{
				mass += this.get(this.size()-1-i).getMass();
				seq.add(factory.getInstance(mass));
			}
		}
		return seq;
	}
	
	/**
	 * Convert this peptide into a sequence of integer masses in the reversed order. 
	 * @param pep
	 * @return sequence of integer masses in the reversed order.
	 */
	public Sequence<IntMass> toSuffixIntMassSequence(IntMassFactory factory)
	{
		Sequence<IntMass> seq = new Sequence<IntMass>();
		for(int i=this.size()-1; i>=0; i--)
			seq.add(factory.getInstance(this.get(i).getMass()));
		return seq;
	}
	
	/**
	 * Sums up the masses of the amino acids plus the mass of a water molecule.
	 * @return the mass in Daltons of the mono isotopic masses plus water.
	 */
	public float getParentMass() {
		return getMass() + (float)Composition.H2O;
	}

	/**
	 * Computes the number of symmetric b/y pairs
	 * @param tolerance tolerance
	 * @param isPPM true if tolerance is in PPM, false otherwise
	 * @return the number of symmetric b/y pairs
	 */
	public int getNumSymmetricPeaks(Tolerance tolerance) {
		ArrayList<Composition> bIons = toCumulativeCompositionSequence(true, new Composition(0,1,0,0,0));
		ArrayList<Composition> yIons = toCumulativeCompositionSequence(false, new Composition(0,3,0,1,0));
		MassListComparator<Composition> comparator = new MassListComparator<Composition>(bIons, yIons);

		return comparator.getMatchedList(tolerance).length;
	}

	/**
	 * Computes the number of symmetric b/y pairs. Use nominal masses
	 * @return the number of symmetric b/y pairs
	 */
	public int getNumSymmetricPeaks() 
	{
		int numSymmPeaks = 0;
		HashSet<Integer> bIons = new HashSet<Integer>();
		int bMass = 1;
		for(int i=0; i<this.size(); i++)
		{
			bMass += this.get(i).getNominalMass();
			bIons.add(bMass);
		}
		int yMass = 19;
		for(int i=this.size()-1; i>=0; i--)
		{
			yMass += this.get(i).getNominalMass();
			if(bIons.contains(yMass))
				numSymmPeaks++;
		}
		return numSymmPeaks;
	}

	/*
  public float getAvgMass()
  {
    float sum = 0.f;
    if(modMass != null)
      for(int i=0; i<this.size(); i++)
        sum += this.get(i).getAvgMass() + modMass[i];
    else
      for(int i=0; i<this.size(); i++)
        sum += this.get(i).getAvgMass();

    return sum;
  }  

  public float getMonoMass(int precision)
  {
    float sum = 0.f;
    if(modMass != null)
      for(int i=0; i<this.size(); i++)
        sum += (float)Math.round((this.get(i).getMonoMass() + modMass[i])*precision)/precision;
    else
      for(int i=0; i<this.size(); i++)
        sum += (float)Math.round(this.get(i).getMonoMass()*precision) / precision;

    return sum;
  }
	 */


	/**
	 * Sums up the integer mass of this Sequence.
	 * @return the integer mass.
	 */
	public int getNominalMass() {
		int sum = 0;
		for(AminoAcid aa : this) {
			sum += aa.getNominalMass();
		}
		return sum;
	}

	/**
	 * Sums up the integer mass of this Sequence.
	 * @return the integer mass.
	 */
	public int getIntMassIndex(IntMassFactory factory) {
		int sum = 0;
		for(AminoAcid aa : this) {
			sum += factory.getMassIndex(aa.getMass());
		}
		return sum;
	}
	
	public Composition getComposition() {
		Composition c = new Composition(0);
		for(AminoAcid aa : this)
			c.add(aa.getComposition());
		return c;
	}

	/*
  public int  getNumOfAA(ArrayList<AminoAcid> aaList)
  {
    int num = 0;
    for(AminoAcid aa : this)
      if(aaList.contains(aa))
        num++;
    return num;
  }
	 */

	public float getProbability()
	{
		float prob = 1;
		for(int i=0; i<this.size(); i++)
		{
			AminoAcid aa = this.get(i);
			prob *= aa.getProbability();
		}
		return prob;
	}


	public float getNumber()
	{
		float number = 1;
		AminoAcid aaL = AminoAcid.getStandardAminoAcid('L');
		AminoAcid aaI = AminoAcid.getStandardAminoAcid('I');
		AminoAcid aaQ = AminoAcid.getStandardAminoAcid('Q');
		AminoAcid aaK = AminoAcid.getStandardAminoAcid('K');
		for(int i=0; i<this.size(); i++)
		{
			AminoAcid aa = this.get(i);
			if(aa == aaL || aa == aaI || aa == aaQ || aa == aaK)
				number *= 2;
		}
		return number;
	} 


	/**
	 * Returns an slice of the current sequence with the given coordinates.
	 * @param from the starting index (inclusive).
	 * @param to the ending index (exclusive).
	 * @return a new Sequence object after the slice operation, null if the 
	 *         ranges yield no sequence.
	 */
	public Peptide slice(int from, int to) {
		from = java.lang.Math.max(0, from);
		to = java.lang.Math.min(this.size(), to);

		ArrayList<AminoAcid> aaList = new ArrayList<AminoAcid>();
		for(int i=from; i<to; i++)
			aaList.add(this.get(i));
		if(aaList.size() > 0) {
			return new Peptide(aaList);
		}
		return null;
	}


	/**
	 * Factory function that creates a sequence from a String.
	 * @param seq the String representing a standard sequence.
	 * @return the Sequence object.
	 */
	public static Peptide getSequence(String seq) {
		ArrayList<AminoAcid> aaList = new ArrayList<AminoAcid>();
		int seqLen = seq.length();
		for(int i=0; i<seqLen; i++) {
			aaList.add(AminoAcid.getStandardAminoAcid(seq.charAt(i)));
		}
		return new Peptide(aaList);
	}

	
	/**
	 * This function checks that the peptide agrees with the given set of masses
	 * @param masses the masses
	 * @return true if the correct, false otherwise
	 */
	public boolean isCorrect(ArrayList<Integer> masses) {
		int cumMass = 0;
		int massIndex = 0;
		int targetMass = masses.get(massIndex++);
		for (AminoAcid aa : this) {
			cumMass += aa.getNominalMass();
			if (cumMass < targetMass) {
				continue;  // move to the next mass
			}

			if (cumMass==targetMass) {
				// we got a match
				if (massIndex < masses.size()) 
					targetMass += masses.get(massIndex++);
				else
					// we matched everything
					return true;
			}
			else {
				// no match
				return false;
			}
		}

		if (massIndex==masses.size()) return true;
		return false;
	}
	
	
	/**
	 * This function checks that a given peptide sequence agrees with the mass 
	 * list. The mass list can expand more than one amino acid.
	 * @param sequence the amino acid letters
	 * @param masses the masses
	 * @param aaSet the amino acid alphabet
	 * @return true if the condition is true, false otherwise
	 */
	public static boolean isCorrect(String sequence, ArrayList<Integer> masses, AminoAcidSet aaSet) {
		int cumMass = 0;
		int massIndex = 0;
		int targetMass = masses.get(massIndex++);
		for (int i=0; i<sequence.length(); i++) {
			cumMass += aaSet.getAminoAcid(sequence.charAt(i)).getNominalMass();
			if (cumMass < targetMass) {
				continue;  // move to the next mass
			}

			if (cumMass==targetMass) {
				// we got a match
				if (massIndex < masses.size()) 
					targetMass += masses.get(massIndex++);
				else
					// we matched everything
					return true;
			}
			else {
				// no match
				return false;
			}
		}

		if (massIndex==masses.size()) return true;
		return false;
	}


	/**
	 * This function checks that a given peptide sequence agrees with the mass 
	 * list. The mass list can expand more than one amino acid.
	 * @param sequence the amino acid letters
	 * @param masses the masses
	 * @return true if the condition is true, false otherwise
	 */
	public static boolean isCorrect(String sequence, ArrayList<Integer> masses) {
		return isCorrect(sequence, masses, AminoAcidSet.getStandardAminoAcidSet());
	}
	
	
	public float[] getPRMMasses(boolean isPrefix, float offset)
	{
		if(isModified) // TODO handle modified peptide
			return null;
		float[] masses = new float[this.size()-1];
		float mass = offset;

		for(int i=0; i<this.size()-1; i++)
		{
			if(isPrefix)
				mass += this.get(i).getMass();
			else
				mass += this.get(this.size()-1-i).getMass();
			masses[i] = mass;
		}
		return masses;
	}

	public boolean isModified() {
		return isModified;
	}


	/**
	 * Given a string of standard amino acids, return the mass of this string.
	 * @param peptide the string of standard amino acids, one amino acid per character.
	 * @return the mass of this peptide in Daltons.
	 */
	public static float getMassFromString(String peptide) {
		float cumMass = 0f;
		for(int i = peptide.length(), j = 0; i > 0; i--, j++) {
			cumMass += AminoAcid.getStandardAminoAcid(peptide.charAt(j)).getMass();

		}
		return cumMass;
	}


	/*
  public ArrayList<Peak> getTheoSpec(boolean isPrefix, int offset)
  { 
    return getTheoSpec(isPrefix, offset, PeakProperty.NORMAL);
  }

  public ArrayList<Peak> getTheoSpec(boolean isPrefix, int offset, PeakProperty property)
  {
    return getTheoSpec(isPrefix, offset, property, 1);
  }


  public ArrayList<Peak> getCharge2TheoSpec(boolean isPrefix, int offset, PeakProperty property)
  {
    return getTheoSpec(isPrefix, offset, property, 2);
  }

  public ArrayList<Peak> getTheoSpec(boolean isPrefix, int offset, PeakProperty property, int charge)
  {
    ArrayList<Peak> theoSpec = new ArrayList<Peak>();

    float mass = offset;

    for(int i=0; i<this.size()-1; i++)
    {
      if(isPrefix)
        mass += this.get(i).getMonoMass() + modMass[i];
      else
        mass += this.get(this.size()-1-i).getMonoMass() + modMass[this.size()-1-i];

      theoSpec.add(new Peak(i, mass/charge, 0, 1, property));
    }
    Collections.sort(theoSpec);
    return theoSpec;
  }

  public int getNumOfSymmetricPeaksInt()
  {
    ArrayList<Integer> bPeaks = new ArrayList<Integer>();
    ArrayList<Integer> yPeaks = new ArrayList<Integer>();

    int bMass = 1, yMass = 19;
    for(int i=0; i<this.size()-1; i++)
    {
      bMass += (int)get(i).getMass();
      bPeaks.add(bMass);
      yMass += (int)get(this.size()-1-i).getMass();
      yPeaks.add(yMass);
    }
    int numSymm = 0;
    for(int i=0; i<bPeaks.size(); i++)
    {
      if(bPeaks.get(i) > (this.getIntMonoMass()+18)/2)
        break;
      for(int j=0; j<yPeaks.size(); j++)
      {
        if((int)bPeaks.get(i) == (int)yPeaks.get(j))
          numSymm++;
      }
    }
    return numSymm;
  }

  public int getNumOfSymmetricPeaks()
  {
    ArrayList<Peak> bPeaks = getTheoSpec(true, 1);
    ArrayList<Peak> yPeaks = getTheoSpec(false, 19);

    return new PeakListComparator(bPeaks, yPeaks).getSharedPeakCount();
  }

  public int[] getIntPRM()
  {
    int[] intMass = new int[this.size()+1];
    int mass = 0;
    intMass[0] = mass;
    for(int i=0; i<this.size(); i++)
    {
      mass += (int)this.get(i).getMonoMass();
      intMass[i+1] = mass;
    }
    return intMass;
  }
  public boolean isPRM(int prm)
  {
    int mass = 0;
    if(prm == mass)
      return true;
    for(int i=0; i<this.size(); i++)
    {
      mass += (int)get(i).getMass();
      if(mass > prm)
        break;
      else if(mass == prm)
        return true;
    }
    return false;
  }
  public ArrayList<Modification> getModifications()
  {
    if(!isModified)
      return null;
    ArrayList<Modification> modList = new ArrayList<Modification>();
    for(int i=0; i<modMass.length; i++)
    {
      if(modMass[i] != 0)
      {
        modList.add(new Modification(this.get(i).getResidue(), modMass[i]));
      }
    }
    return modList;
  }
	 */
	public static void main(String[] a) {
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(System.getProperty("user.home")+"/Research/ToolDistribution/Mods.txt");
		Peptide p = new Peptide("+42.011+15.995MDNKTPVTLAK", aaSet);
		System.out.println(p);
		for(AminoAcid aa : p)
			System.out.println(aa.getResidueStr()+" " + aa.getMass());
		System.out.println(p.getMass());
	}

}
