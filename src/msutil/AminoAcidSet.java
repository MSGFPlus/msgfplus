package msutil;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import parser.BufferedLineReader;

import msutil.Modification.Location;

/**
 * A factory class to instantiate a set of amino acids
 * @author sangtaekim
 *
 */
public class AminoAcidSet implements Iterable<AminoAcid> {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static final AminoAcid[] EMPTY_AA_ARRAY = new AminoAcid[0];

	private HashMap<Location, ArrayList<AminoAcid>> aaListMap;
	
	private static HashMap<Location,Location[]> locMap;
	static {
		locMap = new HashMap<Location,Location[]>();
		locMap.put(Location.Anywhere, new Location[] {Location.Anywhere, Location.N_Term, Location.C_Term, Location.Protein_N_Term, Location.Protein_C_Term});
		locMap.put(Location.N_Term, new Location[] {Location.N_Term, Location.Protein_N_Term});
		locMap.put(Location.C_Term, new Location[] {Location.C_Term, Location.Protein_C_Term});
		locMap.put(Location.Protein_N_Term, new Location[] {Location.Protein_N_Term});
		locMap.put(Location.Protein_C_Term, new Location[] {Location.Protein_C_Term});
	}
	
	// for fast indexing
	private HashMap<Character,AminoAcid> residueMap;	// residue -> aa (residue must be unique)
	private HashMap<AminoAcid, Integer> aa2index;		// aa -> index
	private HashMap<Location, HashMap<Character,AminoAcid[]>> standardResidueAAArrayMap; // std residue -> array of amino acids 
	private HashMap<Location,HashMap<Integer,AminoAcid[]>> nominalMass2aa;	// nominalMass -> array of amino acids
	
	private AminoAcid[] allAminoAcidArr;
	private int maxNumberOfVariableModificationsPerPeptide = 2;
	
	private boolean containsModification;	// true if this contains any variable or terminal (fixed or variable) modification
	private boolean containsNTermModification;	// true if this contains any (fixed or variable) modification specific to N-terminus
	private boolean containsCTermModification;	// true if this contains any (fixed or variable) modification specific to N-terminus
	private boolean containsPhosphorylation;	// true if this contains phosphorylation
	
	private HashSet<Character> modResidueSet = new HashSet<Character>();	// set of symbols used for residues
	private char nextResidue;
	
	// for enzyme
//	private ArrayList<AminoAcid> enzymeAAList;
	private int neighboringAACleavageCredit = 0;
	private int neighboringAACleavagePenalty = 0;
	private int peptideCleavageCredit = 0;
	private int peptideCleavagePenalty = 0;
	private float probCleavageSites = 0;

	AminoAcid lightestAA, heaviestAA;
	
	private AminoAcidSet() // prevents instantiation 
	{
		aaListMap = new HashMap<Location, ArrayList<AminoAcid>>();
		standardResidueAAArrayMap = new HashMap<Location, HashMap<Character,AminoAcid[]>>();
		for(Location location : Location.values())
		{
			aaListMap.put(location, new ArrayList<AminoAcid>());
		}
		nextResidue = 128;
	}	

	/**
	 * Returns the list of amino acids specific to the position.
	 * @return list of intermediate amino acids.
	 */
	public ArrayList<AminoAcid> getAAList(Location location)
	{
		return aaListMap.get(location);
	}
	
	public ArrayList<AminoAcid> getNTermAAList()
	{
		return aaListMap.get(Location.N_Term);
	}

	public ArrayList<AminoAcid> getCTermAAList()
	{
		return aaListMap.get(Location.C_Term);
	}
	
	public ArrayList<AminoAcid> getProtNTermAAList()
	{
		return aaListMap.get(Location.Protein_N_Term);
	}

	public ArrayList<AminoAcid> getProtCTermAAList()
	{
		return aaListMap.get(Location.Protein_N_Term);
	}
	
	/**
	 * Returns the iterator of anywhere amino acids 
	 */
	public Iterator<AminoAcid> iterator() {
		return aaListMap.get(Location.Anywhere).iterator();
	}
	
	/**
	 * Returns the size of amino acid depending on the location.
	 * @param location amino acid location
	 * @return
	 */
	public int size(Location location) {
		return aaListMap.get(location).size();
	}
	
	/**
	 * Reterns the size of anywhere amino acids
	 * @return the size of anywhere amino acids
	 */
	public int size()
	{
		return aaListMap.get(Location.Anywhere).size();
	}
	
	/**
	 * Retrieve an array of amino acids given the specific standard residue. 
	 * @param location amino acid location
	 * @param standardAAResidue the standard residue to look up
	 * @return the array of amino acids or an empty array otherwise
	 */
	public AminoAcid[] getAminoAcids(Location location, char standardAAResidue) {
		AminoAcid[] matches = standardResidueAAArrayMap.get(location).get(standardAAResidue);
		if(matches != null)
			return matches;
		else
			return EMPTY_AA_ARRAY;
	}

	/**
	 * Retrieve an array of amino acids given the specific nominal mass. 
	 * @param location amino acid location
	 * @param nominalMass nominal mass to look up
	 * @return the array of amino acids or an empty list otherwise
	 */
	public AminoAcid[] getAminoAcids(Location location, int nominalMass) {
		AminoAcid[] matches = nominalMass2aa.get(location).get(nominalMass);
		if (matches != null) return matches;
		return EMPTY_AA_ARRAY;
	}
	
	/**
	 * Retrieve an array of amino acids given the specific nominal mass.
	 * @param nominalMass the mass to look up
	 * @return the array of amino acids or an empty list otherwise
	 */
	public AminoAcid[] getAminoAcids(int nominalMass) {
		return getAminoAcids(Location.Anywhere, nominalMass);
	}
	
	/**
	 * Checks whether a residue belongs to this amino acid set
	 * @param residue a residue
	 * @return	true if residue belongs to the amino acid set
	 */
	public boolean contains(char residue) {
		return residueMap.containsKey(residue);  
	}

	/**
	 * Get the amino acid mass of the residue.
	 * @param residue the amino acid mass. Use upper case for standard aa (convention).  
	 *                this method is case sensitive.
	 * @return the amino acid object. null if no aa corresponding to the residue
	 */
	public AminoAcid getAminoAcid(Location location, char residue)	
	{ 
		AminoAcid[] aaArr = getAminoAcids(location, residue);
		for(AminoAcid aa : aaArr)
			if(!aa.isModified())
				return aa;
		return null; 
	}

	/**
	 * Get the amino acid mass of the residue.
	 * @param residue the amino acid mass. Use upper case for standard aa (convention).  
	 *                this method is case sensitive.
	 * @return the amino acid object. null if no aa corresponding to the residue
	 */
	public AminoAcid getAminoAcid(char residue)	{ return residueMap.get(residue); }
	
	/**
	 * Set the number of allowable variable modifications per peptide
	 * @param maxNumberOfVariableModificationsPerPeptide the number of allowable variable modifications per peptide
	 */
	public void setMaxNumberOfVariableModificationsPerPeptide(int maxNumberOfVariableModificationsPerPeptide)
	{
		this.maxNumberOfVariableModificationsPerPeptide = maxNumberOfVariableModificationsPerPeptide;
	}

	/**
	 * Get the number of allowable variable modifications per peptide
	 * @return the number of allowable variable modifications per peptide
	 */
	public int getMaxNumberOfVariableModificationsPerPeptide()
	{
		return this.maxNumberOfVariableModificationsPerPeptide;
	}
	
	/**
	 * Get all amino acids for all locations.
	 * @return an array of all amino acids.
	 */
	public AminoAcid[] getAllAminoAcidArr()
	{
		return this.allAminoAcidArr;
	}
	
	/**
	 * Get the amino acid corresponding to the index
	 * @param index amino acid index
	 * @return amino acid object
	 */
	public AminoAcid getAminoAcid(int index)
	{
		return allAminoAcidArr[index];
	}
	
	/**
	 * Get the index of the aa
	 * @param aa amino acid
	 * @return the index of aa. null if aa does not belong to this amino acid set
	 */
	public int getIndex(AminoAcid aa)
	{
		Integer index = aa2index.get(aa);
		if(index == null)
			index = -1;
		return index;
	}
	
	/**
	 * Get the peptide corresponding to the string sequence. 
	 * @param sequence sequence of the peptide.
	 * @return peptide object of the sequence
	 */
	public Peptide getPeptide(String sequence)
	{
		boolean isModified = false;
		ArrayList<AminoAcid> aaArray = new ArrayList<AminoAcid>();
		for(int i=0; i<sequence.length(); i++)
		{
			char residue = sequence.charAt(i);
			AminoAcid aa = this.getAminoAcid(residue);
			assert(aa != null): sequence + ": " + residue + " is null!";
			if(aa.isModified())
				isModified = true;
			aaArray.add(aa);
		}
		Peptide pep = new Peptide(aaArray);
		pep.setModified(isModified);

		return pep;
	}	

	public int getMaxNominalMass() { return this.heaviestAA.getNominalMass(); }
	public int getMinNominalMass() { return this.lightestAA.getNominalMass(); }
	
	public AminoAcid getLightestAA()	{ return this.lightestAA; }
	public AminoAcid getHeaviestAA()	{ return this.heaviestAA; }

	public boolean containsModification()		{ return this.containsModification; }
	public boolean containsNTermModification()	{ return this.containsNTermModification; }
	public boolean containsCTermModification()	{ return this.containsCTermModification; }
	public boolean containsPhosphorylation()	{ return this.containsPhosphorylation; }
	
	public char getMaxResidue()	{ return nextResidue; }
	
	public void registerEnzyme(Enzyme enzyme)
	{
		if(enzyme == null || enzyme.getResidues() == null || 
				enzyme.getPeptideCleavageEfficiency() == 0 || enzyme.getNeighboringAACleavageEffiency() == 0)
			return;
		
		probCleavageSites = 0;
		for(char residue : enzyme.getResidues())
		{
			AminoAcid aa = this.getAminoAcid(residue);
			if(aa == null)
			{
				System.err.println("Invalid Enzyme cleavage site: " + residue);
				System.exit(-1);
			}
			probCleavageSites += aa.getProbability();
		}
		
		if(probCleavageSites == 0 || probCleavageSites == 1)
		{
			System.err.println("Probability of enzyme residues must be in (0,1)!");
			System.exit(-1);
		}
		
		float peptideCleavageEfficiency = enzyme.getPeptideCleavageEfficiency();
		float neighboringAACleavageEfficiency = enzyme.getNeighboringAACleavageEffiency();
		
		peptideCleavageCredit = (int)Math.round(Math.log(peptideCleavageEfficiency/probCleavageSites));
		peptideCleavagePenalty = (int)Math.round(Math.log((1-peptideCleavageEfficiency)/(1-probCleavageSites)));
		neighboringAACleavageCredit = (int)Math.round(Math.log(neighboringAACleavageEfficiency/probCleavageSites));
		neighboringAACleavagePenalty = (int)Math.round(Math.log((1-neighboringAACleavageEfficiency)/(1-probCleavageSites)));
	}
	
	public int getNeighboringAACleavageCredit()
	{
		return neighboringAACleavageCredit;
	}

	public int getNeighboringAACleavagePenalty()
	{
		return neighboringAACleavagePenalty;
	}
	
	public int getPeptideCleavageCredit()
	{
		return peptideCleavageCredit;
	}
	
	public int getPeptideCleavagePenalty()
	{
		return peptideCleavagePenalty;
	}
	
	public float getProbCleavageSites()
	{
		return probCleavageSites;
	}
	
	public void printAASet()
	{
		System.out.println("NumMods: " + this.getMaxNumberOfVariableModificationsPerPeptide());
		for(Location location : Location.values())
		{
			ArrayList<AminoAcid> aaList = this.getAAList(location);
			System.out.println(location+"\t"+aaList.size());
			for(AminoAcid aa : aaList)
				System.out.println(aa.getResidueStr()+(aa.isModified() ? "*" : "")+"\t"+(int)aa.getResidue()+"\t"+aa.getNominalMass()+"\t"+aa.getMass()+"\t"+aa.getProbability());
		}
	}
	
	// private members to build an amino acid set
	private void addAminoAcid(AminoAcid aa)
	{
		addAminoAcid(aa, Location.Anywhere);
	}

	// private members
	private void addAminoAcid(AminoAcid aa, Location location)
	{
		for(Location loc : locMap.get(location))
		{
			// Debug
//			if(aa.isModified())
//				System.out.println("Debug");
			aaListMap.get(loc).add(aa);
		}
	}
	
	private List<Modification.Instance> modifications;
	
	private void applyModifications(ArrayList<Modification.Instance> mods)
	{
		this.modifications = mods;
		
		// partition modification instances into different types
		HashMap<Modification.Location,ArrayList<Modification.Instance>> fixedMods = new HashMap<Modification.Location,ArrayList<Modification.Instance>>();
		HashMap<Modification.Location,ArrayList<Modification.Instance>> variableMods = new HashMap<Modification.Location,ArrayList<Modification.Instance>>();
		for(Location location : Modification.Location.values())
		{
			fixedMods.put(location, new ArrayList<Modification.Instance>());
			variableMods.put(location, new ArrayList<Modification.Instance>());
		}
		for(Modification.Instance mod : mods)
		{
			if(mod.isFixedModification())
				fixedMods.get(mod.getLocation()).add(mod);
			else
				variableMods.get(mod.getLocation()).add(mod);
		}

		Location[] locArr = new Location[] {
				Location.Anywhere,
				Location.N_Term,
				Location.C_Term,
				Location.Protein_N_Term,
				Location.Protein_C_Term,
		};
		
		// Fixed modifications
		for(Location loc : locArr)
			applyFixedMods(fixedMods, loc);

		// Variable modifications 
		for(Location loc : locArr)
			addVariableMods(variableMods, loc);
		
		// setup containsNTermModification and containsCTermModification
		for(Modification.Instance mod : mods)
		{
			Location location = mod.getLocation();
			if(!containsNTermModification && (location == Location.N_Term || location == Location.Protein_N_Term))
				this.containsNTermModification = true;
			if(!containsNTermModification && (location == Location.C_Term || location == Location.Protein_C_Term))
				this.containsCTermModification = true;
			if(location != Location.Anywhere || !mod.isFixedModification())
				this.containsModification = true;
			if(mod.getModification().getName().equalsIgnoreCase("phosphorylation"))
				this.containsPhosphorylation = true;
		}
	}

	private void applyFixedMods(HashMap<Modification.Location,ArrayList<Modification.Instance>> fixedMods, Location location)
	{
		for(Modification.Instance mod : fixedMods.get(location))
		{
			// residue-specific
			char residue = mod.getResidue();
			if(residue == '*')
				continue;
			
			ArrayList<AminoAcid> oldAAList = this.getAAList(location);
			ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
			
			for(AminoAcid aa : oldAAList)
			{
				if(aa.getUnmodResidue() != residue)
					newAAList.add(aa);
				else
				{
					if(location == Location.Anywhere)
						newAAList.add(aa.getAAWithFixedModification(mod.getModification()));	// replace with a new amino acid
					else
					{
						char modResidue = this.getModifiedResidue(aa.getUnmodResidue());
						// make a new amino acid and add
						ModifiedAminoAcid modAA = new ModifiedAminoAcid(aa, mod, modResidue);
						newAAList.add(modAA);
					}
				}
			}
			
			for(Location loc : locMap.get(location))
				aaListMap.put(loc, new ArrayList<AminoAcid>(newAAList));
		}
		
		// any residue
		for(Modification.Instance mod : fixedMods.get(location))
		{
			char residue = mod.getResidue();
			if(residue != '*')
				continue;
			ArrayList<AminoAcid> oldAAList = this.getAAList(location);
			ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
			
			for(AminoAcid aa : oldAAList)
			{
				if(location == Location.Anywhere)
					newAAList.add(aa.getAAWithFixedModification(mod.getModification()));
				else
				{
					char modResidue = this.getModifiedResidue(aa.getUnmodResidue());
					ModifiedAminoAcid modAA = new ModifiedAminoAcid(aa, mod, modResidue);
					newAAList.add(modAA);
				}
			}
			
			for(Location loc : locMap.get(location))
				aaListMap.put(loc, new ArrayList<AminoAcid>(newAAList));
		}
	}
	
	private void addVariableMods(HashMap<Modification.Location,ArrayList<Modification.Instance>> variableMods, Location location)
	{
		// residue-specific
		for(Location loc : locMap.get(location))
		{
			ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
			ArrayList<AminoAcid> oldAAList = this.getAAList(loc);
			for(AminoAcid targetAA : oldAAList)
			{
				for(Modification.Instance mod : variableMods.get(location))
				{
					char residue = mod.getResidue();
					if(residue == '*')
						continue;
					if(targetAA.getUnmodResidue() == residue)
					{
						if(targetAA.isModified())
						{
							if(mod.getLocation() == Location.Anywhere && targetAA.hasAnywhereVariableMod())	// residue mod
								continue;
							if(mod.getLocation() != Location.Anywhere && targetAA.hasTerminalVariableMod())	// residue mod
								continue;
						}
						char modResidue = this.getModifiedResidue(targetAA.getUnmodResidue());
						ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod, modResidue);
						newAAList.add(modAA);
					}
				}
			}
			for(AminoAcid newAA : newAAList)
				aaListMap.get(loc).add(newAA);
		}
		
		// any residue
		for(Location loc : locMap.get(location))
		{
			ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
			ArrayList<AminoAcid> oldAAList = this.getAAList(loc);
			for(AminoAcid targetAA : oldAAList)
			{
				for(Modification.Instance mod : variableMods.get(location))
				{
					char residue = mod.getResidue();
					if(residue != '*')
						continue;
//					if(location == Location.Anywhere)
//					{
//						System.err.println("Invalid modification: " + mod);
//						System.exit(-1);
//					}
					if(targetAA.isModified())
					{
						if(mod.getLocation() == Location.Anywhere && targetAA.hasAnywhereVariableMod())	// residue mod
							continue;
						if(mod.getLocation() != Location.Anywhere && targetAA.hasTerminalVariableMod())	// residue mod
							continue;
					}
					char modResidue = this.getModifiedResidue(targetAA.getUnmodResidue());
					ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod, modResidue);
					newAAList.add(modAA);
				}
			}
			for(AminoAcid newAA : newAAList)
				aaListMap.get(loc).add(newAA);
		}
	}
	
	private AminoAcidSet finalizeSet()
	{
		standardResidueAAArrayMap = new HashMap<Location, HashMap<Character,AminoAcid[]>>();
		nominalMass2aa = new HashMap<Location,HashMap<Integer,AminoAcid[]>>();
		for(Location location : Location.values())
		{
			standardResidueAAArrayMap.put(location, new HashMap<Character,AminoAcid[]>());
			nominalMass2aa.put(location, new HashMap<Integer,AminoAcid[]>());
		}
		
		// add all amino acids to aaList
		HashSet<AminoAcid> allAASet = new HashSet<AminoAcid>();
		for(Location location : aaListMap.keySet())
		{
			for(AminoAcid aa : aaListMap.get(location))
				allAASet.add(aa);
		}
		
		this.allAminoAcidArr = allAASet.toArray(EMPTY_AA_ARRAY);
		Arrays.sort(allAminoAcidArr);

		// assign index, heaviest and lightest aa
		double minMass = Double.MAX_VALUE;
		int lightIndex = -1;
		double maxMass = Double.MIN_VALUE;
		int heavyIndex = -1;
		aa2index = new HashMap<AminoAcid, Integer>() ;		// aa -> index
		for(int i=0; i<allAminoAcidArr.length; i++)
		{
			aa2index.put(allAminoAcidArr[i], i);
			double mass = allAminoAcidArr[i].getAccurateMass();
			if(mass < minMass)
			{
				lightIndex = i;
				minMass = mass;
			}
			if(mass > maxMass)
			{
				heavyIndex = i;
				maxMass = mass;
			}
		}
		this.heaviestAA = allAminoAcidArr[heavyIndex];
		this.lightestAA = allAminoAcidArr[lightIndex];
		
		// initialize aaList and residueMap
		residueMap = new HashMap<Character, AminoAcid>();
		
		for(AminoAcid aa : allAminoAcidArr)
		{
			assert(residueMap.get(aa.getResidue()) == null): aa.getResidue()+" already exists!";
			residueMap.put(aa.getResidue(), aa);
		}			

		for(Location location : Location.values())
		{
			HashMap<Integer,ArrayList<AminoAcid>> mass2aaList = new HashMap<Integer,ArrayList<AminoAcid>>();
			HashMap<Character,LinkedList<AminoAcid>> stdResidue2aaList = new HashMap<Character,LinkedList<AminoAcid>>();
			
			for (AminoAcid aa : this.getAAList(location)) 
			{
				int thisMass = aa.getNominalMass();
				if (!mass2aaList.containsKey(thisMass)) {
					mass2aaList.put(thisMass, new ArrayList<AminoAcid>());
				}
				mass2aaList.get(thisMass).add(aa);
				
				char stdResidue = aa.getUnmodResidue();
				LinkedList<AminoAcid> aaList = stdResidue2aaList.get(stdResidue);
				if(aaList == null)
					aaList = new LinkedList<AminoAcid>();
				if(!aa.isModified())
					aaList.addFirst(aa);	// unmodified residue is at first
				else
					aaList.addLast(aa);
				stdResidue2aaList.put(stdResidue, aaList);
			}			

			// convert the array back to real arrays
			HashMap<Integer,AminoAcid[]> mass2aaArray = new HashMap<Integer,AminoAcid[]>();
			for (int mass : mass2aaList.keySet()) {
				mass2aaArray.put(mass, mass2aaList.get(mass).toArray(new AminoAcid[0]));
			}
			
			HashMap<Character,AminoAcid[]> stdResidue2aaArray = new HashMap<Character,AminoAcid[]>();
			for(char residue : stdResidue2aaList.keySet()) 
				stdResidue2aaArray.put(residue, stdResidue2aaList.get(residue).toArray(new AminoAcid[0]));

			this.nominalMass2aa.put(location, mass2aaArray);
			this.standardResidueAAArrayMap.put(location, stdResidue2aaArray);
		}
		
		return this;
	}	

	// static members
	private static AminoAcidSet standardAASet = null;
	private static AminoAcidSet standardAASetWithCarbamidomethylatedCys = null;
	private static AminoAcidSet standardAASetWithCarboxyomethylatedCys = null;
	private static AminoAcidSet standardAASetWithCarbamidomethylatedCysWithTerm = null;
	
	public static AminoAcidSet getAminoAcidSetFromModFile(String fileName)
	{
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		int numMods = 2;
		
		// parse modifications
		ArrayList<Modification.Instance> mods = new ArrayList<Modification.Instance>();
		String s;
		int lineNum = 0;
		while((s=in.readLine()) != null)
		{
			lineNum++;
			String[] tokenArr = s.split("#");
			if(tokenArr.length == 0)
				continue;
			s = tokenArr[0].trim();
			if(s.length() == 0)
				continue;
			else if(s.startsWith("NumMods="))
			{
				try {
					numMods = Integer.parseInt(s.split("=")[1].trim());
				} catch (NumberFormatException e)
				{
					System.err.println(fileName + ": Illegal NumMods option at line " + lineNum + ": " + s);
					e.printStackTrace();
					System.exit(-1);
				}
			}
			else
			{
				String[] token = s.split(",");
				if(token.length != 5)
					continue;

				// Mass or Composition
				double modMass = 0;
				String compStr = token[0].trim();
				if(compStr.matches("(C[+-]?\\d*)*(H[+-]?\\d*)*(N[+-]?\\d*)*(O[+-]?\\d*)*(S[+-]?\\d*)*(P[+-]?\\d*)*"))
				{
					modMass = Composition.getMass(compStr);
				}
				else
				{
					try {
						modMass = Double.parseDouble(compStr);
					}
					catch (NumberFormatException e)
					{
						System.err.println(fileName + ": AminoAcidSet: Illegal Mass/Composition at line " + lineNum + ": " + s);
						e.printStackTrace();
						System.exit(-1);
					}
				}
				
				// Residues
				String residueStr = token[1].trim();
				boolean isResidueStrLegitimate = true;
				if(!residueStr.equals("*"))
				{
					if(residueStr.length() > 0)
					{
						for(int i=0; i<residueStr.length(); i++)
						{
							if(!AminoAcid.isStdAminoAcid(residueStr.charAt(i)))
							{
								isResidueStrLegitimate = false;
								break;
							}
						}
					}
					else
						isResidueStrLegitimate = false;
				}
				if(!isResidueStrLegitimate)
				{
					System.err.println(fileName + ": AminoAcidSet: Illegal Residue(s) at line " + lineNum + ": " + s);
					System.exit(-1);
				}
				
				// isFixedModification
				boolean isFixedModification = false;
				if(token[2].trim().equalsIgnoreCase("fix"))
					isFixedModification = true;
				else if(token[2].trim().equalsIgnoreCase("opt"))
					isFixedModification = false;
				else
				{
					System.err.println(fileName + ": AminoAcidSet: Modification must be either fix or opt at line " + lineNum + ": " + s);
					System.exit(-1);
				}
					
				// Location
				Modification.Location location = null;
				String locStr = token[3].split("\\s+")[0].trim();
				if(locStr.equalsIgnoreCase("any"))
					location = Modification.Location.Anywhere;
				else if(locStr.equalsIgnoreCase("N-Term") || locStr.equalsIgnoreCase("NTerm"))
					location = Modification.Location.N_Term;
				else if(locStr.equalsIgnoreCase("C-Term") || locStr.equalsIgnoreCase("CTerm"))
					location = Modification.Location.C_Term;
				else if(locStr.equalsIgnoreCase("Prot-N-Term") || locStr.equalsIgnoreCase("ProtNTerm"))
					location = Modification.Location.Protein_N_Term;
				else if(locStr.equalsIgnoreCase("Prot-C-Term") || locStr.equalsIgnoreCase("ProtCTerm"))
					location = Modification.Location.Protein_C_Term;
				else
				{
					System.err.println(fileName + ": AminoAcidSet: Illegal Location at line " + lineNum + ": " + s);
					System.exit(-1);
				}

				String name = token[4].split("\\s+")[0].trim();
				
				Modification mod = Modification.register(name, modMass);
				
				for(int i=0; i<residueStr.length(); i++)
				{
					char residue = residueStr.charAt(i);
					Modification.Instance modIns = new Modification.Instance(mod, residue, location);
					if(isFixedModification)
						modIns.fixedModification();
					mods.add(modIns);
				}				
			}
		}
		
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(mods);
		aaSet.setMaxNumberOfVariableModificationsPerPeptide(numMods);
		return aaSet;
	}
	
	public static AminoAcidSet getAminoAcidSetFromXMLFile(String fileName)
	{
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		int numMods = 2;
		
		// memorize keywords
		String numModsKey = "<parameter name=\"ptm.mods\">";
		String cysKey = "<parameter name=\"cysteine_protease.cysteine\">";
		String oxidationKey = "<parameter name=\"ptm.OXIDATION\">on</parameter>";
		String lysMetKey = "<parameter name=\"ptm.LYSINE_METHYLATION\">on</parameter>";
		String pyrogluKey = "<parameter name=\"ptm.PYROGLUTAMATE_FORMATION\">on</parameter>";
		String phosphoKey = "<parameter name=\"ptm.PHOSPHORYLATION\">on</parameter>";
		String ntermCarbamyKey = "<parameter name=\"ptm.NTERM_CARBAMYLATION\">on</parameter>";
		String ntermAcetylKey = "<parameter name=\"ptm.NTERM_ACETYLATION\">on</parameter>";
		String ptmKey = "<parameter name=\"ptm.custom_PTM\">";
		String closeKey = "</parameter>";
		
		// parse modifications
		ArrayList<Modification.Instance> mods = new ArrayList<Modification.Instance>();
		String s;
		int lineNum = 0;
		while((s=in.readLine()) != null)
		{
			lineNum++;
			if(s.startsWith(numModsKey))
			{
				try {
					String value = s.substring(numModsKey.length(), s.lastIndexOf(closeKey));
					numMods = Integer.parseInt(value);
				} catch (NumberFormatException e)
				{
					System.err.println(fileName + ": Illegal ptm.mods option at line " + lineNum + ": " + s);
					e.printStackTrace();
					System.exit(-1);
				}
			}
			else if(s.startsWith(cysKey))
			{
				String value = s.substring(cysKey.length(), s.lastIndexOf(closeKey));
				if(value.equals("c57"))
				{
					char residue = 'C';
					Modification mod = Modification.get("Carbamidomethylation");
					Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere).fixedModification();
					mods.add(modIns);
				}
				else if(value.equals("c58"))
				{
					char residue = 'C';
					Modification mod = Modification.get("Carboxymethylation");
					Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere).fixedModification();
					mods.add(modIns);
				}
				else if(value.equals("c99"))
				{
					char residue = 'C';
					Modification mod = Modification.get("NIPCAM");
					Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere).fixedModification();
					mods.add(modIns);
				}
				else if(value.equals("None"))
				{
					// do nothing
				}
				else
				{
					System.err.println(fileName + ": Illegal Cycteine protecting group at line " + lineNum + ": " + s);
					System.exit(-1);
				}
			}
			else if(s.startsWith(ptmKey))	// custom PTM
			{
				String value = s.substring(ptmKey.length(), s.lastIndexOf(closeKey));
				String[] token = value.split(",");
				
				if(token.length != 3)
				{
					System.err.println(fileName + ": Illegal custom ptm option at line " + lineNum + ": " + s);
					System.exit(-1);
				}

				// Mass
				double modMass = 0;
				try {
					modMass = Double.parseDouble(token[0]);
				}
				catch (NumberFormatException e)
				{
					System.err.println(fileName + ": AminoAcidSet: Illegal Mass at line " + lineNum + ": " + s);
					e.printStackTrace();
					System.exit(-1);
				}
				
				// Residues
				String residueStr = token[1];
				boolean isResidueStrLegitimate = true;
				if(!residueStr.equals("*"))
				{
					if(residueStr.length() > 0)
					{
						for(int i=0; i<residueStr.length(); i++)
						{
							if(!AminoAcid.isStdAminoAcid(residueStr.charAt(i)))
							{
								isResidueStrLegitimate = false;
								break;
							}
						}
					}
					else
						isResidueStrLegitimate = false;
				}
				if(!isResidueStrLegitimate)
				{
					System.err.println(fileName + ": AminoAcidSet: Illegal Residue(s) at line " + lineNum + ": " + s);
					System.exit(-1);
				}
				
				boolean isFixedModification = false;
				Modification.Location location = null;
				
				// Type
				String type = token[2];
				if(type.equals("fix"))
				{
					isFixedModification = true;
					location = Location.Anywhere;
				}
				else if(type.equals("opt"))
				{
					isFixedModification = false;
					location = Location.Anywhere;
				}
				else if(type.equals("opt_nterm"))
				{
					isFixedModification = false;
					location = Location.N_Term;
				}
				else if(type.equals("fix_nterm"))
				{
					isFixedModification = true;
					location = Location.N_Term;
				}
				else if(type.equals("opt_cterm"))
				{
					isFixedModification = false;
					location = Location.C_Term;
				}
				else if(type.equals("fix_cterm"))
				{
					isFixedModification = true;
					location = Location.C_Term;
				}
				else
				{
					System.err.println(fileName + ": AminoAcidSet: Illegal Type(s) at line " + lineNum + ": " + s);
					System.exit(-1);
				}

				String name = residueStr + " " + modMass;
				
				Modification mod = Modification.register(name, modMass);
				
				for(int i=0; i<residueStr.length(); i++)
				{
					char residue = residueStr.charAt(i);
					Modification.Instance modIns = new Modification.Instance(mod, residue, location);
					if(isFixedModification)
						modIns.fixedModification();
					mods.add(modIns);
				}				
			}
			else if(s.startsWith(oxidationKey))	// predefined Oxidation
			{
				String residueStr = "M";
				Modification mod = Modification.get("Oxidation");
				for(int i=0; i<residueStr.length(); i++)
				{
					char residue = residueStr.charAt(i);
					Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere);
					mods.add(modIns);
				}				
			}
			else if(s.startsWith(lysMetKey))	// predefined
			{
				String residueStr = "K";
				Modification mod = Modification.get("Methylation");
				for(int i=0; i<residueStr.length(); i++)
				{
					char residue = residueStr.charAt(i);
					Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere);
					mods.add(modIns);
				}				
			}
			else if(s.startsWith(pyrogluKey))	// predefined
			{
				String residueStr = "Q";
				Modification mod = Modification.get("PyrogluQ");
				for(int i=0; i<residueStr.length(); i++)
				{
					char residue = residueStr.charAt(i);
					Modification.Instance modIns = new Modification.Instance(mod, residue, Location.N_Term);
					mods.add(modIns);
				}				
			}
			else if(s.startsWith(phosphoKey))	// predefined
			{
				String residueStr = "STY";
				Modification mod = Modification.get("Phosphorylation");
				for(int i=0; i<residueStr.length(); i++)
				{
					char residue = residueStr.charAt(i);
					Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere);
					mods.add(modIns);
				}				
			}
			else if(s.startsWith(ntermCarbamyKey))	// predefined
			{
				String residueStr = "*";
				Modification mod = Modification.get("Carbamylation");
				for(int i=0; i<residueStr.length(); i++)
				{
					char residue = residueStr.charAt(i);
					Modification.Instance modIns = new Modification.Instance(mod, residue, Location.N_Term);
					mods.add(modIns);
				}				
			}
			else if(s.startsWith(ntermAcetylKey))	// predefined
			{
				String residueStr = "*";
				Modification mod = Modification.get("Acetylation");
				for(int i=0; i<residueStr.length(); i++)
				{
					char residue = residueStr.charAt(i);
					Modification.Instance modIns = new Modification.Instance(mod, residue, Location.N_Term);
					mods.add(modIns);
				}				
			}
		}
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(mods);
		aaSet.setMaxNumberOfVariableModificationsPerPeptide(numMods);
		return aaSet;
	}
	
	public List<Modification.Instance> getModifications()
	{
		return 	modifications;
	}
	
	/**
	 * Gets standard amino acids from file
	 * @param fileName amino acid set file name.
	 * @return amino acid set object. 
	 */
	public static AminoAcidSet getAminoAcidSet(String fileName)
	{
		AminoAcidSet aaSet = new AminoAcidSet();
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String s;
		int lineNum = 0;
		int fileType = 0;	// 0: G,Glycine,57.021464   1: G=57.021463723
		while((s=in.readLine()) != null)
		{
			lineNum++;
			if(s.startsWith("#") || s.length() == 0)
				continue;
			if(fileType == 0 && Character.isDigit(s.charAt(0)))
			{
				fileType = 1;
				continue;
			}

			AminoAcid aa;
			if(fileType == 0)	// composition is available e.g. G, Glycine, C2H3N1O1
			{
				String[] token = s.split(",");
				if(token.length != 3)
					continue;
				String residueStr = token[0].trim();
				if(residueStr.length() != 1)
				{
					System.err.println("Illegal amino acid file format: " + fileName);
					System.err.println("Residue must be a single character: " + s);
					System.exit(-1);
				}
				char residue = residueStr.charAt(0);
				if(!Character.isUpperCase(residue))
				{
					System.err.println("Illegal amino acid file format: " + fileName);
					System.err.println("Residue must be an upper case letter: " + s);
					System.exit(-1);
				}
				String name = token[1].trim();

				// composition is available
				if(token[2].matches("(C\\d+)*(H\\d+)*(N\\d+)*(O\\d+)*(S\\d+)*"))
				{
					String compositionStr = token[2].trim();	// e.g. C5H9N1O1S1
					Composition composition = new Composition(compositionStr);
					aa = AminoAcid.getAminoAcid(residue, name, composition);
				}
				else	// composition is not available, there should be a mass
				{
					double mass = -1;
					try {
						mass = Double.parseDouble(token[2]);
					} catch (NumberFormatException e)
					{
						System.err.println("Illigal AASet File format at line "+ lineNum + ": " + s);
						System.exit(-1);
					} 
					aa = AminoAcid.getCustomAminoAcid(residue, name, mass);
				}
			}
			else	// fileType == 1, only masses (and probabilities) are available (e.g. D=115 or D=115,0.0467)
			{
				String[] token = s.split("=");
				if(token.length != 2 || token[0].length() != 1 || !Character.isLetter(token[0].charAt(0)))
				{
					System.err.println("Illegal AASet File format at line" + lineNum + ": " + s);
					System.exit(-1);
				}
				char residue = token[0].charAt(0);
				String name = token[0];
				float mass = -1;
				float prob = 0.05f;
				try {
					if(!token[1].contains(","))
						mass = Float.parseFloat(token[1]);
					else
					{
						mass = Float.parseFloat(token[1].split(",")[0]);
						prob = Float.parseFloat(token[1].split(",")[1]);
					}
				} catch (NumberFormatException e)
				{
					System.err.println("Illegal AASet File format at line" + lineNum + ": " + s);
					System.exit(-1);
				}
				if(mass <= 0)
				{
					System.err.println("Illegal AASet File format at line" + lineNum + ": " + s);
					System.exit(-1);
				}
				aa = AminoAcid.getCustomAminoAcid(residue, name, mass).setProbability(prob);
			}
			aaSet.addAminoAcid(aa);
		}
		aaSet.finalizeSet();
		try {
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return aaSet;
	}

	public static AminoAcidSet getStandardAminoAcidSet()	
	{
		if(standardAASet == null)
		{
			standardAASet = new AminoAcidSet();
			for(AminoAcid aa : AminoAcid.getStandardAminoAcids())
				standardAASet.addAminoAcid(aa);
			standardAASet.finalizeSet();
		}
		return standardAASet;
	}

	public static AminoAcidSet getStandardAminoAcidSetWithFixedCarbamidomethylatedCys()
	{
		if(standardAASetWithCarbamidomethylatedCys == null)
		{
			ArrayList<Modification.Instance> mods = new ArrayList<Modification.Instance>();
			mods.add(new Modification.Instance(Modification.get("Carbamidomethylation"), 'C').fixedModification());
			standardAASetWithCarbamidomethylatedCys = AminoAcidSet.getAminoAcidSet(mods); 
		}
		return standardAASetWithCarbamidomethylatedCys;
	}

	public static AminoAcidSet getStandardAminoAcidSetWithFixedCarboxymethylatedCys()
	{
		if(standardAASetWithCarboxyomethylatedCys == null)
		{
			ArrayList<Modification.Instance> mods = new ArrayList<Modification.Instance>();
			mods.add(new Modification.Instance(Modification.get("Carboxymethylation"), 'C').fixedModification());
			standardAASetWithCarboxyomethylatedCys = AminoAcidSet.getAminoAcidSet(mods);
		}
		return standardAASetWithCarboxyomethylatedCys;
	}
	
	/**
	 * Creates an alternative amino acid set with the terminal amino acid also 
	 * encoded.
	 * @return the AminoAcidSet with C+57 and X with an arbitrary mass.
	 */
	public static AminoAcidSet getStandardAminoAcidSetWithFixedCarbamidomethylatedCysWithTerm() {
		if(standardAASetWithCarbamidomethylatedCysWithTerm == null) {
			Modification.Instance[] mods = { 
					new Modification.Instance(Modification.get("Carbamidomethylation"), 'C').fixedModification()
			};

			HashMap<Character,Modification.Instance> modTable = new HashMap<Character,Modification.Instance>();
			for(Modification.Instance mod : mods)
			{
				if(mod.isFixedModification()) // variable modifications will be ignored
					modTable.put(mod.getResidue(), mod);
			}
			AminoAcidSet aaSet = new AminoAcidSet();
			for(AminoAcid aa : AminoAcid.getStandardAminoAcids())
			{
				Modification.Instance mod = modTable.get(aa);
				if(mod == null)
					aaSet.addAminoAcid(aa);
				else
					aaSet.addAminoAcid(aa.getAAWithFixedModification(mod.getModification()));
			}
			// terminal has 60 has mass, this is arbitrary
//			aaSet.registerAminoAcid(new AminoAcid('X', "STOP", new Composition(2,6,1,1,0)));
			// modified by Sangtae
			aaSet.addAminoAcid(AminoAcid.getCustomAminoAcid('X', new Composition(2,6,1,1,0).getMass()));

			standardAASetWithCarbamidomethylatedCysWithTerm = aaSet.finalizeSet();
		}
		return standardAASetWithCarbamidomethylatedCysWithTerm;
	}	
	
	public static AminoAcidSet getAminoAcidSet(ArrayList<Modification.Instance> mods)
	{
		AminoAcidSet aaSet = new AminoAcidSet();
		for(AminoAcid aa : getStandardAminoAcidSet())
			aaSet.addAminoAcid(aa);
		
		aaSet.applyModifications(mods);
		aaSet.finalizeSet();
		
		return aaSet;
	}

	public static AminoAcidSet getAminoAcidSet(AminoAcidSet baseAASet, ArrayList<Modification.Instance> mods)
	{
		AminoAcidSet aaSet = new AminoAcidSet();
		for(AminoAcid aa : baseAASet)
			aaSet.addAminoAcid(aa);
		
		aaSet.applyModifications(mods);
		aaSet.finalizeSet();
		
		return aaSet;
	}

	public static AminoAcidSet getAminoAcidSetFromModAAList(AminoAcidSet baseAASet, ArrayList<AminoAcid> modAAList)
	{
		AminoAcidSet aaSet = new AminoAcidSet();
		for(AminoAcid aa : baseAASet)
			aaSet.addAminoAcid(aa);

		for(AminoAcid aa : modAAList)
			aaSet.addAminoAcid(aa);
		
		aaSet.finalizeSet();
		
		return aaSet;
	}
	
	// returns a new residue for modified amino acid
	private char getModifiedResidue(char unmodResidue)
	{
		if(!Character.isUpperCase(unmodResidue))
		{
			System.err.println("Invalid unmodified residue: " + unmodResidue);
			System.exit(-1);
		}
		// if lower case letter is available
		char lowerCaseR = Character.toLowerCase(unmodResidue);
		if(!modResidueSet.contains(lowerCaseR))
		{
			modResidueSet.add(lowerCaseR);
			return lowerCaseR;
		}
		
		// if not, use char value >= 128
		char symbol = this.nextResidue;
		nextResidue++;
		if(nextResidue > Character.MAX_VALUE)
		{
			System.err.println("Too many modifications!");
			System.exit(-1);
		}
		return symbol;
	}

	public static void main(String argv[])
	{
//		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile("/home/sangtaekim/Test/Matt/MSGFDB_Mods.txt");
//		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile("/home/sangtaekim/Research/Data/IPRG2012/params.xml");
//		DBScanner.setAminoAcidProbabilities("/home/sangtaekim/Research/Data/CommonContaminants/IPI_human_3.79_withContam.fasta", aaSet);
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile("/home/sangtaekim/Test/June/params.xml");
		
		aaSet.printAASet();
//		for(AminoAcid aa : aaSet.getAminoAcids(Location.N_Term, 'E'))
//			System.out.println(aa.getResidueStr()+"\t"+aa.getMass());
	}
}
