package msutil;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;

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
	
	// for fast indexing
	private HashMap<Character,AminoAcid> residueMap;	// residue -> aa (residue must be unique)
	private HashMap<AminoAcid, Integer> aa2index;		// aa -> index
	private HashMap<Location, HashMap<Character,AminoAcid[]>> standardResidueAAArrayMap; // std residue -> array of amino acids 
	private HashMap<Location,HashMap<Integer,AminoAcid[]>> nominalMass2aa;	// nominalMass -> array of amino acids
	
	private AminoAcid[] allAminoAcidArr;
	private int maxNumberOfVariableModificationsPerPeptide = 2;
	
	private HashSet<Character> modResidueSet = new HashSet<Character>();	// set of symbols used for residues
	private char nextResidue;
	
	AminoAcid lightestAA, heaviestAA;
	
	private AminoAcidSet() // prevents instantiation 
	{
		aaListMap = new HashMap<Location, ArrayList<AminoAcid>>();
		standardResidueAAArrayMap = new HashMap<Location, HashMap<Character,AminoAcid[]>>();
		for(Location location : Location.values())
			aaListMap.put(location, new ArrayList<AminoAcid>());
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
	@Override
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
	
//	/**
//	 * Returns the largest char residue code plus 1 (if the value is smaller than 128, return 128)
//	 * @return largest char residue code + 1
//	 */
//	public int getMaxResidue()
//	{
//		return maxResidue;
//	}
	
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
			if(aa.isVariableModification())
				isModified = true;
			aaArray.add(aa);
		}
		Peptide pep = new Peptide(aaArray);
		pep.setModified(isModified);

		return pep;
	}	

	public int getMaxNominalMass() { return this.lightestAA.getNominalMass(); }
	public int getMinNominalMass() { return this.heaviestAA.getNominalMass(); }
	
	public AminoAcid getLightestAA()	{ return this.lightestAA; }
	public AminoAcid getHeaviestAA()	{ return this.heaviestAA; }
	
	// private members to build an amino acid set
	private void addAminoAcid(AminoAcid aa)
	{
		addAminoAcid(aa, Location.Anywhere);
	}

	// private members
	private void addAminoAcid(AminoAcid aa, Location location)
	{
		if(location == Location.Anywhere)
		{
			for(Location loc : Location.values())
				this.aaListMap.get(loc).add(aa);
		}
		else if(location == Location.N_Term)
		{
			this.aaListMap.get(Location.N_Term).add(aa);
			this.aaListMap.get(Location.Protein_N_Term).add(aa);
		}
		else if(location == Location.C_Term)
		{
			this.aaListMap.get(Location.C_Term).add(aa);
			this.aaListMap.get(Location.Protein_C_Term).add(aa);
		}
		else if(location == Location.Protein_N_Term)
		{
			this.aaListMap.get(Location.Protein_N_Term).add(aa);
		}
		else if(location == Location.Protein_C_Term)
		{
			this.aaListMap.get(Location.Protein_C_Term).add(aa);
		}
	}
	
	private void applyModifications(ArrayList<Modification.Instance> mods)
	{
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

		// Apply anywhere fixed modifications: simply change the amino acid
		for(Modification.Instance mod : fixedMods.get(Location.Anywhere))
		{
			char residue = mod.getResidue();
			AminoAcid aa = getAminoAcid(residue);
			if(aa == null)
			{
				System.err.println("Invalid modification: " + mod);
				System.exit(-1);
			}
			aa.applyFixedModification(mod.getModification());
		}
		
		// Apply anywhere variable modifications
		addVariableMods(variableMods, Location.Anywhere);
		
		// Variable N-term mods (residue-specific)
		addVariableMods(variableMods, Location.N_Term);
		
		// Variable C-term mods (residue-specific)
		addVariableMods(variableMods, Location.C_Term);
		
		// Fixed N-term
		addFixedTerminalMods(fixedMods, Location.N_Term);
		
		// Fixed C-term
		addFixedTerminalMods(fixedMods, Location.C_Term);

		// Fixed N-term
		addFixedTerminalMods(fixedMods, Location.Protein_N_Term);
		
		// Fixed C-term
		addFixedTerminalMods(fixedMods, Location.Protein_C_Term);
		
	}

	private void addVariableMods(HashMap<Modification.Location,ArrayList<Modification.Instance>> variableMods, Location location)
	{
		// residue-specific
		for(Modification.Instance mod : variableMods.get(location))
		{
			char residue = mod.getResidue();
			if(residue == '*')
				continue;
			for(AminoAcid targetAA : getAminoAcids(location, residue))
			{
				char modResidue = this.getModifiedResidue(residue);
				ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod.getModification(), modResidue);
				this.addAminoAcid(modAA, location);
			}
		}
		
		// any residue
		for(Modification.Instance mod : variableMods.get(location))
		{
			char residue = mod.getResidue();
			if(residue != '*')
				continue;
			for(AminoAcid targetAA : this.aaListMap.get(location))
			{
				char modResidue = this.getModifiedResidue(residue);
				ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod.getModification(), modResidue);
				this.addAminoAcid(modAA, location);
			}
		}
	}
	
	private void addFixedTerminalMods(HashMap<Modification.Location,ArrayList<Modification.Instance>> fixedMods, Location location)
	{
		// residue-specific
		for(Modification.Instance mod : fixedMods.get(location))
		{
			char residue = mod.getResidue();
			if(residue == '*')
				continue;
			ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
			for(AminoAcid targetAA : getAminoAcids(location, residue))
			{
				char modResidue = this.getModifiedResidue(residue);
				ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod.getModification(), modResidue);
				newAAList.add(modAA);
			}
			aaListMap.put(location, newAAList);
		}
		
		

		// any residue
		for(Modification.Instance mod : fixedMods.get(location))
		{
			char residue = mod.getResidue();
			if(residue != '*')
				continue;
			ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
			for(AminoAcid targetAA : this.aaListMap.get(location))
			{
				char modResidue = this.getModifiedResidue(residue);
				ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod.getModification(), modResidue);
				newAAList.add(modAA);
			}
			aaListMap.put(location, newAAList);
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
		ArrayList<AminoAcid> allAAList = new ArrayList<AminoAcid>();
		for(Location location : aaListMap.keySet())
		{
			for(AminoAcid aa : aaListMap.get(location))
				allAAList.add(aa);
		}
		
		Collections.sort(allAAList);
		this.allAminoAcidArr = allAAList.toArray(EMPTY_AA_ARRAY);

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
		
		for(AminoAcid aa : allAAList)
		{
			assert(residueMap.get(aa.getResidue()) == null): aa.getResidue()+" already exists!";
			residueMap.put(aa.getResidue(), aa);
		}			

		for(Location location : Location.values())
		{
			HashMap<Integer,ArrayList<AminoAcid>> mass2aaList = new HashMap<Integer,ArrayList<AminoAcid>>();
			HashMap<Character,LinkedList<AminoAcid>> stdResidue2aaList = new HashMap<Character,LinkedList<AminoAcid>>();
			
			for (AminoAcid aa : this.aaListMap.get(location)) 
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
				if(!aa.isVariableModification())
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
			if(s.startsWith("#") || s.trim().length() == 0)
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
		aaSet.finalizeSet();
		
		aaSet.applyModifications(mods);
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
//		testReadingFile();
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile("/home/sangtaekim/Research/ToolDistribution/Mods.txt");
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		System.out.println("Standard aa and Anyware variable modifications");
//		for(AminoAcid aa : aaSet)
//			System.out.println(aaSet.getIndex(aa)+"\t"+aa.getResidueStr()+"\t"+aa.getResidue()+"\t"+aa.getMass()+"\t"+aa.getName()+"\t"+aa.isModified());
//		System.out.println("Terminal aa");
//		for(AminoAcid aa : aaSet.getAminoAcids('E'))
//			System.out.println(aa.getResidue());
//		for(AminoAcid aa : aaSet.getNTermAminoAcids('E'))
//			System.out.println(aa.getResidue());
	}
}
