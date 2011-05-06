package msutil;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

import parser.BufferedLineReader;

import msutil.Modification.Location;

/**
 * A factory class to instantiate a set of amino acids
 * @author sangtaekim
 *
 */

public class AminoAcidSetWithTerm implements Iterable<AminoAcid> {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static final AminoAcid[] EMPTY_AA_ARRAY = new AminoAcid[0];

	ArrayList<AminoAcid> aaList;	// contains all amino acids
	
	// for fast indexing
	private HashMap<Character,AminoAcid> residueMap;	// residue -> aa (residue must be unique)
	private HashMap<AminoAcid,Integer> indexMap;	// aa -> index (index is unique)
	private HashMap<Integer,AminoAcid[]> nominalMass2aa;	// nominalMass -> array of amino acids
	private int maxMass, minMass;

	// partitioned into positions
	private HashMap<Location, ArrayList<AminoAcid>> aaListMap;
	private HashMap<Location, HashMap<Character,AminoAcid[]>> standardResidueAAArrayMap; // std residue -> array of amino acids 
	
	// terminal fixed modification
	
	
	private int maxNumberOfVariableModificationsPerPeptide = 2;
	
	private HashSet<Character> modResidueSet = new HashSet<Character>();	// set of symbols used for residues
	
	private AminoAcidSetWithTerm() // prevents instantiation 
	{
		aaListMap = new HashMap<Location, ArrayList<AminoAcid>>();
		standardResidueAAArrayMap = new HashMap<Location, HashMap<Character,AminoAcid[]>>();
		for(Location location : Location.values())
		{
			aaListMap.put(location, new ArrayList<AminoAcid>());
			standardResidueAAArrayMap.put(location, new HashMap<Character,AminoAcid[]>());
		}
		aaList = null;
	}	

	/**
	 * Returns an amino acid object with a given index
	 * @param index amino acid index
	 * @return amino acid object
	 */
	public AminoAcid getAminoAcid(int index)	{ return aaList.get(index); }

	/**
	 * Retrieve an array of amino acids given the specific nominal mass.
	 * @param nominalMass the mass to look up
	 * @return the array of amino acids or an empty list otherwise
	 */
	public AminoAcid[] getAminoAcids(int nominalMass) {
		if (nominalMass2aa.containsKey(nominalMass)) return nominalMass2aa.get(nominalMass);
		return EMPTY_AA_ARRAY;
	}

	/**
	 * Returns the list of amino acids specific to the position.
	 * @return list of intermediate amino acids.
	 */
	public ArrayList<AminoAcid> getAAList(Location location)
	{
		return aaListMap.get(location);
	}
	
	/**
	 * Returns the iterator of 
	 */
	@Override
	public Iterator<AminoAcid> iterator() {
		return aaListMap.get(Location.Anywhere).iterator();
	}
	
	/**
	 * Retrieve an array of amino acids given the specific standard residue. 
	 * @param standardAAResidue the standard residue to look up
	 * @return the array of amino acids or null otherwise
	 */
	public AminoAcid[] getAminoAcids(Location location, char standardAAResidue) {
		return standardResidueAAArrayMap.get(location).get(standardAAResidue);
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
	 * Get the index of a given amino acid object
	 * @param aa amino acid
	 * @return non-negative index if aa belongs to this amino acid set, -1 otherwise.
	 */
	public int getIndex(AminoAcid aa) { Integer index = indexMap.get(aa); if(index==null) return -1; else return index; }
	
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
	 * Get the index of the corresponding residue
	 * @param residue
	 * @return
	 */
	public int getIndex(char residue) { return getIndex(getAminoAcid(residue)); }

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
			if(aa.isModified())
				isModified = true;
			aaArray.add(aa);
		}
		Peptide pep = new Peptide(aaArray);
		pep.setModified(isModified);

		return pep;
	}	

	private void addAminoAcid(AminoAcid aa)
	{
		addAminoAcid(aa, Location.Anywhere);
	}

	// private members
	private void addAminoAcid(AminoAcid aa, Location location)
	{
		aaList.add(aa);
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
		// Group modification instances into types
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

		// Apply anywhere fixed modifications
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
			AminoAcid targetAA = getAminoAcid(residue);
			if(targetAA == null)
			{
				System.err.println("Invalid modification: " + mod);
				System.exit(-1);
			}
			char modResidue = this.getModifiedResidue(residue);
			ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod.getModification(), modResidue);
			this.addAminoAcid(modAA, location);
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
		for(Modification.Instance mod : fixedMods.get(location))
		{
			char residue = mod.getResidue();
			if(residue != '*')	
			{
				// residue-specific fixed terminal mods are not supported
				System.err.println("Invalid modification: " + mod);
				System.exit(-1);
			}
			ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
			for(AminoAcid targetAA : this.aaListMap.get(location))
			{
				char modResidue = this.getModifiedResidue(residue);
				ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod.getModification(), modResidue);
				aaList.add(modAA);
				newAAList.add(modAA);
			}
			aaListMap.put(location, newAAList);
		}
	}
	
	private AminoAcidSet finalizeSet()
	{
		Collections.sort(aaList);
		int index=-1;

		// initialize indexMap and residueMap
		residueMap = new HashMap<Character, AminoAcid>();
		indexMap = new HashMap<AminoAcid, Integer>();
		
		for(AminoAcid aa : aaList)
		{
			index++;
			indexMap.put(aa, index);
			assert(residueMap.get(aa.getResidue()) == null): aa.getResidue()+" already exists!";
			residueMap.put(aa.getResidue(), aa);
		}			

		// initialize max and min masses
		this.minMass = Integer.MAX_VALUE;
		this.maxMass = Integer.MIN_VALUE;

		// initialize the mass2aaa array and standardResidueMap
		this.nominalMass2aa = new HashMap<Integer,AminoAcid[]>();
		standardResidueAAArrayMap = new HashMap<Character,AminoAcid[]>();
		standardResidueAAArrayMapNTerm = new HashMap<Character,AminoAcid[]>();
		
		HashMap<Integer,ArrayList<AminoAcid>> mass2aaArray = new HashMap<Integer,ArrayList<AminoAcid>>();
		HashMap<Character,LinkedList<AminoAcid>> stdResidue2aaArray = new HashMap<Character,LinkedList<AminoAcid>>();
		HashMap<Character,LinkedList<AminoAcid>> stdResidue2NTermaaArray = new HashMap<Character,LinkedList<AminoAcid>>();
		
		for (AminoAcid aa : aaList) 
		{
			int thisMass = aa.getNominalMass();
			if (!mass2aaArray.containsKey(thisMass)) {
				mass2aaArray.put(thisMass, new ArrayList<AminoAcid>());
			}
			mass2aaArray.get(thisMass).add(aa);
			
			char stdResidue;
			if(!aa.isModified())
				stdResidue = aa.getResidue();
			else
				stdResidue = ((ModifiedAminoAcid)aa).getUnmodResidue();
			
			LinkedList<AminoAcid> aaList = stdResidue2aaArray.get(stdResidue);
			LinkedList<AminoAcid> aaListNTerm = stdResidue2NTermaaArray.get(stdResidue);
			if(aaList == null)
			{
				aaList = new LinkedList<AminoAcid>();
				aaListNTerm = new LinkedList<AminoAcid>();
			}
			if(!aa.isModified())
			{
				aaList.addFirst(aa);	// unmodified residue is at first
				aaListNTerm.addFirst(aa);
			}
			else
			{
				aaList.addLast(aa);
				aaListNTerm.addLast(aa);
			}
			stdResidue2aaArray.put(stdResidue, aaList);
			stdResidue2NTermaaArray.put(stdResidue, aaListNTerm);
					
			if (thisMass > this.maxMass) this.maxMass = thisMass;
			if (thisMass < this.minMass) this.minMass = thisMass;
		}
		
		
		// N-term
		for(ModifiedAminoAcid aa : this.nTermVariableModAA)
		{
			char stdResidue = aa.getUnmodResidue();
			if(!Character.isUpperCase(stdResidue))	
				continue;	// ignore non residue-specific aa
			LinkedList<AminoAcid> aaList = stdResidue2NTermaaArray.get(stdResidue);
			aaList.addLast(aa);
			stdResidue2NTermaaArray.put(stdResidue, aaList);
		}
		
		// convert the array back to real arrays
		for (int mass : mass2aaArray.keySet()) {
			nominalMass2aa.put(mass, mass2aaArray.get(mass).toArray(new AminoAcid[0]));
		}
		
		for(char residue : stdResidue2aaArray.keySet()) 
			this.standardResidueAAArrayMap.put(residue, stdResidue2aaArray.get(residue).toArray(new AminoAcid[0]));

		for(char residue : stdResidue2NTermaaArray.keySet()) 
			this.standardResidueAAArrayMapNTerm.put(residue, stdResidue2NTermaaArray.get(residue).toArray(new AminoAcid[0]));
		
		return this;
	}	

	// static members
	private static AminoAcidSet standardAASet = null;
	private static AminoAcidSet standardAASetWithCarbamidomethylatedCys = null;
	private static AminoAcidSet standardAASetWithCarboxyomethylatedCys = null;
	private static AminoAcidSet standardAASetWithCarbamidomethylatedCysWithTerm = null;

	public int getMaxMass() { return this.maxMass; }
	public int getMinMass() { return this.minMass; }
	
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
		// if lower case is available
		char lowerCaseR = Character.toLowerCase(unmodResidue);
		if(!modResidueSet.contains(lowerCaseR))
		{
			modResidueSet.add(lowerCaseR);
			return lowerCaseR;
		}
		
		String nonStandardAAChars = "bjouxz@#$%&";
		for(int i=0; i<nonStandardAAChars.length(); i++)
		{
			char c = nonStandardAAChars.charAt(i);
			if(!modResidueSet.contains(c))
			{
				modResidueSet.add(c);
				return c;
			}
		}

		char symbol = (char)128;
		while(modResidueSet.contains(symbol))
			symbol++;
		
		if(symbol >= (1 << Character.SIZE))
		{
			System.err.println("Too many modifications!");
			System.exit(-1);
		}
		modResidueSet.add(symbol);
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
