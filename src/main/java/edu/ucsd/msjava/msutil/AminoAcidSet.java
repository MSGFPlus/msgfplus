package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.msutil.Modification.Location;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.ui.MSGFPlus;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

/**
 * A factory class to instantiate a set of amino acids
 *
 * @author sangtaekim
 */
public class AminoAcidSet implements Iterable<AminoAcid> {
    /**
     *
     */
    private static final AminoAcid[] EMPTY_AA_ARRAY = new AminoAcid[0];

    private HashMap<Location, ArrayList<AminoAcid>> aaListMap;

    private static HashMap<Location, Location[]> locMap;

    /**
     * This tracks any default mods that the user has defined
     * Keys are mod names and values are the mod mass that the user defined for this modification
     * This list is used to warn users of non-standard mod masses for default mods
     */
    private static Hashtable<String, Double> defaultModUsage = new Hashtable<>();

    static {
        locMap = new HashMap<>();
        locMap.put(Location.Anywhere, new Location[]{Location.Anywhere, Location.N_Term, Location.C_Term, Location.Protein_N_Term, Location.Protein_C_Term});
        locMap.put(Location.N_Term, new Location[]{Location.N_Term, Location.Protein_N_Term});
        locMap.put(Location.C_Term, new Location[]{Location.C_Term, Location.Protein_C_Term});
        locMap.put(Location.Protein_N_Term, new Location[]{Location.Protein_N_Term});
        locMap.put(Location.Protein_C_Term, new Location[]{Location.Protein_C_Term});
    }

    // for fast indexing
    private HashMap<Character, AminoAcid> residueMap;    // residue -> aa (residue must be unique)
    private HashMap<AminoAcid, Integer> aa2index;        // aa -> index
    private HashMap<Location, HashMap<Character, AminoAcid[]>> standardResidueAAArrayMap; // std residue -> array of amino acids
    private HashMap<Location, HashMap<Integer, AminoAcid[]>> nominalMass2aa;    // nominalMass -> array of amino acids

    private AminoAcid[] allAminoAcidArr;
    private int maxNumberOfVariableModificationsPerPeptide = 2;

    private boolean containsModification;    // true if this contains any variable or terminal (fixed or variable) modification
    private boolean containsNTermModification;    // true if this contains any (fixed or variable) modification specific to N-terminus
    private boolean containsCTermModification;    // true if this contains any (fixed or variable) modification specific to N-terminus
    private boolean containsPhosphorylation;    // true if this contains phosphorylation
    private boolean containsITRAQ;    // true if this contains iTRAQ
    private boolean containsTMT;    // true if this contains iTRAQ

    private HashSet<Character> modResidueSet = new HashSet<>();    // set of symbols used for residues
    private char nextResidue;

    // for enzyme
//	private ArrayList<AminoAcid> enzymeAAList;
    private int neighboringAACleavageCredit = 0;
    private int neighboringAACleavagePenalty = 0;
    private int peptideCleavageCredit = 0;
    private int peptideCleavagePenalty = 0;
    private float probCleavageSites = 0;

    AminoAcid lightestAA, heaviestAA;

    /**
     * This tracks user-friendly descriptions of the modifications in use
     */
    private ArrayList<String> modificationsInUse = new ArrayList<>();

    private AminoAcidSet() // prevents instantiation
    {
        aaListMap = new HashMap<Location, ArrayList<AminoAcid>>();
        standardResidueAAArrayMap = new HashMap<Location, HashMap<Character, AminoAcid[]>>();
        for (Location location : Location.values()) {
            aaListMap.put(location, new ArrayList<AminoAcid>());
        }
        nextResidue = 128;
    }

    /**
     * Returns the list of amino acids specific to the position.
     *
     * @return list of intermediate amino acids.
     */
    public ArrayList<AminoAcid> getAAList(Location location) {
        return aaListMap.get(location);
    }

    public ArrayList<AminoAcid> getNTermAAList() {
        return aaListMap.get(Location.N_Term);
    }

    public ArrayList<AminoAcid> getCTermAAList() {
        return aaListMap.get(Location.C_Term);
    }

    public ArrayList<AminoAcid> getProtNTermAAList() {
        return aaListMap.get(Location.Protein_N_Term);
    }

    public ArrayList<AminoAcid> getProtCTermAAList() {
        return aaListMap.get(Location.Protein_N_Term);
    }

    public ArrayList<String> getModificationsInUse() {
        return modificationsInUse;
    }

    /**
     * Returns the iterator of anywhere amino acids
     */
    public Iterator<AminoAcid> iterator() {
        return aaListMap.get(Location.Anywhere).iterator();
    }

    /**
     * Returns the size of amino acid depending on the location.
     *
     * @param location amino acid location
     * @return
     */
    public int size(Location location) {
        return aaListMap.get(location).size();
    }

    /**
     * Reterns the size of anywhere amino acids
     *
     * @return the size of anywhere amino acids
     */
    public int size() {
        return aaListMap.get(Location.Anywhere).size();
    }

    /**
     * Retrieve an array of amino acids given the specific standard residue.
     *
     * @param location          amino acid location
     * @param standardAAResidue the standard residue to look up
     * @return the array of amino acids or an empty array otherwise
     */
    public AminoAcid[] getAminoAcids(Location location, char standardAAResidue) {
        AminoAcid[] matches = standardResidueAAArrayMap.get(location).get(standardAAResidue);
        if (matches != null)
            return matches;
        else
            return EMPTY_AA_ARRAY;
    }

    /**
     * Retrieve an array of amino acids given the specific nominal mass.
     *
     * @param location    amino acid location
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
     *
     * @param nominalMass the mass to look up
     * @return the array of amino acids or an empty list otherwise
     */
    public AminoAcid[] getAminoAcids(int nominalMass) {
        return getAminoAcids(Location.Anywhere, nominalMass);
    }

    /**
     * Checks whether a residue belongs to this amino acid set
     *
     * @param residue a residue
     * @return true if residue belongs to the amino acid set
     */
    public boolean contains(char residue) {
        return residueMap.containsKey(residue);
    }

    /**
     * Returns a list of all residues without mods
     *
     * @return
     */
    public ArrayList<Character> getResidueListWithoutMods() {
        ArrayList<Character> residues = new ArrayList<Character>();
        for (Map.Entry<Character, AminoAcid> aa : residueMap.entrySet()) {
            char residue = aa.getValue().getUnmodResidue();
            if (!residues.contains(residue)) {
                residues.add(residue);
            }
        }
        return residues;
    }

    /**
     * Returns a list of all residues, including modified residues
     *
     * @return
     */
    public ArrayList<Character> getResidueList() {
        return new ArrayList<Character>(residueMap.keySet());
    }

    /**
     * Get the amino acid mass of the residue.
     *
     * @param residue the amino acid mass. Use uppercase for standard aa (convention).
     *                this method is case sensitive.
     * @return the amino acid object. null if no aa corresponding to the residue
     */
    public AminoAcid getAminoAcid(Location location, char residue) {
        AminoAcid[] aaArr = getAminoAcids(location, residue);
        for (AminoAcid aa : aaArr)
            if (!aa.isModified())
                return aa;
        return null;
    }

    /**
     * Get the amino acid mass of the residue.
     *
     * @param residue the amino acid mass. Use uppercase for standard aa (convention).
     *                this method is case sensitive.
     * @return the amino acid object. null if no aa corresponding to the residue
     */
    public AminoAcid getAminoAcid(char residue) {
        return residueMap.get(residue);
    }

    /**
     * Set the number of allowable variable modifications per peptide
     *
     * @param maxNumberOfVariableModificationsPerPeptide the number of allowable variable modifications per peptide
     */
    public void setMaxNumberOfVariableModificationsPerPeptide(int maxNumberOfVariableModificationsPerPeptide) {
        this.maxNumberOfVariableModificationsPerPeptide = maxNumberOfVariableModificationsPerPeptide;
    }

    /**
     * Get the number of allowable variable modifications per peptide
     *
     * @return the number of allowable variable modifications per peptide
     */
    public int getMaxNumberOfVariableModificationsPerPeptide() {
        return this.maxNumberOfVariableModificationsPerPeptide;
    }

    /**
     * Get all amino acids for all locations.
     *
     * @return an array of all amino acids.
     */
    public AminoAcid[] getAllAminoAcidArr() {
        return this.allAminoAcidArr;
    }

    /**
     * Get the amino acid corresponding to the index
     *
     * @param index amino acid index
     * @return amino acid object
     */
    public AminoAcid getAminoAcid(int index) {
        return allAminoAcidArr[index];
    }

    /**
     * Get the index of the aa
     *
     * @param aa amino acid
     * @return the index of aa. null if aa does not belong to this amino acid set
     */
    public int getIndex(AminoAcid aa) {
        Integer index = aa2index.get(aa);
        if (index == null)
            index = -1;
        return index;
    }

    /**
     * Get the peptide corresponding to the string sequence.
     *
     * @param sequence sequence of the peptide.
     * @return peptide object of the sequence
     */
    public Peptide getPeptide(String sequence) {
        boolean isModified = false;
        ArrayList<AminoAcid> aaArray = new ArrayList<AminoAcid>();
        for (int i = 0; i < sequence.length(); i++) {
            char residue = sequence.charAt(i);
            AminoAcid aa = this.getAminoAcid(residue);
            if (aa == null) {
                System.out.println(sequence + ": " + residue + " is null!");
            }
            assert (aa != null) : sequence + ": " + residue + " is null!";
            if (aa.isModified())
                isModified = true;
            aaArray.add(aa);
        }
        Peptide pep = new Peptide(aaArray);
        pep.setModified(isModified);

        return pep;
    }

    public int getMaxNominalMass() {
        return this.heaviestAA.getNominalMass();
    }

    public int getMinNominalMass() {
        return this.lightestAA.getNominalMass();
    }

    public AminoAcid getLightestAA() {
        return this.lightestAA;
    }

    public AminoAcid getHeaviestAA() {
        return this.heaviestAA;
    }

    public boolean containsModification() {
        return this.containsModification;
    }

    public boolean containsNTermModification() {
        return this.containsNTermModification;
    }

    public boolean containsCTermModification() {
        return this.containsCTermModification;
    }

    public boolean containsPhosphorylation() {
        return this.containsPhosphorylation;
    }

    public boolean containsITRAQ() {
        return this.containsITRAQ;
    }

    public boolean containsTMT() {
        return this.containsTMT;
    }

    public char getMaxResidue() {
        return nextResidue;
    }

    public void registerEnzyme(Enzyme enzyme) {
        if (enzyme == null || enzyme.getResidues() == null ||
                enzyme.getPeptideCleavageEfficiency() == 0 || enzyme.getNeighboringAACleavageEfficiency() == 0)
            return;

        probCleavageSites = 0;
        for (char residue : enzyme.getResidues()) {
            AminoAcid aa = this.getAminoAcid(residue);
            if (aa == null) {
                System.err.println("Invalid Enzyme cleavage site: " + residue);
                System.exit(-1);
            }
            probCleavageSites += aa.getProbability();
        }

        if (probCleavageSites == 0 || probCleavageSites == 1) {
            System.err.println("Probability of enzyme residues must be in (0,1)!");
            System.exit(-1);
        }

        float peptideCleavageEfficiency = enzyme.getPeptideCleavageEfficiency();
        float neighboringAACleavageEfficiency = enzyme.getNeighboringAACleavageEfficiency();

        peptideCleavageCredit = (int) Math.round(Math.log(peptideCleavageEfficiency / probCleavageSites));
        peptideCleavagePenalty = (int) Math.round(Math.log((1 - peptideCleavageEfficiency) / (1 - probCleavageSites)));
        neighboringAACleavageCredit = (int) Math.round(Math.log(neighboringAACleavageEfficiency / probCleavageSites));
        neighboringAACleavagePenalty = (int) Math.round(Math.log((1 - neighboringAACleavageEfficiency) / (1 - probCleavageSites)));
    }

    public int getNeighboringAACleavageCredit() {
        return neighboringAACleavageCredit;
    }

    public int getNeighboringAACleavagePenalty() {
        return neighboringAACleavagePenalty;
    }

    public int getPeptideCleavageCredit() {
        return peptideCleavageCredit;
    }

    public int getPeptideCleavagePenalty() {
        return peptideCleavagePenalty;
    }

    public float getProbCleavageSites() {
        return probCleavageSites;
    }

    public void printAASet() {
        System.out.println("NumMods: " + this.getMaxNumberOfVariableModificationsPerPeptide());
        for (Location location : Location.values()) {
            ArrayList<AminoAcid> aaList = this.getAAList(location);
            System.out.println(location + "\t" + aaList.size());
            for (AminoAcid aa : aaList)
                System.out.println(aa.getResidueStr() + (aa.isModified() ? "*" : "") + "\t" + (int) aa.getResidue() + "\t" + aa.getNominalMass() + "\t" + aa.getMass() + "\t" + aa.getProbability());
        }
    }

    // private members to build an amino acid set
    private void addAminoAcid(AminoAcid aa) {
        addAminoAcid(aa, Location.Anywhere);
    }

    // private members
    private void addAminoAcid(AminoAcid aa, Location location) {
        for (Location loc : locMap.get(location)) {
            // Debug
//			if(aa.isModified())
//				System.out.println("Debug");
            aaListMap.get(loc).add(aa);
        }
    }

    private List<Modification.Instance> modifications;

    /**
     * Add a dynamic or static modification that applies to a residue or the N- or C-terminus
     * @param modFileName Mod file name
     * @param lineNum Line number
     * @param dataLine Text from this line in the mod file
     * @param mods Existing mod instances
     * @param modIns New mod instance
     * @return True if successful, false if the same modification is defined for the same residue twice
     */
    private static boolean addModInstance(
            String modFileName, int lineNum, String dataLine,
            ArrayList<Modification.Instance> mods, Modification.Instance modIns) {

        for (Modification.Instance comparisonItem : mods) {
            if (modIns.getResidue() == comparisonItem.getResidue() &&
                    modIns.getLocation() == comparisonItem.getLocation() &&
                    modIns.getModification().getName().equals(comparisonItem.getModification().getName())) {
                System.err.println(
                        "Error: The same modification is defined for the same residue twice; \n" +
                        "the duplicate definition is on line " + lineNum +
                        " in file " + modFileName + ": " + dataLine);

                return false;
            }
        }

        mods.add(modIns);
        return true;
    }

    private void applyModifications(ArrayList<Modification.Instance> mods) {
        this.modifications = mods;

        modificationsInUse.clear();

        if (mods.size() == 0) {
            return;
        }

        // partition modification instances into different types
        HashMap<Modification.Location, ArrayList<Modification.Instance>> fixedMods = new HashMap<>();
        HashMap<Modification.Location, ArrayList<Modification.Instance>> variableMods = new HashMap<>();

        for (Location location : Modification.Location.values()) {
            fixedMods.put(location, new ArrayList<>());
            variableMods.put(location, new ArrayList<>());
        }

        for (Modification.Instance mod : mods) {
            if (mod.isFixedModification())
                fixedMods.get(mod.getLocation()).add(mod);
            else
                variableMods.get(mod.getLocation()).add(mod);
        }

        Location[] locArr = new Location[]{
                Location.Anywhere,
                Location.N_Term,
                Location.C_Term,
                Location.Protein_N_Term,
                Location.Protein_C_Term,
        };

        // Fixed modifications
        for (Location loc : locArr)
            applyFixedMods(fixedMods, loc);

        // Variable modifications
        for (Location loc : locArr)
            addVariableMods(variableMods, loc);

        // setup containsNTermModification and containsCTermModification
        for (Modification.Instance mod : mods) {
            Location location = mod.getLocation();
            if (!containsNTermModification && (location == Location.N_Term || location == Location.Protein_N_Term))
                this.containsNTermModification = true;
            if (!containsCTermModification && (location == Location.C_Term || location == Location.Protein_C_Term))
                this.containsCTermModification = true;
            if (location != Location.Anywhere || !mod.isFixedModification())
                this.containsModification = true;
            if (mod.getModification().getName().toLowerCase().startsWith("phospho"))
                this.containsPhosphorylation = true;
            if (mod.getModification().getName().toLowerCase().startsWith("itraq"))
                this.containsITRAQ = true;
            if (mod.getModification().getName().toLowerCase().startsWith("tmt"))
                this.containsTMT = true;

            String modType = mod.isFixedModification() ? "Fixed (static)" : "Variable (dynamic)";
            String modLocation;

            switch (mod.getLocation()) {
                case Anywhere:
                    modLocation = "";
                    break;
                case N_Term:
                    modLocation = " at the peptide N-terminus";
                    break;
                case C_Term:
                    modLocation = " at the peptide C-terminus";
                    break;
                case Protein_N_Term:
                    modLocation = " at the protein N-terminus";
                    break;
                case Protein_C_Term:
                    modLocation = " at the protein C-terminus";
                    break;
                default:
                    modLocation = " at ???";
                    break;
            }

            String modInfo = modType + ": " + mod.getModification().getName() + " on " + mod.getResidue() + modLocation;

            modificationsInUse.add(modInfo);
        }
    }

    private void applyFixedMods(HashMap<Modification.Location, ArrayList<Modification.Instance>> fixedMods, Location location) {
        for (Modification.Instance mod : fixedMods.get(location)) {
            // residue-specific
            char residue = mod.getResidue();
            if (residue == '*')
                continue;

            ArrayList<AminoAcid> oldAAList = this.getAAList(location);
            ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();

            for (AminoAcid aa : oldAAList) {
                if (aa.getUnmodResidue() != residue)
                    newAAList.add(aa);
                else {
                    if (location == Location.Anywhere)
                        newAAList.add(aa.getAAWithFixedModification(mod.getModification()));    // replace with a new amino acid
                    else {
//						char modResidue = this.getModifiedResidue(aa.getUnmodResidue());
//						ModifiedAminoAcid modAA = new ModifiedAminoAcid(aa, mod, modResidue);
                        ModifiedAminoAcid modAA = getModifiedAminoAcid(aa, mod);
                        newAAList.add(modAA);
                    }
                }
            }

            for (Location loc : locMap.get(location))
                aaListMap.put(loc, new ArrayList<AminoAcid>(newAAList));
        }

        // any residue
        for (Modification.Instance mod : fixedMods.get(location)) {
            char residue = mod.getResidue();
            if (residue != '*')
                continue;
            ArrayList<AminoAcid> oldAAList = this.getAAList(location);
            ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();

            for (AminoAcid aa : oldAAList) {
                if (location == Location.Anywhere)
                    newAAList.add(aa.getAAWithFixedModification(mod.getModification()));
                else {
//					char modResidue = this.getModifiedResidue(aa.getUnmodResidue());
//					ModifiedAminoAcid modAA = new ModifiedAminoAcid(aa, mod, modResidue);
                    ModifiedAminoAcid modAA = getModifiedAminoAcid(aa, mod);
                    newAAList.add(modAA);
                }
            }

            for (Location loc : locMap.get(location))
                aaListMap.put(loc, new ArrayList<AminoAcid>(newAAList));
        }
    }

    private void addVariableMods(HashMap<Modification.Location, ArrayList<Modification.Instance>> variableMods, Location location) {
        // residue-specific
        for (Location loc : locMap.get(location)) {
            ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
            ArrayList<AminoAcid> oldAAList = this.getAAList(loc);
            for (AminoAcid targetAA : oldAAList) {
                for (Modification.Instance mod : variableMods.get(location)) {
                    char residue = mod.getResidue();
                    if (residue == '*')
                        continue;
                    if (targetAA.getUnmodResidue() == residue) {
                        if (targetAA.isModified()) {
//							if(mod.getLocation() == Location.Anywhere && targetAA.hasResidueSpecificVariableMod())	// residue mod
//								continue;
//							if(mod.getLocation() != Location.Anywhere && targetAA.hasTerminalVariableMod())	// residue mod
//								continue;
                            if (targetAA.hasResidueSpecificVariableMod())
                                continue;
                        }
//						char modResidue = this.getModifiedResidue(targetAA.getUnmodResidue());
//						ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod, modResidue);
                        ModifiedAminoAcid modAA = getModifiedAminoAcid(targetAA, mod);
                        newAAList.add(modAA);
                    }
                }
            }
            for (AminoAcid newAA : newAAList)
                aaListMap.get(loc).add(newAA);
        }

        // any residue
        for (Location loc : locMap.get(location)) {
            ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
            ArrayList<AminoAcid> oldAAList = this.getAAList(loc);
            for (AminoAcid targetAA : oldAAList) {
                for (Modification.Instance mod : variableMods.get(location)) {
                    char residue = mod.getResidue();
                    if (residue != '*')
                        continue;
//					if(location == Location.Anywhere)
//					{
//						System.err.println("Invalid modification: " + mod);
//						System.exit(-1);
//					}
                    if (targetAA.isModified()) {
//						if(mod.getLocation() == Location.Anywhere && targetAA.hasResidueSpecificVariableMod())	// residue mod
//							continue;
//						if(mod.getLocation() != Location.Anywhere && targetAA.hasTerminalVariableMod())	// residue mod
//							continue;
                        if (targetAA.hasTerminalVariableMod())
                            continue;
                    }
//					char modResidue = this.getModifiedResidue(targetAA.getUnmodResidue());
//					ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod, modResidue);
                    ModifiedAminoAcid modAA = getModifiedAminoAcid(targetAA, mod);
                    newAAList.add(modAA);
                }
            }
            for (AminoAcid newAA : newAAList)
                aaListMap.get(loc).add(newAA);
        }
    }

    private AminoAcidSet finalizeSet() {
        standardResidueAAArrayMap = new HashMap<Location, HashMap<Character, AminoAcid[]>>();
        nominalMass2aa = new HashMap<Location, HashMap<Integer, AminoAcid[]>>();
        for (Location location : Location.values()) {
            standardResidueAAArrayMap.put(location, new HashMap<Character, AminoAcid[]>());
            nominalMass2aa.put(location, new HashMap<Integer, AminoAcid[]>());
        }

        // add all amino acids to aaList
        HashSet<AminoAcid> allAASet = new HashSet<AminoAcid>();
        for (Location location : aaListMap.keySet()) {
            for (AminoAcid aa : aaListMap.get(location))
                allAASet.add(aa);
        }

        this.allAminoAcidArr = allAASet.toArray(EMPTY_AA_ARRAY);
        Arrays.sort(allAminoAcidArr);

        // assign index, heaviest and lightest aa
        double minMass = Double.MAX_VALUE;
        int lightIndex = -1;
        double maxMass = Double.MIN_VALUE;
        int heavyIndex = -1;
        aa2index = new HashMap<AminoAcid, Integer>();        // aa -> index
        for (int i = 0; i < allAminoAcidArr.length; i++) {
            aa2index.put(allAminoAcidArr[i], i);
            double mass = allAminoAcidArr[i].getAccurateMass();
            if (mass < minMass) {
                lightIndex = i;
                minMass = mass;
            }
            if (mass > maxMass) {
                heavyIndex = i;
                maxMass = mass;
            }
        }
        this.heaviestAA = allAminoAcidArr[heavyIndex];
        this.lightestAA = allAminoAcidArr[lightIndex];

        // initialize aaList and residueMap
        residueMap = new HashMap<Character, AminoAcid>();

        for (AminoAcid aa : allAminoAcidArr) {
            assert (residueMap.get(aa.getResidue()) == null) : aa.getResidue() + " already exists!";
            residueMap.put(aa.getResidue(), aa);
        }

        for (Location location : Location.values()) {
            HashMap<Integer, ArrayList<AminoAcid>> mass2aaList = new HashMap<Integer, ArrayList<AminoAcid>>();
            HashMap<Character, LinkedList<AminoAcid>> stdResidue2aaList = new HashMap<Character, LinkedList<AminoAcid>>();

            for (AminoAcid aa : this.getAAList(location)) {
                int thisMass = aa.getNominalMass();
                if (!mass2aaList.containsKey(thisMass)) {
                    mass2aaList.put(thisMass, new ArrayList<AminoAcid>());
                }
                mass2aaList.get(thisMass).add(aa);

                char stdResidue = aa.getUnmodResidue();
                LinkedList<AminoAcid> aaList = stdResidue2aaList.get(stdResidue);
                if (aaList == null)
                    aaList = new LinkedList<AminoAcid>();
                if (!aa.isModified())
                    aaList.addFirst(aa);    // unmodified residue is at first
                else
                    aaList.addLast(aa);
                stdResidue2aaList.put(stdResidue, aaList);
            }

            // convert the array back to real arrays
            HashMap<Integer, AminoAcid[]> mass2aaArray = new HashMap<Integer, AminoAcid[]>();
            for (int mass : mass2aaList.keySet()) {
                mass2aaArray.put(mass, mass2aaList.get(mass).toArray(new AminoAcid[0]));
            }

            HashMap<Character, AminoAcid[]> stdResidue2aaArray = new HashMap<Character, AminoAcid[]>();
            for (char residue : stdResidue2aaList.keySet())
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

    public static AminoAcidSet getAminoAcidSetFromModFile(String modFilePath, ParamManager paramManager) {
        BufferedLineReader reader = null;

        File modFile = new File(modFilePath);

        try {
            reader = new BufferedLineReader(modFile.getPath());
        } catch (IOException e) {
            System.err.println("Error opening modification file " + modFile.getPath());
            e.printStackTrace();
            System.exit(-1);
        }

        // parse modifications
        ArrayList<Modification.Instance> mods = new ArrayList<>();
        ArrayList<AminoAcid> customAA = new ArrayList<>();
        String dataLine;
        String sourceFileName = modFile.getName();
        int lineNum = 0;
        int maxNumMods = paramManager.getMaxNumModsPerPeptide();
        ModificationMetadata modMetadata = new ModificationMetadata(maxNumMods);

        while ((dataLine = reader.readLine()) != null) {
            lineNum++;
            boolean success = parseConfigEntry(sourceFileName, lineNum, dataLine, mods, customAA, modMetadata);
            if (!success) {
                System.exit(-1);
            }
        }

        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(mods, customAA);

        if (numMods != modMetadata.getNumMods()) {
            aaSet.setMaxNumberOfVariableModificationsPerPeptide(modMetadata.getNumMods());
            paramManager.setMaxNumMods(modMetadata.getNumMods());
        }

        try {
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return aaSet;
    }

    /**
     * Associate modifications in rawMods with amino acids
     * @param modConfigFilePath
     * @param modsByLine Hashtable where keys are the line number in the MSGF+ parameter file and values are the text from the given line
     * @param numMods Max Number of Dynamic (Variable) Modifications per peptide
     * @return
     */
    public static AminoAcidSet getAminoAcidSetFromList(String modConfigFilePath, Hashtable<Integer, String> modsByLine, int numMods) {
        BufferedLineReader reader = null;

        // parse modifications
        ArrayList<Modification.Instance> mods = new ArrayList<>();
        ArrayList<AminoAcid> customAA = new ArrayList<>();
        int maxNumMods = paramManager.getMaxNumModsPerPeptide();
        ModificationMetadata modMetadata = new ModificationMetadata(maxNumMods);

        modsByLine.forEach((lineNum, dataLine) -> {
            boolean success = parseConfigEntry(modConfigFilePath, lineNum, dataLine, mods, customAA, modMetadata);
            if (!success) {
                System.exit(-1);
            }
         });

        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(mods, customAA);
        aaSet.setMaxNumberOfVariableModificationsPerPeptide(modMetadata.getNumMods());
        return aaSet;
    }

    private static boolean parseConfigEntry(
            String sourceFilePath, int lineNum, String dataLine,
            ArrayList<Modification.Instance> mods,
            ArrayList<AminoAcid> customAA,
            ModificationMetadata modMetadata) {

        String[] tokenArray = dataLine.split("#");
        if (tokenArray.length == 0)
            return true;

        String modSetting = tokenArray[0].trim();
        if (modSetting.length() == 0) {
            return true;
        }

        if (modSetting.toLowerCase().startsWith("nummods=")) {
            try {
                int numMods = Integer.parseInt(modSetting.split("=")[1].trim());
                modMetadata.setMaxNumModsPerPeptide(numMods);
            } catch (NumberFormatException e) {
                System.err.println("Error: Invalid NumMods option at line " + lineNum +
                        " in file " + sourceFilePath + ": " + modSetting);
                e.printStackTrace();
                return false;
            }
        } else {
            String[] modInfo = modSetting.split(",");
            if (modInfo.length < 5) {
                System.out.println("Ignoring line " + lineNum +
                        " in file " + sourceFilePath +
                        " since does not have 5 parts separated by commas: " + modSetting);
                return true;
            }

            // Mass or Composition
            double modMass = 0;
            String compStr = modInfo[0].trim();

            // First try to parse compStr as an empirical formula
            // Supports C, H, N, O, S, P, Br, Cl, Fe, and Se

            Double mass = Composition.getMass(compStr);
            if (mass != null)
                modMass = mass;
            else {
                try {
                    modMass = Double.parseDouble(compStr);
                } catch (NumberFormatException e) {
                    System.err.println("Error: Invalid Mass/Composition at line " + lineNum +
                            " in file " + sourceFilePath + ": " + modSetting);
                    e.printStackTrace();
                    return false;
                }
            }

            String customAAResidues = modMetadata.getCustomAAResidues();

            // Residues
            String residueStr = modInfo[1].trim();
            boolean isResidueStrLegitimate = true;
            boolean matchesCustomAA = false;
            if (!residueStr.equals("*")) {
                if (residueStr.length() > 0) {
                    for (int i = 0; i < residueStr.length(); i++) {
                        boolean matchesCustom = customAAResidues.indexOf(residueStr.charAt(i)) > -1;
                        if (matchesCustom) {
                            matchesCustomAA = true;
                        }
                        if (!matchesCustom && !AminoAcid.isStdAminoAcid(residueStr.charAt(i))) {
                            isResidueStrLegitimate = false;
                            break;
                        }
                    }
                } else
                    isResidueStrLegitimate = false;
            }

            // isFixedModification
            boolean isFixedModification = false;
            boolean isCustomAminoAcid = false;
            boolean modTypeParseFailed = false;

            if (modInfo[2].trim().equalsIgnoreCase("fix"))
                isFixedModification = true;
            else if (modInfo[2].trim().equalsIgnoreCase("opt"))
                isFixedModification = false;
            else if (modInfo[2].trim().equalsIgnoreCase("custom"))
                isCustomAminoAcid = true;
            else
                modTypeParseFailed = true;

            if ((!isResidueStrLegitimate && !isCustomAminoAcid) || (isCustomAminoAcid && matchesCustomAA)) {
                System.err.println("Error: Invalid Residue(s) at line " + lineNum +
                        " in file " + sourceFilePath + ": " + modSetting);
                return false;
            }
            if (isCustomAminoAcid && (residueStr.length() > 1 || !residueStr.toLowerCase().matches("[bjouxz]"))) {
                System.err.println("Error: Invalid Residue(s) at line " + lineNum +
                        " in file " + sourceFilePath + ": " + modSetting);
                System.err.println("Custom Amino acids are only allowed using B, J, O, U, X, or Z as the custom symbol.");
                return false;
            }
            if (isCustomAminoAcid && !compStr.matches("([CHNOS][0-9]{0,3})+")) {
                System.err.println("Error: Invalid composition/mass at line " + lineNum +
                        " in file " + sourceFilePath + ": " + modSetting);
                System.err.println("Custom Amino acids must supply a composition string, and must not use elements other than C H N O S.");
                return false;
            }
            if (modTypeParseFailed) {
                System.err.println("Error: Modification must be either fix, opt, or custom at line " + lineNum +
                        " in file " + sourceFilePath + ": " + modSetting);
                return false;
            }

            // Location
            Modification.Location location = null;
            String customResidueBase = "";
            String locStr = modInfo[3].trim().split("\\s+")[0].trim();
            if (locStr.equalsIgnoreCase("any"))
                location = Modification.Location.Anywhere;
            else if (locStr.equalsIgnoreCase("N-Term") || locStr.equalsIgnoreCase("NTerm"))
                location = Modification.Location.N_Term;
            else if (locStr.equalsIgnoreCase("C-Term") || locStr.equalsIgnoreCase("CTerm"))
                location = Modification.Location.C_Term;
            else if (locStr.equalsIgnoreCase("Prot-N-Term") || locStr.equalsIgnoreCase("ProtNTerm"))
                location = Modification.Location.Protein_N_Term;
            else if (locStr.equalsIgnoreCase("Prot-C-Term") || locStr.equalsIgnoreCase("ProtCTerm"))
                location = Modification.Location.Protein_C_Term;
            else if (isCustomAminoAcid)
                customResidueBase = locStr;
            else {
                System.err.println("Error: Invalid Location at line " + lineNum +
                        " in file " + sourceFilePath + ": " + modSetting);
                return false;
            }

            String name = modInfo[4].trim().split("\\s+")[0].trim();
            if (!isCustomAminoAcid) {
                if (isModConflict(sourceFilePath, lineNum, modSetting, name, modMass)) {
                    return false;
                }

                Modification mod = Modification.register(name, modMass);

                for (int i = 0; i < residueStr.length(); i++) {
                    char residue = residueStr.charAt(i);
                    Modification.Instance modIns = new Modification.Instance(mod, residue, location);
                    if (isFixedModification) {
                        modIns.fixedModification();
                    }

                    if (!addModInstance(sourceFilePath, lineNum, modSetting, mods, modIns)) {
                        return false;
                    }
                }
            } else {
                char customAminoAcidSymbol = residueStr.charAt(0);

                AminoAcid aa = new AminoAcid(customAminoAcidSymbol, name, new Composition(compStr));
                if (customAAResidues.contains(Character.toString(customAminoAcidSymbol))) {
                    System.err.println(
                            "Error: Duplicate custom amino acid symbol; \n" +
                                    "the duplicate definition is on line " + lineNum +
                                    " in file " + sourceFilePath + ": " + modSetting);
                    return false;
                }
                modMetadata.addCustomAminoAcidSymbol(customAminoAcidSymbol);
                customAA.add(aa);
            }
        }

        return true;
    }

    public static AminoAcidSet getAminoAcidSetFromXMLFile(String modFilePath) {

        File modFile = new File(modFilePath);

        BufferedLineReader reader = null;
        try {
            reader = new BufferedLineReader(modFile.getPath());
        } catch (IOException e) {
            System.err.println("Error opening modification file " + modFile.getPath());
            e.printStackTrace();
            System.exit(-1);
        }

        int numMods = 3;

        // Define keywords
        String numModsKey = "<parameter name=\"ptm.mods\">";
        String cysKey = "<parameter name=\"cysteine_protease.cysteine\">";
        String oxidationKey = "<parameter name=\"ptm.OXIDATION\">on</parameter>";
        String lysMetKey = "<parameter name=\"ptm.LYSINE_METHYLATION\">on</parameter>";
        String pyrogluKey = "<parameter name=\"ptm.PYROGLUTAMATE_FORMATION\">on</parameter>";
        String phosphoKey = "<parameter name=\"ptm.PHOSPHORYLATION\">on</parameter>";
        String ntermCarbamylKey = "<parameter name=\"ptm.NTERM_CARBAMYLATION\">on</parameter>";
        String ntermAcetylKey = "<parameter name=\"ptm.NTERM_ACETYLATION\">on</parameter>";
        String ptmKey = "<parameter name=\"ptm.custom_PTM\">";
        String closeKey = "</parameter>";

        // parse modifications
        ArrayList<Modification.Instance> mods = new ArrayList<>();
        String dataLine;
        int lineNum = 0;
        while ((dataLine = reader.readLine()) != null) {
            lineNum++;
            if (dataLine.startsWith(numModsKey)) {
                try {
                    String value = dataLine.substring(numModsKey.length(), dataLine.lastIndexOf(closeKey));
                    numMods = Integer.parseInt(value);
                } catch (NumberFormatException e) {
                    System.err.println("Error: Invalid ptm.mods option at line " + lineNum +
                            " in file " + modFile.getName() + ": " + dataLine);
                    e.printStackTrace();
                    System.exit(-1);
                }
            } else if (dataLine.startsWith(cysKey)) {
                String value = dataLine.substring(cysKey.length(), dataLine.lastIndexOf(closeKey));
                if (value.equalsIgnoreCase("c57")) {
                    char residue = 'C';
                    Modification mod = Modification.Carbamidomethyl;
                    Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere).fixedModification();
                    if (!addModInstance(modFile.getName(), lineNum, dataLine, mods, modIns)) {
                        System.exit(-1);
                    }
                } else if (value.equalsIgnoreCase("c58")) {
                    char residue = 'C';
                    Modification mod = Modification.Carboxymethyl;
                    Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere).fixedModification();
                    mods.add(modIns);
                } else if (value.equalsIgnoreCase("c99")) {
                    char residue = 'C';
                    Modification mod = Modification.NIPCAM;
                    Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere).fixedModification();
                    mods.add(modIns);
                } else if (value.equalsIgnoreCase("None")) {
                    // do nothing
                } else {
                    System.err.println("Error: Invalid Cysteine protecting group at line " + lineNum +
                            " in file " + modFile.getName() + ": "+ dataLine);
                    System.exit(-1);
                }
            } else if (dataLine.startsWith(ptmKey))    // custom PTM
            {
                String value = dataLine.substring(ptmKey.length(), dataLine.lastIndexOf(closeKey));
                String[] token = value.split(",");

                if (token.length != 3) {
                    System.err.println("Error: Invalid custom ptm option at line " + lineNum +
                            " in file " + modFile.getName() + ": "+ dataLine);
                    System.exit(-1);
                }

                // Mass
                double modMass = 0;
                try {
                    modMass = Double.parseDouble(token[0]);
                } catch (NumberFormatException e) {
                    System.err.println("Error: Invalid Mass at line " + lineNum +
                            " in file " + modFile.getName() + ": "+ dataLine);
                    e.printStackTrace();
                    System.exit(-1);
                }

                // Residues
                String residueStr = token[1];
                boolean isResidueStrLegitimate = true;
                if (!residueStr.equals("*")) {
                    if (residueStr.length() > 0) {
                        for (int i = 0; i < residueStr.length(); i++) {
                            if (!AminoAcid.isStdAminoAcid(residueStr.charAt(i))) {
                                isResidueStrLegitimate = false;
                                break;
                            }
                        }
                    } else
                        isResidueStrLegitimate = false;
                }
                if (!isResidueStrLegitimate) {
                    System.err.println("Error: Invalid Residue(s) at line " + lineNum +
                            " in file " + modFile.getName() + ": " + dataLine);
                    System.exit(-1);
                }

                // Location
                Modification.Location location = null;
                boolean isFixedModification = false;
                String locStr = token[2];

                if (locStr.equalsIgnoreCase("fix")) {
                    isFixedModification = true;
                    location = Location.Anywhere;
                } else if (locStr.equalsIgnoreCase("opt")) {
                    isFixedModification = false;
                    location = Location.Anywhere;
                } else if (locStr.equalsIgnoreCase("opt_nterm")) {
                    isFixedModification = false;
                    location = Location.N_Term;
                } else if (locStr.equalsIgnoreCase("fix_nterm")) {
                    isFixedModification = true;
                    location = Location.N_Term;
                } else if (locStr.equalsIgnoreCase("opt_cterm")) {
                    isFixedModification = false;
                    location = Location.C_Term;
                } else if (locStr.equalsIgnoreCase("fix_cterm")) {
                    isFixedModification = true;
                    location = Location.C_Term;
                } else {
                    System.err.println("Error: Invalid custom_PTM location at line " + lineNum +
                            " in file " + modFile.getName() + ": " + dataLine);
                    System.exit(-1);
                }

                String name = residueStr + " " + modMass;

                if (isModConflict(modFile.getName(), lineNum, dataLine, name, modMass)) {
                    System.exit(-1);
                }

                Modification mod = Modification.register(name, modMass);

                for (int i = 0; i < residueStr.length(); i++) {
                    char residue = residueStr.charAt(i);
                    Modification.Instance modIns = new Modification.Instance(mod, residue, location);
                    if (isFixedModification)
                        modIns.fixedModification();
                    if (!addModInstance(modFile.getName(), lineNum, dataLine, mods, modIns)) {
                        System.exit(-1);
                    }
                }
            } else if (dataLine.startsWith(oxidationKey))    // predefined Oxidized methionine
            {
                String residueStr = "M";
                Modification mod = Modification.Oxidation;
                for (int i = 0; i < residueStr.length(); i++) {
                    char residue = residueStr.charAt(i);
                    Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere);
                    if (!addModInstance(modFile.getName(), lineNum, dataLine, mods, modIns)) {
                        System.exit(-1);
                    }
                }
            } else if (dataLine.startsWith(lysMetKey))    // predefined lysine methylation
            {
                String residueStr = "K";
                Modification mod = Modification.Methyl;
                for (int i = 0; i < residueStr.length(); i++) {
                    char residue = residueStr.charAt(i);
                    Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere);
                    if (!addModInstance(modFile.getName(), lineNum, dataLine, mods, modIns)) {
                        System.exit(-1);
                    }
                }
            } else if (dataLine.startsWith(pyrogluKey))    // predefined pyro glu Q
            {
                String residueStr = "Q";
                Modification mod = Modification.PyroGluQ;
                for (int i = 0; i < residueStr.length(); i++) {
                    char residue = residueStr.charAt(i);
                    Modification.Instance modIns = new Modification.Instance(mod, residue, Location.N_Term);
                    if (!addModInstance(modFile.getName(), lineNum, dataLine, mods, modIns)) {
                        System.exit(-1);
                    }
                }
            } else if (dataLine.startsWith(phosphoKey))    // predefined STY phosphorylation
            {
                String residueStr = "STY";
                Modification mod = Modification.Phospho;
                for (int i = 0; i < residueStr.length(); i++) {
                    char residue = residueStr.charAt(i);
                    Modification.Instance modIns = new Modification.Instance(mod, residue, Location.Anywhere);
                    if (!addModInstance(modFile.getName(), lineNum, dataLine, mods, modIns)) {
                        System.exit(-1);
                    }
                }
            } else if (dataLine.startsWith(ntermCarbamylKey))    // predefined N-terminal carbamylation
            {
                String residueStr = "*";
                Modification mod = Modification.Carbamyl;
                for (int i = 0; i < residueStr.length(); i++) {
                    char residue = residueStr.charAt(i);
                    Modification.Instance modIns = new Modification.Instance(mod, residue, Location.N_Term);
                    if (!addModInstance(modFile.getName(), lineNum, dataLine, mods, modIns)) {
                        System.exit(-1);
                    }
                }
            } else if (dataLine.startsWith(ntermAcetylKey))    // predefined N-terminal acetylation
            {
                String residueStr = "*";
                Modification mod = Modification.Acetyl;
                for (int i = 0; i < residueStr.length(); i++) {
                    char residue = residueStr.charAt(i);
                    Modification.Instance modIns = new Modification.Instance(mod, residue, Location.N_Term);
                    if (!addModInstance(modFile.getName(), lineNum, dataLine, mods, modIns)) {
                        System.exit(-1);
                    }
                }
            }
        }
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(mods);
        aaSet.setMaxNumberOfVariableModificationsPerPeptide(numMods);

        try {
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return aaSet;
    }

    public List<Modification.Instance> getModifications() {
        return modifications;
    }

    /**
     * Gets standard amino acids from file
     *
     * @param aaFilePath amino acid set file name.
     * @return amino acid set object.
     */
    public static AminoAcidSet getAminoAcidSet(String aaFilePath) {
        AminoAcidSet aaSet = new AminoAcidSet();
        BufferedLineReader reader = null;

        File aaFile = new File(aaFilePath);

        try {
            reader = new BufferedLineReader(aaFile.getPath());
        } catch (IOException e) {
            e.printStackTrace();
        }

        String dataLine;
        int lineNum = 0;
        int fileType = 0;    // 0: G,Glycine,57.021464   1: G=57.021463723
        while ((dataLine = reader.readLine()) != null) {
            lineNum++;
            if (dataLine.startsWith("#") || dataLine.length() == 0)
                continue;

            if (fileType == 0 && Character.isDigit(dataLine.charAt(0))) {
                fileType = 1;
                continue;
            }

            AminoAcid aa;
            if (fileType == 0) {
                // Composition is available, e.g.
                // G, Glycine, C2H3N1O1

                String[] token = dataLine.split(",");
                if (token.length != 3) {
                    System.out.println("Ignoring line " + lineNum +
                            " in file " + aaFile.getName() + " since not 3 comma separated fields");
                    continue;
                }

                String residueStr = token[0].trim();
                if (residueStr.length() != 1) {
                    System.err.println("Error: Invalid AASet file format at line " + lineNum +
                            " in file " + aaFile.getName() + " (residue must be a single character): " + dataLine);
                    System.exit(-1);
                }

                char residue = residueStr.charAt(0);
                if (!Character.isUpperCase(residue)) {
                    System.err.println("Error: Invalid AASet file format at line " + lineNum +
                            " in file " + aaFile.getName() + " (residue must be an uppercase letter): " + dataLine);
                    System.exit(-1);
                }
                String name = token[1].trim();

                if (token[2].matches("(C\\d+)*(H\\d+)*(N\\d+)*(O\\d+)*(S\\d+)*")) {
                    // Defined via a composition, e.g. C5H9N1O1S1
                    String compositionStr = token[2].trim();
                    Composition composition = new Composition(compositionStr);
                    aa = AminoAcid.getAminoAcid(residue, name, composition);
                } else {
                    // Not a composition; should be a mass
                    double mass = -1;
                    try {
                        mass = Double.parseDouble(token[2]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid AASet file format at line " + lineNum +
                                " in file " + aaFile.getName() +
                                " (should be a composition like C5H7NO3 or a mass): " + dataLine);
                        System.exit(-1);
                    }
                    aa = AminoAcid.getCustomAminoAcid(residue, name, mass);
                }
            } else {
                // fileType == 1, only masses (and probabilities) are available (e.g. D=115 or D=115,0.0467)
                String[] token = dataLine.split("=");
                if (token.length != 2) {
                    System.err.println("Error: Invalid AASet file format at line " + lineNum +
                            " in file " + aaFile.getName() + " (splitting on = should give 2 items): " + dataLine);
                    System.exit(-1);
                }

                if (token[0].length() != 1) {
                    System.err.println("Error: Invalid AASet file format at line " + lineNum +
                            " in file " + aaFile.getName() + " (amino acid symbol must be a single character): " + dataLine);
                    System.exit(-1);
                }

                if (!Character.isLetter(token[0].charAt(0))) {
                    System.err.println("Error: Invalid AASet file format at line " + lineNum +
                            " in file " + aaFile.getName() + " (amino acid symbol must be a letter): " + dataLine);
                    System.exit(-1);
                }

                char residue = token[0].charAt(0);
                String name = token[0];
                float mass = -1;
                float prob = 0.05f;
                String probabilityAddon = "";

                try {
                    if (!token[1].contains(","))
                        mass = Float.parseFloat(token[1]);
                    else {
                        probabilityAddon = " or probability";
                        mass = Float.parseFloat(token[1].split(",")[0]);
                        prob = Float.parseFloat(token[1].split(",")[1]);
                    }
                } catch (NumberFormatException e) {
                    System.err.println("Invalid AASet file format at line " + lineNum +
                            " in file " + aaFile.getName() +
                            " (NumberFormatException parsing the mass" + probabilityAddon + "): " + dataLine);
                    System.exit(-1);
                }
                if (mass <= 0) {
                    System.err.println("Invalid AASet file format at line " + lineNum +
                            " in file " + aaFile.getName() +
                            " (could not parse the mass" + probabilityAddon + "): " + dataLine);
                    System.exit(-1);
                }
                aa = AminoAcid.getCustomAminoAcid(residue, name, mass).setProbability(prob);
            }
            aaSet.addAminoAcid(aa);
        }
        aaSet.finalizeSet();

        try {
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return aaSet;
    }

    public static AminoAcidSet getStandardAminoAcidSet() {
        if (standardAASet == null) {
            standardAASet = new AminoAcidSet();
            for (AminoAcid aa : AminoAcid.getStandardAminoAcids())
                standardAASet.addAminoAcid(aa);
            standardAASet.finalizeSet();
        }
        return standardAASet;
    }

    public static AminoAcidSet getStandardAminoAcidSetWithFixedCarbamidomethylatedCys() {
        if (standardAASetWithCarbamidomethylatedCys == null) {
            ArrayList<Modification.Instance> mods = new ArrayList<>();
            mods.add(new Modification.Instance(Modification.Carbamidomethyl, 'C').fixedModification());
            standardAASetWithCarbamidomethylatedCys = AminoAcidSet.getAminoAcidSet(mods);
        }
        return standardAASetWithCarbamidomethylatedCys;
    }

    public static AminoAcidSet getStandardAminoAcidSetWithFixedCarboxymethylatedCys() {
        if (standardAASetWithCarboxyomethylatedCys == null) {
            ArrayList<Modification.Instance> mods = new ArrayList<>();
            mods.add(new Modification.Instance(Modification.Carboxymethyl, 'C').fixedModification());
            standardAASetWithCarboxyomethylatedCys = AminoAcidSet.getAminoAcidSet(mods);
        }
        return standardAASetWithCarboxyomethylatedCys;
    }

    /**
     * Creates an alternative amino acid set with the terminal amino acid also
     * encoded.
     *
     * @return the AminoAcidSet with C+57 and X with an arbitrary mass.
     */
    public static AminoAcidSet getStandardAminoAcidSetWithFixedCarbamidomethylatedCysWithTerm() {
        if (standardAASetWithCarbamidomethylatedCysWithTerm == null) {
            Modification.Instance[] mods = {
                    new Modification.Instance(Modification.Carbamidomethyl, 'C').fixedModification()
            };

            HashMap<Character, Modification.Instance> modTable = new HashMap<Character, Modification.Instance>();
            for (Modification.Instance mod : mods) {
                if (mod.isFixedModification()) // variable modifications will be ignored
                    modTable.put(mod.getResidue(), mod);
            }
            AminoAcidSet aaSet = new AminoAcidSet();
            for (AminoAcid aa : AminoAcid.getStandardAminoAcids()) {
                Modification.Instance mod = modTable.get(aa);
                if (mod == null)
                    aaSet.addAminoAcid(aa);
                else
                    aaSet.addAminoAcid(aa.getAAWithFixedModification(mod.getModification()));
            }
            // terminal has 60 has mass, this is arbitrary
//			aaSet.registerAminoAcid(new AminoAcid('X', "STOP", new Composition(2,6,1,1,0)));

            // modified by Sangtae
            aaSet.addAminoAcid(AminoAcid.getCustomAminoAcid('X', new Composition(2, 6, 1, 1, 0).getMass()));

            standardAASetWithCarbamidomethylatedCysWithTerm = aaSet.finalizeSet();
        }
        return standardAASetWithCarbamidomethylatedCysWithTerm;
    }

    public static AminoAcidSet getAminoAcidSet(ArrayList<Modification.Instance> mods) {
        AminoAcidSet aaSet = new AminoAcidSet();
        for (AminoAcid aa : getStandardAminoAcidSet())
            aaSet.addAminoAcid(aa);

        aaSet.applyModifications(mods);
        aaSet.finalizeSet();

        return aaSet;
    }

    public static AminoAcidSet getAminoAcidSet(ArrayList<Modification.Instance> mods, ArrayList<AminoAcid> customAminoAcids) {
        AminoAcidSet aaSet = new AminoAcidSet();
        for (AminoAcid aa : getStandardAminoAcidSet())
            aaSet.addAminoAcid(aa);

        for (AminoAcid aa : customAminoAcids)
            aaSet.addAminoAcid(aa);

        aaSet.applyModifications(mods);
        aaSet.finalizeSet();

        return aaSet;
    }

    public static AminoAcidSet getAminoAcidSet(AminoAcidSet baseAASet, ArrayList<Modification.Instance> mods) {
        AminoAcidSet aaSet = new AminoAcidSet();
        for (AminoAcid aa : baseAASet)
            aaSet.addAminoAcid(aa);

        aaSet.applyModifications(mods);
        aaSet.finalizeSet();

        return aaSet;
    }

    public static AminoAcidSet getAminoAcidSetFromModAAList(AminoAcidSet baseAASet, ArrayList<AminoAcid> modAAList) {
        AminoAcidSet aaSet = new AminoAcidSet();
        for (AminoAcid aa : baseAASet)
            aaSet.addAminoAcid(aa);

        for (AminoAcid aa : modAAList)
            aaSet.addAminoAcid(aa);

        aaSet.finalizeSet();

        return aaSet;
    }

    // returns a new residue for modified amino acid
    private char getModifiedResidue(char unmodifiedResidue) {
        if (!Character.isUpperCase(unmodifiedResidue)) {
            System.err.println("Invalid unmodified residue: " + unmodifiedResidue);
            System.exit(-1);
        }
        // if lowercase letter is available
        char lowerCaseR = Character.toLowerCase(unmodifiedResidue);
        if (!modResidueSet.contains(lowerCaseR)) {
            modResidueSet.add(lowerCaseR);
            return lowerCaseR;
        }

        // if not, use char value >= 128
        char symbol = this.nextResidue;
        nextResidue++;
        if (nextResidue > Character.MAX_VALUE) {
            System.err.println("Too many modifications!");
            System.exit(-1);
        }
        return symbol;
    }

    /**
     * Checks for a conflicting mod definition by modification name
     * @param modFileName Mod file name
     * @param lineNum Line number
     * @param dataLine Text from this line in the mod file
     * @param modName Modification name (case-sensitive)
     * @param modMass Monoisotopic mass
     * @return True if an existing mod is defined with this name but a different mass
     */
    private static boolean isModConflict(
            String modFileName, int lineNum, String dataLine,
            String modName, double modMass) {

        if (!Modification.isModConflict(modName, modMass)) {
            return false;
        }

        // Conflicting mod
        Modification existingMod = Modification.getModByName(modName);

        // Is the user overriding one of the default mods?
        Double existingOverrideMass = defaultModUsage.get(modName);
        if (existingOverrideMass != null) {
            // The mass has already been overridden and a warning has already been shown
            // Make sure the new mass is close to existingOverrideMass
            if (Math.abs(existingOverrideMass.doubleValue() - modMass) <= Modification.MOD_MASS_COMPARISON_THRESHOLD) {
                // Similar masses; no issue
                return false;
            }
        } else {

            for (Modification defaultMod : Modification.getDefaultModList()) {
                if (defaultMod.getName().equals(modName)) {
                    // Warn the user
                    System.out.println(
                            "Warning: Non-standard modification mass defined on line " + lineNum +
                                    " in file " + modFileName + ": " + dataLine);

                    System.out.println("Modification " + modName + " typically has mass " + existingMod.getAccurateMass());
                    System.out.println("Overriding with user-defined value of " + modMass);

                    defaultModUsage.put(modName, modMass);
                    return false;
                }
            }
        }

        System.err.println(
                "Error: Two modifications are defined with the same name but different masses; \n" +
                "the duplicate definition is on line " + lineNum +
                        " in file " + modFileName + ": " + dataLine);

        System.err.println("Modification " + modName + " is already defined with mass " + existingMod.getAccurateMass());
        System.err.println("The duplicate definition has mass " + modMass);
        return true;
    }

    private List<ModifiedAminoAcid> modAAList = new ArrayList<ModifiedAminoAcid>();

    private ModifiedAminoAcid getModifiedAminoAcid(AminoAcid targetAA, Modification.Instance mod) {
        for (ModifiedAminoAcid modAA : modAAList) {
            if (modAA.getTargetAA() == targetAA && modAA.getModification() == mod.getModification())
                return modAA;
        }

        char modResidue = this.getModifiedResidue(targetAA.getUnmodResidue());
        ModifiedAminoAcid modAA = new ModifiedAminoAcid(targetAA, mod, modResidue);
        modAAList.add(modAA);

        return modAA;
    }

    public static void main(String argv[]) {
        ParamManager paramManager = new ParamManager("MS-GF+ AminoAcidSet", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "n/a");
        Path modFilePath = Paths.get(System.getProperty("user.home") + "Research", "Data", "Debug", "mods.txt");
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFilePath.toString(), paramManager);
        aaSet.printAASet();
    }

    private static class ModificationMetadata {
        public ModificationMetadata(int maxNumModsPerPeptide) {
            this.maxNumModsPerPeptide = maxNumModsPerPeptide;
            this.customAAResidues = "";
        }

        public void addCustomAminoAcidSymbol(char customAminoAcidSymbol) {
            customAAResidues += customAminoAcidSymbol;
        }

        public void setMaxNumModsPerPeptide(int newModCount) { maxNumModsPerPeptide = newModCount; }

        // Unused: public void setCustomAAResidues(String residues) { customAAResidues = residues; }

        public int getMaxNumModsPerPeptide() {
            return maxNumModsPerPeptide;
        }

        public String getCustomAAResidues() {
            return customAAResidues;
        }

        int maxNumModsPerPeptide;
        String customAAResidues;

    }
}
