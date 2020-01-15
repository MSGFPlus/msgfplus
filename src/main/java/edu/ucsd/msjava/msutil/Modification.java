package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.msgf.NominalMass;

import java.util.Comparator;
import java.util.HashMap;

/**
 * A class representing a modification.
 *
 * @author sangtaekim
 */
public class Modification {
    /**
     * Threshold to use when determining if two modifications have the same mass
     */
    public static final double MOD_MASS_COMPARISON_THRESHOLD = 0.01;

    private final String name;
    private final double mass;
    private final int nominalMass;
    private String modId = "";

    /**
     * Empirical formula or modification mass of this modification
     * This is null in certain instances (e.g. custom amino acid residue or non-standard modifications)
     */
    private Composition composition;

    private Modification(String name, Composition composition) {
        this.name = name;
        this.mass = composition.getAccurateMass();
        this.nominalMass = composition.getNominalMass();
        this.composition = composition;
    }

    private Modification(String name, double mass) {
        this.name = name;
        this.mass = mass;
        this.nominalMass = NominalMass.toNominalMass((float) mass);
    }

    /**
     * Modification name
     * @return
     */
    public String getName() {
        return name;
    }

    /**
     * Modification mass (as a float)
     * @return
     */
    public float getMass() {
        return (float) mass;
    }

    /**
     * Modification mass (as a double)
     * @return
     */
    public double getAccurateMass() {
        return mass;
    }

    /**
     * Modification mass (as an integer)
     * @return
     */
    public int getNominalMass() {
        return nominalMass;
    }

    /**
     * Modification identifier (used in mzid output)
     * @return
     */
    public String getModId() {
        return modId;
    }

    /**
     * Empirical formula or modification mass of this modification
     * This is null in certain instances (e.g. custom amino acid residue or non-standard modifications)
     */
    public Composition getComposition() {
        if (composition == null)
            return null;

        return composition;
    }

    /**
     * List of default modifications
     * @return
     */
    public static Modification[] getDefaultModList() {
        return defaultModList;
    }

    /**
     * Looks for an existing mod with the given name
     * @param name Modification name (case-sensitive); getAminoAcidSetFromXMLFile uses 'residueStr + " " + modMass'
     * @param mass Monoisotopic mass
     * @return True if an existing mod exists, and the mass is different (by more than 0.001 Da); otherwise false
     */
    public static boolean isModConflict(String name, double mass) {
        return isModConflict(name, mass, MOD_MASS_COMPARISON_THRESHOLD);
    }

    /**
     * Looks for an existing mod with the given name
     * @param name Modification name (case-sensitive); getAminoAcidSetFromXMLFile uses 'residueStr + " " + modMass'
     * @param mass Monoisotopic mass
     * @return True if an existing mod exists, and the mass is different (by more than massTolerance Da); otherwise false
     */
    public static boolean isModConflict(String name, double mass, double massTolerance) {
        Modification existingMod = modTable.get(name);

        if (existingMod == null)
            return false;

        if (Math.abs(existingMod.mass - mass) > massTolerance)
            return true;

        return false;
    }

    /**
     * Looks for an existing mod with the given name
     * @param name Modification name (case-sensitive)
     * @param composition Modification empirical formula
     * @return True if an existing mod exists, and the mass is different (by more than 0.001 Da); otherwise false
     */
    public static boolean isModConflict(String name, Composition composition) {
        return isModConflict(name, composition.getAccurateMass(), MOD_MASS_COMPARISON_THRESHOLD);
    }

    /**
     * Looks for an existing mod with the given name
     * @param name Modification name (case-sensitive)
     * @param composition Modification empirical formula
     * @return True if an existing mod exists, and the mass is different (by more than massTolerance Da); otherwise false
     */
    public static boolean isModConflict(String name, Composition composition, double massTolerance) {
        return isModConflict(name, composition.getAccurateMass(), massTolerance);
    }

    /**
     * Register a modification by name and mass
     * @param modName Modification name (though getAminoAcidSetFromXMLFile uses 'residueStr + " " + modMass')
     * @param mass Monoisotopic mass
     * @return
     */
    public static Modification register(String modName, double mass) {
        Modification mod = new Modification(modName, mass);
        setModIdentifier(mod);
        modTable.put(modName, mod);
        return mod;
    }

    /**
     * Register a modification by name and composition
     * @param name Modification name
     * @param composition Modification empirical formula
     * @return
     */
    public static Modification register(String name, Composition composition) {
        Modification mod = new Modification(name, composition);
        setModIdentifier(mod);
        modTable.put(name, mod);
        return mod;
    }

    /**
     * Set the mod identifiers for any mods that do not have one.
     * This allows user-specified modifications to take precedence over built-in default modifications
     */
    public static void setModIdentifiers() {
        for (Modification mod : modTable.values()) {
            if (mod.getModId().equals("")) {
                setModIdentifier(mod);
            }
        }
    }

    /**
     * Generate a unique identifier for the modification to be used in mzid output (in peptide and peptideEvidence IDs)
     * @param mod 
     */
    private static void setModIdentifier(Modification mod) {
        double mass = mod.getAccurateMass();
        String baseId = "";
        if (mass >= 0) {
            baseId += "+";
        }
        baseId += Math.round(mod.getAccurateMass());
        String id = baseId;
        int count = 0;
        while (true) {
            boolean foundConflict = false;
            for (Modification existing : modTable.values()) {
                if (existing.modId.equals(id)) {
                    // massMatch: if composition is not null, match on composition; otherwise, match on double-precision mass.
                    boolean massMatch = Composition.equals(existing.composition, mod.composition);
                    if (existing.composition == null) {
                        massMatch = existing.mass == mod.mass;
                    }

                    // If a modification has the same name and composition (or modification mass), give it the same identifier
                    boolean isFullMassMatch = existing.name.equals(mod.name) && massMatch;
                    if (!isFullMassMatch) {
                        foundConflict = true;
                        break;
                    }
                }
            }

            if (!foundConflict) {
                break;
            }

            id = baseId + "#" + (++count);
        }

        mod.modId = id;
    }

    public static Modification getModByName(String name) { return modTable.get(name); }
    public static final Modification Carbamidomethyl = new Modification("Carbamidomethyl", new Composition(2, 3, 1, 1, 0));
    public static final Modification Carboxymethyl = new Modification("Carboxymethyl", new Composition(2, 2, 2, 0, 0));
    public static final Modification NIPCAM = new Modification("NIPCAM", new Composition(5, 9, 1, 1, 0));
    public static final Modification Oxidation = new Modification("Oxidation", new Composition(0, 0, 0, 1, 0));
    public static final Modification Phospho = new Modification("Phospho", Composition.getMass("HO3P"));
    public static final Modification Methyl = new Modification("Methyl", new Composition(1, 2, 0, 0, 0));
    public static final Modification PyroGluQ = new Modification("Gln->pyro-Glu", Composition.getMass("H-3N-1"));    // Pyro-glu from Q
    public static final Modification PyroGluE = new Modification("Glu->pyro-Glu", Composition.getMass("H-2O-1"));    // Pyro-glu from E
    public static final Modification Carbamyl = new Modification("Carbamyl", new Composition(1, 1, 1, 1, 0));
    public static final Modification Acetyl = new Modification("Acetyl", new Composition(2, 2, 0, 1, 0));
    public static final Modification PyroCarbamidomethyl = new Modification("Pyro-carbamidomethyl", Composition.getMass("H-3N-1"));

    // static member
    private static final Modification[] defaultModList =
            {
                    Carbamidomethyl,
                    Carboxymethyl,
                    NIPCAM,
                    Oxidation,
                    Phospho,
                    Methyl,
                    PyroGluQ,
                    PyroGluE,
                    Carbamyl,
                    Acetyl,
                    PyroCarbamidomethyl
            };

    /**
     * Keys are modification names
     * Values are modification details
     */
    private static final HashMap<String, Modification> modTable;

    static {
        modTable = new HashMap<>();
        for (Modification mod : defaultModList) {
            modTable.put(mod.getName(), mod);
        }
    }

    public enum Location {
        Anywhere,
        N_Term,
        C_Term,
        Protein_N_Term,
        Protein_C_Term,
    }

    /**
     * A class representing the modification instance.
     *
     * @author sangtaekim
     */
    public static class Instance {
        private final Modification mod;
        private final char residue;    // if null, no amino acid specificity
        private Location location;    // N_Term, C_Term, Anywhere
        private boolean isFixedModification = false;

        public Instance(Modification mod, char residue, Location location) {
            this.mod = mod;
            this.residue = residue;
            this.location = location;
        }

        public Instance(Modification mod, char residue) {
            this(mod, residue, Location.Anywhere);
        }

        public Instance fixedModification() {
            isFixedModification = true;
            return this;
        }

        public Modification getModification() {
            return mod;
        }

        public char getResidue() {
            return residue;
        }

        public Location getLocation() {
            return location;
        }

        public boolean isFixedModification() {
            return isFixedModification;
        }

        public String toString() {
            return mod.getName() + " " + residue + " " + location + " " + (isFixedModification ? "Fixed (static)" : "Variable (dynamic)");
        }

        @Override
        public boolean equals(Object obj) {
            if (obj instanceof Instance) {
                Instance other = (Instance) obj;
                return mod == other.mod && this.residue == other.residue && this.location == other.location && this.isFixedModification == other.isFixedModification;
            }
            return false;
        }

        @Override
        public int hashCode() {
            return mod.getName().hashCode() + new Character(residue).hashCode() + location.hashCode() + new Boolean(isFixedModification).hashCode();
        }
    }

    public static class MassComparator implements Comparator<Modification> {
        @Override
        public int compare(Modification a, Modification b) {
            return Double.compare(a.getAccurateMass(), b.getAccurateMass());
        }
    }
}
