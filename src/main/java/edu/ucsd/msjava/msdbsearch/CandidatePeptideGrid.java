package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Modification.Location;

public class CandidatePeptideGrid {
    private static final int STANDARD_RESIDUE_MAX_RESIDUE = 128;

    private final AminoAcidSet aaSet;
    private final Enzyme enzyme;
    private final int maxPeptideLength;
    private final int numMaxMods;
    private final int maxNumVariantsPerPeptide;
    private final int maxNumMissedCleavages;
    private final int[] nMissedCleavages;
    private final char[] residues;
    private final boolean enzymeIsNonSpecific;

    private int[][] nominalPRM;
    private double[][] prm;
    private int[][] numMods;
    private StringBuffer[] peptide;

    // caching amino acid set for fast search

    // anywhere aa (including modified aa)
    private int[][] aaNominalMass; // residue -> mass list
    private double[][] aaMass;
    private char[][] aaResidue;

    // N-term aa set
    private int[][] nTermAANominalMass; // residue -> mass list
    private double[][] nTermAAMass;
    private char[][] nTermAAResidue;

    // C-term aa set
    private int[][] cTermAANominalMass; // residue -> mass list
    private double[][] cTermAAMass;
    private char[][] cTermAAResidue;

    // Protein N-term aa set
    private int[][] protNTermAANominalMass; // residue -> mass list
    private double[][] protNTermAAMass;
    private char[][] protNTermAAResidue;

    // Protein C-term aa set
    private int[][] protCTermAANominalMass; // residue -> mass list
    private double[][] protCTermAAMass;
    private char[][] protCTermAAResidue;

    // Protein N-term Met cleavage
    private int length;
    private int[] size;

//	public CandidatePeptideGrid(AminoAcidSet aaSet, int maxPeptideLength)
//	{
//		this(aaSet, maxPeptideLength, Constants.NUM_VARIANTS_PER_PEPTIDE);
//	}

    public CandidatePeptideGrid(AminoAcidSet aaSet, Enzyme enzyme, int maxPeptideLength, int maxNumVariantsPerPeptide, int maxMissedCleavages) {
        this.numMaxMods = aaSet.getMaxNumberOfVariableModificationsPerPeptide();
        this.maxPeptideLength = maxPeptideLength;
        this.maxNumVariantsPerPeptide = maxNumVariantsPerPeptide;
        this.maxNumMissedCleavages=maxMissedCleavages;
        this.aaSet = aaSet;
        this.enzyme = enzyme;
        this.enzymeIsNonSpecific = enzyme.getName().equals("UnspecificCleavage");

        cacheAASet();

        nominalPRM = new int[maxNumVariantsPerPeptide][maxPeptideLength + 1];
        prm = new double[maxNumVariantsPerPeptide][maxPeptideLength + 1];
        numMods = new int[maxNumVariantsPerPeptide][maxPeptideLength + 1];
        peptide = new StringBuffer[maxNumVariantsPerPeptide];
        size = new int[maxPeptideLength + 1];
        nMissedCleavages = new int[maxPeptideLength + 1];
        residues = new char[maxPeptideLength + 1];

        initializeNTerm();
    }

    private void initializeNTerm() {
        for (int i = 0; i < maxNumVariantsPerPeptide; i++) {
            nominalPRM[i][0] = 0;
            prm[i][0] = 0.;
            numMods[i][0] = 0;
            peptide[i] = new StringBuffer();
        }
        size[0] = 1;
        nMissedCleavages[0]=0;
        residues[0] = '_';
        length = 0;
    }

    public int[] getNominalPRMGrid(int index) {
        return this.nominalPRM[index];
    }

    public double[] getPRMGrid(int index) {
        return this.prm[index];
    }

    public int size() {
        return size[length];
    }

    public float getPeptideMass(int index) {
        return (float) prm[index][length];
    }

    public int getNominalPeptideMass(int index) {
        return nominalPRM[index][length];
    }

    public String getPeptideSeq(int index) {
        return peptide[index].toString();
    }
    
    /**
     * Test whether the peptide currently represented by the grid contains more
     * than the maximum number of allowed missed cleavages.
     * 
     * @param index This parameter is unused, but is necessary because of how
     * this class is extended by CandidatePeptideGridConsideringMetCleavage,
     * which uses the index to route the call to one of two different grids.
     * 
     * @return true for over the maximum number of allowed missed cleavages, 
     * false otherwise.
     * 
     * @see CandidatePeptideGridConsideringMetCleavage
     */
    public boolean gridIsOverMaxMissedCleavages(int index) {
        return maxNumMissedCleavages != -1 && nMissedCleavages[length] > maxNumMissedCleavages;
    }
    
    /**
     * Return the number of missed cleavages in the peptides the grid is
     * representing.
     * 
     * @param index This parameter is unused, but is necessary because of how
     * this class is extended by CandidatePeptideGridConsideringMetCleavage,
     * which uses the index to route the call to one of two different grids.
     * 
     * @return The number of missed cleavages in the current grid peptide 
     * sequence.
     * 
     * @see CandidatePeptideGridConsideringMetCleavage
     */
    public int getPeptideNumMissedCleavages(int index) {
        return nMissedCleavages[length];
    }
    
    public int getNumMods(int index) {
        return numMods[index][length];
    }

//	public boolean addResidue(char residue)
//	{
//		return addResidue(length+1, residue);
//	}

    // if residue is not a standard residue, return false
    public boolean addResidue(int length, char residue) {
        double[] aaMassArr = aaMass[residue];
        if (aaMassArr == null || length > maxPeptideLength)
            return false;

        int[] aaNominalMassArr = aaNominalMass[residue];
        char[] aaResidueArr = aaResidue[residue];

        return addResidue(aaMassArr, aaNominalMassArr, aaResidueArr, length);
    }

    public boolean addProtNTermResidue(char residue) {
        double[] aaMassArr = protNTermAAMass[residue];
        if (aaMassArr == null)
            return false;

        int[] aaNominalMassArr = protNTermAANominalMass[residue];
        char[] aaResidueArr = protNTermAAResidue[residue];

        return addResidue(aaMassArr, aaNominalMassArr, aaResidueArr, 1);
    }

    public boolean addNTermResidue(char residue) {
        double[] aaMassArr = nTermAAMass[residue];
        if (aaMassArr == null)
            return false;

        int[] aaNominalMassArr = nTermAANominalMass[residue];
        char[] aaResidueArr = nTermAAResidue[residue];

        return addResidue(aaMassArr, aaNominalMassArr, aaResidueArr, 1);
    }

    public boolean addProtCTermResidue(int length, char residue) {
        double[] aaMassArr = protCTermAAMass[residue];
        if (aaMassArr == null)
            return false;

        int[] aaNominalMassArr = protCTermAANominalMass[residue];
        char[] aaResidueArr = protCTermAAResidue[residue];

        return addResidue(aaMassArr, aaNominalMassArr, aaResidueArr, length);
    }

    public boolean addCTermResidue(int length, char residue) {
        double[] aaMassArr = cTermAAMass[residue];
        if (aaMassArr == null)
            return false;

        int[] aaNominalMassArr = cTermAANominalMass[residue];
        char[] aaResidueArr = cTermAAResidue[residue];

        return addResidue(aaMassArr, aaNominalMassArr, aaResidueArr, length);
    }

    public boolean isNTermMetCleaved(int index) {
        return false;
    }

    private boolean addResidue(double[] aaMassArr, int[] aaNominalMassArr, char[] aaResidueArr, int length) {
        int parentSize = size[length - 1];
        for (int parentIndex = 0; parentIndex < parentSize; parentIndex++) {
            nominalPRM[parentIndex][length] = nominalPRM[parentIndex][length - 1] + aaNominalMassArr[0];
            prm[parentIndex][length] = prm[parentIndex][length - 1] + aaMassArr[0];
            numMods[parentIndex][length] = numMods[parentIndex][length - 1];
            peptide[parentIndex].setLength(length - 1);
            peptide[parentIndex].append(aaResidueArr[0]);
        }
        size[length] = parentSize;
        
        // modified residue: copy prms up to length-1 into new array
        if (aaMassArr.length > 1 && parentSize < maxNumVariantsPerPeptide) {
            int newIndex = parentSize;
            for (int parentIndex = 0; parentIndex < parentSize; parentIndex++) {
                int numModParent = numMods[parentIndex][length - 1];
                if (numModParent < numMaxMods) {
                    for (int j = 1; j < aaMassArr.length; j++) {
                        for (int k = 1; k < length; k++) {
                            nominalPRM[newIndex][k] = nominalPRM[parentIndex][k];
                            prm[newIndex][k] = prm[parentIndex][k];
                        }
                        peptide[newIndex] = new StringBuffer(peptide[parentIndex].substring(0, length - 1));
                        nominalPRM[newIndex][length] = nominalPRM[newIndex][length - 1] + aaNominalMassArr[j];
                        prm[newIndex][length] = prm[newIndex][length - 1] + aaMassArr[j];
                        numMods[newIndex][length] = numModParent + 1;
                        peptide[newIndex].append(aaResidueArr[j]);
                        newIndex++;
                        if (newIndex >= maxNumVariantsPerPeptide)
                            break;
                    }
                }
                if (newIndex >= maxNumVariantsPerPeptide)
                    break;
            }
            size[length] = newIndex;
        }
        this.length = length;
        
        /* If we are imposing a limit on the maximum number of missed cleavages
         * allowed on candidate peptides.
         */
        if(maxNumMissedCleavages != -1 && !enzymeIsNonSpecific) {
            /* If enzyme cleaves before the amino acid (N-term enzyme), and this
             * is not the first amino acid of the peptide, then it is a missed 
             * cleavage.
             * 
             * E.g., AspN cleaves before D, so peptide YYD has a missed cleavage
             * at position 3, but peptide DYY has no missed cleavages.
             */
            if(enzyme.isCleavable(aaResidueArr[0]) && enzyme.isNTerm() && length > 1) {
                nMissedCleavages[length] = nMissedCleavages[length-1] + 1;
            }

            /* For C-term enzymes, we need to look backward one residue to
             * determine if adding this residue creates a missed cleavage.
             * 
             * E.g., for Trpysin, if the previous residue is K but we are 
             * extending the peptide with another amino acid, the new peptide 
             * has 1 missed cleavage at position length-1 because the K did not
             * cleave.
             */
            else if(enzyme.isCTerm() && enzyme.isCleavable(residues[length-1])) {
                nMissedCleavages[length] = nMissedCleavages[length-1] + 1;
            }

            /* Otherwise, the number of missed cleavages stays the same as the
             * previous peptide. */
            else {
                nMissedCleavages[length] = nMissedCleavages[length-1];
            }
            
            /* Store the look back residue to avoid repeated String parsing */
            residues[length]=aaResidueArr[0];

            /* Return false if the new peptide is over the maximum numer of
             * missed cleavages */
            if(nMissedCleavages[length] > maxNumMissedCleavages) return false;
        }
        
        return true;
    }

    private void cacheAASet() {
        for (Location location : Location.values())
            cacheAASet(location);
    }

    private void cacheAASet(Location location) {
        int[][] stdResidue2NominalMasses = null;
        double[][] stdResidue2Masses = null;
        char[][] stdResidue2Residues = null;

        if (location == Location.Anywhere) {
            stdResidue2NominalMasses = aaNominalMass = new int[STANDARD_RESIDUE_MAX_RESIDUE][];
            stdResidue2Masses = aaMass = new double[STANDARD_RESIDUE_MAX_RESIDUE][];
            stdResidue2Residues = aaResidue = new char[STANDARD_RESIDUE_MAX_RESIDUE][];
        } else if (location == Location.N_Term) {
            stdResidue2NominalMasses = nTermAANominalMass = new int[STANDARD_RESIDUE_MAX_RESIDUE][];
            stdResidue2Masses = nTermAAMass = new double[STANDARD_RESIDUE_MAX_RESIDUE][];
            stdResidue2Residues = nTermAAResidue = new char[STANDARD_RESIDUE_MAX_RESIDUE][];
        } else if (location == Location.C_Term) {
            stdResidue2NominalMasses = cTermAANominalMass = new int[STANDARD_RESIDUE_MAX_RESIDUE][];
            stdResidue2Masses = cTermAAMass = new double[STANDARD_RESIDUE_MAX_RESIDUE][];
            stdResidue2Residues = cTermAAResidue = new char[STANDARD_RESIDUE_MAX_RESIDUE][];
        } else if (location == Location.Protein_N_Term) {
            stdResidue2NominalMasses = protNTermAANominalMass = new int[STANDARD_RESIDUE_MAX_RESIDUE][];
            stdResidue2Masses = protNTermAAMass = new double[STANDARD_RESIDUE_MAX_RESIDUE][];
            stdResidue2Residues = protNTermAAResidue = new char[STANDARD_RESIDUE_MAX_RESIDUE][];
        } else if (location == Location.Protein_C_Term) {
            stdResidue2NominalMasses = protCTermAANominalMass = new int[STANDARD_RESIDUE_MAX_RESIDUE][];
            stdResidue2Masses = protCTermAAMass = new double[STANDARD_RESIDUE_MAX_RESIDUE][];
            stdResidue2Residues = protCTermAAResidue = new char[STANDARD_RESIDUE_MAX_RESIDUE][];
        }

        //for(AminoAcid aa : AminoAcidSet.getStandardAminoAcidSet())
        for (Character aa : aaSet.getResidueListWithoutMods()) {
            //char residue = aa.getResidue();
            char residue = aa.charValue();
            AminoAcid[] aaArr = aaSet.getAminoAcids(location, residue);
            stdResidue2NominalMasses[residue] = new int[aaArr.length];
            stdResidue2Masses[residue] = new double[aaArr.length];
            stdResidue2Residues[residue] = new char[aaArr.length];
            for (int i = 0; i < aaArr.length; i++) {
                stdResidue2NominalMasses[residue][i] = aaArr[i].getNominalMass();
                stdResidue2Masses[residue][i] = aaArr[i].getAccurateMass();
                stdResidue2Residues[residue][i] = aaArr[i].getResidue();
            }
        }
    }
}
