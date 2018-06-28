package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;

public class CandidatePeptideGridConsideringMetCleavage extends CandidatePeptideGrid {

    private final CandidatePeptideGrid candidatePepGridMetCleaved;        // For peptides with Met cleaved
    boolean isProteinNTermWithHeadingMet = false;

//	public CandidatePeptideGridConsideringMetCleavage(AminoAcidSet aaSet, int maxPeptideLength)
//	{
//		this(aaSet, maxPeptideLength, Constants.NUM_VARIANTS_PER_PEPTIDE);
//	}

    public CandidatePeptideGridConsideringMetCleavage(AminoAcidSet aaSet, Enzyme enzyme, int maxPeptideLength, int maxNumVariantsPerPeptide, int maxNumMissedCleavages) {
        super(aaSet, enzyme, maxPeptideLength, maxNumVariantsPerPeptide, maxNumMissedCleavages);
        candidatePepGridMetCleaved = new CandidatePeptideGrid(aaSet, enzyme, maxPeptideLength, maxNumVariantsPerPeptide, maxNumMissedCleavages);
    }

    @Override
    public boolean addProtNTermResidue(char residue) {
        isProteinNTermWithHeadingMet = residue == 'M';
        return super.addProtNTermResidue(residue);
    }

    @Override
    public boolean addNTermResidue(char residue) {
        isProteinNTermWithHeadingMet = false;
        return super.addNTermResidue(residue);
    }

    @Override
    public boolean addResidue(int length, char residue) {
        /* Because of the way the algorithm nests enumerating peptides with
         * and without methionine cleaved, we must consider the case where
         * adding a residue causes more missed cleavages in the peptide
         * that retains the N-term methionine. E.g., if the enzyme is AspN
         * and we have two grids: 'M' and '', and add D to both we get 'MD' and
         * 'D' where the grid with 'MD' now has a missed cleavage and the
         * other with 'D' does not.
         */
        boolean op1 = super.addResidue(length, residue);
        boolean op2 = false;

        if (isProteinNTermWithHeadingMet) {
            if (length == 2)        // Second aa after M (e.g. _.M'G')
                op2 = candidatePepGridMetCleaved.addProtNTermResidue(residue);
            else
                op2 = candidatePepGridMetCleaved.addResidue(length - 1, residue);
        }

        /* Fail once both grids are rejecting extension */
        return op1 || op2;
    }

    @Override
    public boolean addProtCTermResidue(int length, char residue) {
        if (!super.addProtCTermResidue(length, residue))
            return false;

        if (isProteinNTermWithHeadingMet) {
            return candidatePepGridMetCleaved.addProtCTermResidue(length - 1, residue);
        } else
            return true;
    }

    @Override
    public boolean addCTermResidue(int length, char residue) {
        if (!super.addCTermResidue(length, residue))
            return false;

        if (isProteinNTermWithHeadingMet) {
            return candidatePepGridMetCleaved.addCTermResidue(length - 1, residue);
        } else
            return true;
    }

    @Override
    public int size() {
        if (!isProteinNTermWithHeadingMet)
            return super.size();
        else
            return super.size() + candidatePepGridMetCleaved.size();
    }

    @Override
    public boolean isNTermMetCleaved(int index) {
        int sizeNormPep = super.size();
        return index >= sizeNormPep;
    }

    @Override
    public int[] getNominalPRMGrid(int index) {
        if (!isProteinNTermWithHeadingMet)
            return super.getNominalPRMGrid(index);
        int sizeNormPep = super.size();
        if (index < sizeNormPep)
            return super.getNominalPRMGrid(index);
        else
            return candidatePepGridMetCleaved.getNominalPRMGrid(index - sizeNormPep);
    }

    @Override
    public double[] getPRMGrid(int index) {
        if (!isProteinNTermWithHeadingMet)
            return super.getPRMGrid(index);
        int sizeNormPep = super.size();
        if (index < sizeNormPep)
            return super.getPRMGrid(index);
        else
            return candidatePepGridMetCleaved.getPRMGrid(index - sizeNormPep);
    }

    @Override
    public float getPeptideMass(int index) {
        if (!isProteinNTermWithHeadingMet)
            return super.getPeptideMass(index);
        int sizeNormPep = super.size();
        if (index < sizeNormPep)
            return super.getPeptideMass(index);
        else
            return candidatePepGridMetCleaved.getPeptideMass(index - sizeNormPep);
    }

    @Override
    public int getNominalPeptideMass(int index) {
        if (!isProteinNTermWithHeadingMet)
            return super.getNominalPeptideMass(index);
        int sizeNormPep = super.size();
        if (index < sizeNormPep)
            return super.getNominalPeptideMass(index);
        else
            return candidatePepGridMetCleaved.getNominalPeptideMass(index - sizeNormPep);
    }

    @Override
    public String getPeptideSeq(int index) {
        if (!isProteinNTermWithHeadingMet)
            return super.getPeptideSeq(index);
        int sizeNormPep = super.size();
        if (index < sizeNormPep)
            return super.getPeptideSeq(index);
        else
            return candidatePepGridMetCleaved.getPeptideSeq(index - sizeNormPep);
    }

    @Override
    public int getNumMods(int index) {
        if (!isProteinNTermWithHeadingMet)
            return super.getNumMods(index);
        int sizeNormPep = super.size();
        if (index < sizeNormPep)
            return super.getNumMods(index);
        else
            return candidatePepGridMetCleaved.getNumMods(index - sizeNormPep);
    }

    @Override
    public boolean gridIsOverMaxMissedCleavages(int index) {
        /* Protein sequence did not start with methionine */
        if (!isProteinNTermWithHeadingMet)
            return super.gridIsOverMaxMissedCleavages(index);

        /* Protein sequence did begin with methionine, so route the test to the
         * appropriate grid based on the argument index.
         */
        int sizeNormPep = super.size();
        if (index < sizeNormPep)
            return super.gridIsOverMaxMissedCleavages(index);
        else
            return candidatePepGridMetCleaved.gridIsOverMaxMissedCleavages(index - sizeNormPep);
    }

    @Override
    public int getPeptideNumMissedCleavages(int index) {
        /* Protein sequence did not start with methionine */
        if (!isProteinNTermWithHeadingMet)
            return super.getPeptideNumMissedCleavages(index);

        /* Protein sequence did begin with methionine, so route the test to the
         * appropriate grid based on the argument index.
         */
        int sizeNormPep = super.size();
        if (index < sizeNormPep)
            return super.getPeptideNumMissedCleavages(index);
        else
            return candidatePepGridMetCleaved.getPeptideNumMissedCleavages(index - sizeNormPep);
    }
}
