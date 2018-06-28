package msgfplus;

import edu.ucsd.msjava.msdbsearch.CandidatePeptideGrid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;

import java.io.File;

import static org.junit.Assert.*;

import org.junit.Test;


public class TestCandidatePeptideGrid {

    private void printCandidatePeptideGrid(CandidatePeptideGrid candidatePepGrid) {
        System.out.printf("-------GRID--------\n");
        for (int j = 0; j < candidatePepGrid.size(); j++) {
            System.out.printf("%d : %s\n", j, candidatePepGrid.getPeptideSeq(j));
        }
    }

    @Test
    public void testCandidatePeptideGrid_No_Modified_Residues() {
        System.out.println("Test Unmodified Residues");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aminoAcidSet, Enzyme.TRYPSIN, 3, 8, 1);

        candidatePepGrid.addNTermResidue('A');
        assertEquals("No modifications, so size should stay 1", 1, candidatePepGrid.size());

        candidatePepGrid.addResidue(2, 'A');
        assertEquals("No modifications, so size should stay 1", 1, candidatePepGrid.size());

        candidatePepGrid.addResidue(3, 'A');
        assertEquals("No modifications, so size should stay 1", 1, candidatePepGrid.size());

        assertEquals("Should contain only the peptide AAA", "AAA", candidatePepGrid.getPeptideSeq(0));
    }

    @Test
    public void testCandidatePeptideGrid_Modified_Residues() {
        System.out.println("Test Modified Residues");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aminoAcidSet, Enzyme.TRYPSIN, 3, 8, 1);

        candidatePepGrid.addNTermResidue('S');
        assertEquals("1 variably modified residue, grid size should be 2", 2, candidatePepGrid.size());


        candidatePepGrid.addResidue(2, 'T');
        assertEquals("2 variably modified residues, grid size should be 4", 4, candidatePepGrid.size());

        candidatePepGrid.addResidue(3, 'Y');
        assertEquals("3 variably modified residues, grid size should be 8", 8, candidatePepGrid.size());

        assertEquals("The peptide in position 0 should be the unmodified sequence", "STY", candidatePepGrid.getPeptideSeq(0));
    }

    @Test
    public void testCandidatePeptideGrid_Modified_and_Unmodified_Residues() {
        System.out.println("Test Mixture of Modified and Unmodified Residues");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aminoAcidSet, Enzyme.TRYPSIN, 3, 8, 1);

        candidatePepGrid.addNTermResidue('S');
        assertEquals("1 variably modified residue, grid size should be 2", 2, candidatePepGrid.size());


        candidatePepGrid.addResidue(2, 'A');
        assertEquals("1 variably modified residue, and one unmodified residue, grid size should be 2", 2, candidatePepGrid.size());

        candidatePepGrid.addResidue(3, 'Y');
        assertEquals("2 variably modified residues, and one unmodified residue, grid size should be 4", 4, candidatePepGrid.size());

        assertEquals("The peptide in position 0 should be the unmodified sequence", "SAY", candidatePepGrid.getPeptideSeq(0));
    }

    @Test
    public void testCandidatePeptideGrid_Size_Reset() {
        System.out.println("Test Reusing the Grid for a New Peptide");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aminoAcidSet, Enzyme.TRYPSIN, 3, 8, 1);

        candidatePepGrid.addNTermResidue('S');
        candidatePepGrid.addResidue(2, 'A');
        candidatePepGrid.addResidue(3, 'Y');

        candidatePepGrid.addNTermResidue('A');
        assertEquals("Reusing grid, size should be 1", 1, candidatePepGrid.size());
        assertEquals("Reusing grid, peptide should be 'A'", "A", candidatePepGrid.getPeptideSeq(0));
    }

    @Test
    public void testCandidatePeptideGrid_Missed_Cleavages_CTerm_Enzyme() {
        System.out.println("Test Missed Cleavages - C-term Enzyme");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aminoAcidSet, Enzyme.TRYPSIN, 3, 8, 1);

        candidatePepGrid.addNTermResidue('A');
        assertEquals("First amino acid A when cleaving with Trypsin should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));

        candidatePepGrid.addNTermResidue('K');
        assertEquals("First amino acid K when cleaving with Trypsin should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        candidatePepGrid.addResidue(2, 'R');
        assertEquals("Second amino acid R when cleaving with Trypsin should report 1 missed cleavage for peptide KR", 1, candidatePepGrid.getPeptideNumMissedCleavages(0));

        boolean result = candidatePepGrid.addResidue(3, 'A');
        assertEquals("grid should return false trying to add 'A' to 'KR' because peptide KRA exceeds max 2 missed clavages", false, result);

        result = candidatePepGrid.gridIsOverMaxMissedCleavages(0);
        assertEquals("grid should return true that the peptide it represents exceeds max 2 missed clavages", true, result);
    }

    @Test
    public void testCandidatePeptideGrid_Missed_Cleavages_NTerm_Enzyme() {
        System.out.println("Test Missed Cleavages - N-term Enzyme");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aminoAcidSet, Enzyme.AspN, 3, 8, 1);

        candidatePepGrid.addNTermResidue('D');
        assertEquals("First amino acid D when cleaving with AspN should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));

        candidatePepGrid.addNTermResidue('A');
        assertEquals("First amino acid A when cleaving with AspN should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));

        candidatePepGrid.addResidue(2, 'D');
        assertEquals("Second amino acid D when cleaving with AspN should report 1 missed cleavage for AD", 1, candidatePepGrid.getPeptideNumMissedCleavages(0));

        candidatePepGrid.addResidue(3, 'A');
        assertEquals("Third amino acid A when cleaving with AspN should report 1 missed cleavage for ADA", 1, candidatePepGrid.getPeptideNumMissedCleavages(0));

        boolean result = candidatePepGrid.addResidue(4, 'D');
        assertEquals("grid should return false trying to add 'D' to 'ADA' because it exceeds max 2 missed clavages", false, result);
    }

    @Test
    public void testCandidatePeptideGrid_Missed_Cleavages_NoCleavage_Enzyme() {
        System.out.println("Test Missed Cleavages - NoCleavage");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aminoAcidSet, Enzyme.NoCleavage, 3, 8, 1);

        candidatePepGrid.addNTermResidue('A');
        assertEquals("First amino acid A with no-cleave enzyme should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));

        candidatePepGrid.addNTermResidue('A');
        assertEquals("Second amino acid A with no-cleave enzyme should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));

        candidatePepGrid.addNTermResidue('A');
        assertEquals("Third amino acid A with no-cleave enzyme should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
    }

    @Test
    public void testCandidatePeptideGrid_Missed_Cleavages_Unspecific_Enzyme() {
        System.out.println("Test Missed Cleavages - Unspecific Enzyme");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aminoAcidSet, Enzyme.UnspecificCleavage, 3, 8, 1);

        candidatePepGrid.addNTermResidue('A');
        assertEquals("First amino acid A with unspecific enzyme should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));

        candidatePepGrid.addNTermResidue('A');
        assertEquals("Second amino acid A with unspecific enzyme should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));

        candidatePepGrid.addNTermResidue('A');
        assertEquals("Third amino acid A with unspecific enzyme should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
    }

    @Test
    public void testCandidatePeptideGrid_Missed_Cleavages_Reuse() {
        System.out.println("Test Missed Cleavages When Reusing the Grid - Trypsin");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aminoAcidSet, Enzyme.TRYPSIN, 3, 8, 1);

        /* Use till ti returns false */
        candidatePepGrid.addNTermResidue('K');
        candidatePepGrid.addResidue(2, 'R');
        candidatePepGrid.addResidue(3, 'A');

        /* Reuse, in the middle */
        candidatePepGrid.addResidue(2, 'R');
        assertEquals("grid should return 1 missed cleavages on reuse", 1, candidatePepGrid.getPeptideNumMissedCleavages(0));

        /* Reuse, in the middle */
        candidatePepGrid.addResidue(2, 'R');
        assertEquals("grid should return 1 missed cleavages on reuse", 1, candidatePepGrid.getPeptideNumMissedCleavages(0));

        /* Reuse */
        candidatePepGrid.addNTermResidue('A');
        assertEquals("grid should return 0 missed cleavages on reuse", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));

    }

    @Test
    public void testCandidatePeptideGrid_Missed_Cleavages_No_Limit() {
        System.out.println("Test Missed Cleavages - No Limit on Maximum");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");

        /* Passing -1 for max missed cleavages specifies 'unlimitted' */
        CandidatePeptideGrid candidatePepGrid = new CandidatePeptideGrid(aminoAcidSet, Enzyme.TRYPSIN, 3, 8, -1);

        candidatePepGrid.addNTermResidue('A');
        assertEquals("First amino acid A when cleaving with Trypsin should report 0 missed cleavages", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));

        /* Generate two missed cleavages and test result is still true */
        candidatePepGrid.addNTermResidue('K');
        candidatePepGrid.addResidue(2, 'R');
        boolean result = candidatePepGrid.addResidue(3, 'A');
        assertEquals("grid should return true trying to add 'A' to 'KR' because no limit on number of missed cleavages", true, result);
        result = candidatePepGrid.gridIsOverMaxMissedCleavages(0);
        assertEquals("grid should always return that it is under the max number of allowed missed cleavages", false, result);
    }

}
