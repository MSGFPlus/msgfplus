package msgfplus;

import edu.ucsd.msjava.msdbsearch.CandidatePeptideGridConsideringMetCleavage;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;

import java.io.File;

import static org.junit.Assert.*;

import org.junit.Test;


public class TestCandidatePeptideGridConsideringMetCleavage {

    private void printCandidatePeptideGridConsideringMetCleavage(CandidatePeptideGridConsideringMetCleavage candidatePepGrid) {
        System.out.printf("-------GRID--------\n");
        for (int j = 0; j < candidatePepGrid.size(); j++) {
            System.out.printf("%d : %s\n", j, candidatePepGrid.getPeptideSeq(j));
        }
    }

    /* Test the expected grid sizes when no modified residues are considered */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_No_Modified_Residues() {
        System.out.println("Test Unmodified Residues");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.TRYPSIN, 4, 8, 1);

        /* Add a methionine, so the size should be 2 when the grid instantiates
         * one grid for generating peptides with methionine and one for ones
         * with methionine cleaved */
        candidatePepGrid.addProtNTermResidue('M');
        assertEquals("Methioinine should cause two grids to be instantiated with initial size 2", 2, candidatePepGrid.size());

        candidatePepGrid.addResidue(2, 'A');
        assertEquals("No modifications, so size should stay 2", 2, candidatePepGrid.size());

        candidatePepGrid.addResidue(3, 'A');
        assertEquals("No modifications, so size should stay 2", 2, candidatePepGrid.size());

        candidatePepGrid.addResidue(4, 'A');
        assertEquals("No modifications, so size should stay 2", 2, candidatePepGrid.size());

        assertEquals("Should contain only the peptide MAAA", "MAAA", candidatePepGrid.getPeptideSeq(0));
        assertEquals("Should contain only the peptide AAA", "AAA", candidatePepGrid.getPeptideSeq(1));
    }

    /* Test the expected grid sizes when only modified residues are considered */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_Modified_Residues() {
        System.out.println("Test Modified Residues");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.TRYPSIN, 4, 8, 1);

        /* Add a methionine, so the size should be 2 when the grid instantiates
         * one grid for generating peptides with methionine and one for ones
         * with methionine cleaved */
        candidatePepGrid.addProtNTermResidue('M');
        assertEquals("Methioinine should cause two grids to be instantiated with initial size 2", 2, candidatePepGrid.size());

        candidatePepGrid.addResidue(2, 'S');
        assertEquals("1 variably modified residue, grid size should be 4", 4, candidatePepGrid.size());

        candidatePepGrid.addResidue(3, 'T');
        assertEquals("2 variably modified residues, grid size should be 8", 8, candidatePepGrid.size());

        candidatePepGrid.addResidue(4, 'Y');
        assertEquals("3 variably modified residues, grid size should be 16", 16, candidatePepGrid.size());

        assertEquals("The peptide in position 0 should be the unmodified sequence", "MSTY", candidatePepGrid.getPeptideSeq(0));
        assertEquals("The peptide in position 8 should be the unmodified sequence", "STY", candidatePepGrid.getPeptideSeq(8));
    }

    /* Test the expected grid sizes when both modified and unmodified residues
     * are considered */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_Modified_and_Unmodified_Residues() {
        System.out.println("Test Mixture of Modified and Unmodified Residues");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.TRYPSIN, 4, 8, 1);

        /* Add a methionine, so the size should be 2 when the grid instantiates
         * one grid for generating peptides with methionine and one for ones
         * with methionine cleaved */
        candidatePepGrid.addProtNTermResidue('M');
        assertEquals("Methioinine should cause two grids to be instantiated with initial size 2", 2, candidatePepGrid.size());

        candidatePepGrid.addResidue(2, 'S');
        assertEquals("1 variably modified residue, grid size should be 4", 4, candidatePepGrid.size());

        candidatePepGrid.addResidue(3, 'A');
        assertEquals("1 variably modified residue, and one unmodified residue, grid size should be 4", 4, candidatePepGrid.size());

        candidatePepGrid.addResidue(4, 'Y');
        assertEquals("2 variably modified residues, and one unmodified residue, grid size should be 8", 8, candidatePepGrid.size());

        assertEquals("The peptide in position 0 should be the unmodified sequence", "MSAY", candidatePepGrid.getPeptideSeq(0));
        assertEquals("The peptide in position 5 should be the unmodified sequence", "SAY", candidatePepGrid.getPeptideSeq(4));
    }

    /* Test that the grid size resets as expected when re-using the grid */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_Size_Reset() {
        System.out.println("Test Reusing the Grid for a New Peptide");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.TRYPSIN, 3, 8, 1);

        candidatePepGrid.addProtNTermResidue('M');
        candidatePepGrid.addNTermResidue('S');
        candidatePepGrid.addResidue(2, 'A');
        candidatePepGrid.addResidue(3, 'Y');

        candidatePepGrid.addProtNTermResidue('M');
        assertEquals("Reusing grid, size should be 2", 2, candidatePepGrid.size());
        assertEquals("Reusing grid, peptide at index 0 should be 'M'", "M", candidatePepGrid.getPeptideSeq(0));
        assertEquals("Reusing grid, peptide at index 1 should be ''", "", candidatePepGrid.getPeptideSeq(1));
    }

    /* Test missed cleavage detection and reporting for the grids including and
     * excluding methionine when using a C-term cleaving enzyme.
     */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_Missed_Cleavages_CTerm_Enzyme() {
        System.out.println("Test Missed Cleavages - C-term Enzyme");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.TRYPSIN, 4, 8, 1);

        candidatePepGrid.addProtNTermResidue('M');

        /* Start out adding a non-cleaving amino acid to verify it returns 0
         * missed cleavages */
        candidatePepGrid.addResidue(2, 'A');
        assertEquals("Adding amino acid A to 'M' when cleaving with Trypsin should report 0 missed cleavages for [M]A", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid A to '' when cleaving with Trypsin should report 0 missed cleavages for A", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));

        /* Start over adding a cleaving amino acid to verify it returns 0
         * missed cleavages */
        candidatePepGrid.addResidue(2, 'K');
        assertEquals("Adding amino acid K to 'M' when cleaving with Trypsin should report 0 missed cleavages for [M]K", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid K to '' when cleaving with Trypsin should report 0 missed cleavages for K", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));

        /* Add another cleaving amino acid, which should turn the previous K
         * into a missed cleavage */
        candidatePepGrid.addResidue(3, 'R');
        assertEquals("Adding amino acid R to 'MK' when cleaving with Trypsin should report 1 missed cleavage for peptides MKR", 1, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid R to K when cleaving with Trypsin should report 1 missed cleavage for peptides KR", 1, candidatePepGrid.getPeptideNumMissedCleavages(1));

        /* Test detection of over max rejecting addition and explict tests for
         * over-max of the methionine and no-methionine grids */
        boolean result = candidatePepGrid.addResidue(4, 'A');
        assertEquals("grid should return false trying to add 'A' to '[M]KR' because peptides [M]KRA exceed max 2 missed cleavages (both grids reject the addition)", false, result);

        result = candidatePepGrid.gridIsOverMaxMissedCleavages(0);
        assertEquals("grid including methionine should return true for overMax after adding 'A' to 'MKR' because peptide MKRA exceeds max 2 missed clavages", true, result);

        result = candidatePepGrid.gridIsOverMaxMissedCleavages(1);
        assertEquals("grid excluding methionine should return true for overMax after adding 'A' to 'KR' because peptide KRA exceeds max 2 missed clavages", true, result);
    }

    /* Test missed cleavage detection and reporting for the grids including and
     * excluding methionine when using an N-term cleaving enzyme.
     */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_Missed_Cleavages_NTerm_Enzyme() {
        System.out.println("Test Missed Cleavages - N-term Enzyme");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.AspN, 5, 8, 1);

        candidatePepGrid.addProtNTermResidue('M');

        /* Start out adding a non-cleaving amino acid to verify it returns 0
         * missed cleavages */
        candidatePepGrid.addResidue(2, 'A');
        assertEquals("Adding amino acid A to 'M' when cleaving with AspN should report 0 missed cleavages for MA", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid A to '' when cleaving with AspN should report 0 missed cleavages for A", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));

        /* Start over adding a cleaving amino acid to verify the grid that
         * includes methionine reports 1 missed cleavage but the grid that
         * excludes methionine reports 0 missed cleavages */
        candidatePepGrid.addResidue(2, 'D');
        assertEquals("Adding amino acid D to 'M' when cleaving with AspN should report 1 missed cleavage for MD", 1, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid D to '' when cleaving with AspN should report 0 missed cleavage for D", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));

        /* Test the success of adding another 'D', and internal divergence of
         * over-max for methionine and non-methionine grids */
        boolean result = candidatePepGrid.addResidue(3, 'D');
        assertEquals("Adding D to should return true because it is under max missed clavages for 'DD' but not for 'MDD' (methionine grid rejected addition, the methionine cleaving grid accepted it)", true, result);
        assertEquals("Adding amino acid D to 'MD' when cleaving with AspN should report 2 missed cleavages for MDD", 2, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid D to 'D' when cleaving with AspN should report 1 missed cleavage for DD", 1, candidatePepGrid.getPeptideNumMissedCleavages(1));

        result = candidatePepGrid.gridIsOverMaxMissedCleavages(0);
        assertEquals("grid including methionine should report that it is over the max number of missed cleavages", true, result);

        result = candidatePepGrid.gridIsOverMaxMissedCleavages(1);
        assertEquals("grid excluding methionine should report that it is NOT over the max number of missed cleavages", false, result);

        /* Test adding an additional missed cleavage triggers rejection by both
         * grids */
        result = candidatePepGrid.addResidue(4, 'D');
        assertEquals("grid should return false trying to add 'D' because both 'MDDD' and 'MDD' exceed max 2 missed cleavages", false, result);
    }

    /* Test missed cleavage detection and reporting for the grids including and
     * excluding methionine when using an unspecific cleaving enzyme.
     */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_Missed_Cleavages_Unspecific_Enzyme() {
        System.out.println("Test Missed Cleavages - Unspecific Enzyme");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.UnspecificCleavage, 5, 8, 1);

        candidatePepGrid.addProtNTermResidue('M');

        /* First amino acid should report 0 missed cleavages */
        candidatePepGrid.addResidue(2, 'A');
        assertEquals("Adding amino acid A to 'M' when cleaving with unspecific enzyme should report 0 missed cleavages for MA", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid A to '' when cleaving with unspecific enzyme should report 0 missed cleavages for A", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));

        /* Second amino acid should report 0 missed cleavages */
        candidatePepGrid.addResidue(3, 'A');
        assertEquals("Adding amino acid A to 'MA' when cleaving with unspecific enzyme should report 0 missed cleavages for MAA", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid A to 'A' when cleaving with unspecific enzyme should report 0 missed cleavages for AA", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));

        /* Third amino acid should report 0 missed cleavages */
        candidatePepGrid.addResidue(3, 'A');
        assertEquals("Adding amino acid A to 'MAA' when cleaving with unspecific enzyme should report 0 missed cleavages for MAAA", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid A to 'AA' when cleaving with unspecific enzyme should report 0 missed cleavages for AAA", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));
    }

    /* Test missed cleavage detection and reporting for the grids including and
     * excluding methionine when using an unspecific cleaving enzyme.
     */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_Missed_Cleavages_NoCleavage_Enzyme() {
        System.out.println("Test Missed Cleavages - NoCleavage Enzyme");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.NoCleavage, 5, 8, 1);

        candidatePepGrid.addProtNTermResidue('M');

        /* First amino acid should report 0 missed cleavages */
        candidatePepGrid.addResidue(2, 'A');
        assertEquals("Adding amino acid A to 'M' when cleaving with no-cleave enzyme should report 0 missed cleavages for MA", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid A to '' when cleaving with no-cleave enzyme should report 0 missed cleavages for A", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));

        /* Second amino acid should report 0 missed cleavages */
        candidatePepGrid.addResidue(3, 'A');
        assertEquals("Adding amino acid A to 'MA' when cleaving with no-cleave enzyme should report 0 missed cleavages for MAA", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid A to 'A' when cleaving with no-cleave enzyme should report 0 missed cleavages for AA", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));

        /* Third amino acid should report 0 missed cleavages */
        candidatePepGrid.addResidue(3, 'A');
        assertEquals("Adding amino acid A to 'MAA' when cleaving with no-cleave enzyme should report 0 missed cleavages for MAAA", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("Adding amino acid A to 'AA' when cleaving with no-cleave enzyme should report 0 missed cleavages for AAA", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));
    }

    /* The grids are instantiated once and reused many times. Test that
     * shortening the peptide in the grid correctly rewinds the number of missed
     * cleavages */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_Missed_Cleavages_Reuse() {
        System.out.println("Test Missed Cleavages When Reusing the Grid");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.TRYPSIN, 3, 8, 1);

        candidatePepGrid.addProtNTermResidue('M');

        /* Use till it returns false */
        candidatePepGrid.addResidue(2, 'K');
        candidatePepGrid.addResidue(3, 'R');
        candidatePepGrid.addResidue(4, 'A');

        /* Reuse, at begining to give 0 missed cleavages */
        candidatePepGrid.addResidue(2, 'R');
        assertEquals("methionine grid should return 0 missed cleavages on reuse", 0, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("grid should return 0 missed cleavages on reuse", 0, candidatePepGrid.getPeptideNumMissedCleavages(1));

        /* Add residue after R to trigger missed cleavage */
        candidatePepGrid.addResidue(3, 'A');
        assertEquals("methionine grid should return 1 missed cleavages on reuse", 1, candidatePepGrid.getPeptideNumMissedCleavages(0));
        assertEquals("grid should return 1 missed cleavages on reuse", 1, candidatePepGrid.getPeptideNumMissedCleavages(1));
    }

    /* Specifying -1 for max missed clavages specifies 'unlimited' which can
     * be used as default behavior for backward compatibility */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_Missed_Cleavages_No_Limit() {
        System.out.println("Test Missed Cleavages - No Limit on Maximum");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");

        /* Passing -1 for max missed cleavages specified 'unlimited' */
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.TRYPSIN, 4, 8, -1);

        /* Generate two missed cleavages and test result is still true */
        candidatePepGrid.addProtNTermResidue('M');
        candidatePepGrid.addResidue(2, 'K');
        candidatePepGrid.addResidue(3, 'R');
        boolean result = candidatePepGrid.addResidue(4, 'A');
        assertEquals("grid should return true trying to add 'A' to 'KR' because no limit on number of missed cleavages", true, result);

        result = candidatePepGrid.gridIsOverMaxMissedCleavages(0);
        assertEquals("methionine grid should always return that it is under the max number of allowed missed cleavages", false, result);

        result = candidatePepGrid.gridIsOverMaxMissedCleavages(1);
        assertEquals("grid should always return that it is under the max number of allowed missed cleavages", false, result);
    }

    /* Specifying -1 for max missed clavages specifies 'unlimited' which can
     * be used as default behavior for backward compatibility */
    @Test
    public void testCandidatePeptideGridConsideringMetCleavage_No_Missed_Cleavages_Allowed() {
        System.out.println("Test Missed Cleavages - No Limit on Maximum");
        File dir = new File("src/test/resources");
        AminoAcidSet aminoAcidSet = AminoAcidSet.getAminoAcidSetFromModFile(dir.getPath() + File.separator + "mods/TestCandidatePeptideGrid.txt");

        /* Passing -1 for max missed cleavages specified 'unlimited' */
        CandidatePeptideGridConsideringMetCleavage candidatePepGrid = new CandidatePeptideGridConsideringMetCleavage(aminoAcidSet, Enzyme.TRYPSIN, 4, 8, 0);

        /* Generate two missed cleavages and test result is still true */
        candidatePepGrid.addProtNTermResidue('M');
        boolean result = candidatePepGrid.addResidue(2, 'K');
        assertEquals("grid should return true trying to add 'K' to '[M]' because [M]K has no missed cleavages", true, result);

        result = candidatePepGrid.addResidue(3, 'A');
        assertEquals("grid should return false trying to add 'A' to '[M]KA' because [M]KA has one missed cleavage", false, result);

        result = candidatePepGrid.gridIsOverMaxMissedCleavages(0);
        assertEquals("methionine grid should always return that it is over max number of allowed missed cleavages", true, result);

        result = candidatePepGrid.gridIsOverMaxMissedCleavages(1);
        assertEquals("grid should always return that it is over the max number of allowed missed cleavages", true, result);
    }

}
