package msgfplus;

import java.io.File;
import java.net.URISyntaxException;

import org.junit.Test;

import edu.ucsd.msjava.msdbsearch.CompactFastaSequence;
import edu.ucsd.msjava.msdbsearch.DBScanner;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.suffixarray.SuffixArray;
import edu.ucsd.msjava.suffixarray.SuffixArraySequence;

public class TestSA {
    @Test
    public void getAAProbabilities() throws URISyntaxException {
        File dbFile = new File(TestSA.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI());
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
        DBScanner.setAminoAcidProbabilities(dbFile.getPath(), aaSet);
        for(AminoAcid aa : aaSet)
        {
            System.out.println(aa.getResidue()+"\t"+aa.getProbability());
        }
    }
    
    @Test
    public void getNumCandidatePeptides() throws URISyntaxException {
        File dbFile = new File(TestSA.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI());
        SuffixArraySequence sequence = new SuffixArraySequence(dbFile.getPath());
        SuffixArray sa = new SuffixArray(sequence);
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(new File(TestSA.class.getClassLoader().getResource("Mods.txt").toURI()).getAbsolutePath());
        System.out.println("NumPeptides: " + sa.getNumCandidatePeptides(aaSet, 2364.981689453125f, new Tolerance(10, true)));
    }

    
    @Test
    public void testRedundantProteins() throws URISyntaxException {
        File databaseFile = new File(TestSA.class.getClassLoader().getResource("ecoli.fasta").toURI());
        
        CompactFastaSequence fastaSequence = new CompactFastaSequence(databaseFile.getPath());
        float ratioUniqueProteins = fastaSequence.getRatioUniqueProteins();
        if(ratioUniqueProteins < 0.5f)
        {
            fastaSequence.printTooManyDuplicateSequencesMessage(databaseFile.getName(), "MS-GF+", ratioUniqueProteins);
            System.exit(-1);
        }
        
        float fractionDecoyProteins = fastaSequence.getFractionDecoyProteins();
        if(fractionDecoyProteins < 0.4f || fractionDecoyProteins > 0.6f)
        {
            System.err.println("Error while reading: " + databaseFile.getName() + " (fraction of decoy proteins: "+ fractionDecoyProteins+ ")");
            System.err.println("Delete " + databaseFile.getName() + " and run MS-GF+ again.");
            System.exit(-1);
        }
        
    }
    
}
