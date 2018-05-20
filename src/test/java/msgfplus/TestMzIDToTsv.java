package msgfplus;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.List;

import org.junit.Ignore;
import org.junit.Test;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.mzid.Unimod;
import edu.ucsd.msjava.mzid.UnimodComposition;
import edu.ucsd.msjava.parser.TSVParser;
import edu.ucsd.msjava.ui.MzIDToTsv;

public class TestMzIDToTsv {

    @Test
    public void testConversionError() throws URISyntaxException, IOException {

        File mzidFile = new File(TestMzIDToTsv.class.getClassLoader().getResource("test.mzid").toURI());

        File tsvFile = File.createTempFile("test", ".mgf");
        
        String[] argv = {"-i", mzidFile.getPath(), "-o", tsvFile.getPath()};
        MzIDToTsv.main(argv);

        tsvFile.deleteOnExit();
    }

    @Test
    public void testReadingUnimodCompositions()
    {
        String deltaComposition = Unimod.getUnimod().getDeltaComposition("UNIMOD:1379");
        System.out.println(deltaComposition);
        
        UnimodComposition comp = new UnimodComposition();
        comp.add(AminoAcidSet.getStandardAminoAcidSet().getAminoAcid('K').getComposition());
        comp.add(deltaComposition);
        comp.add(956.322026);
        System.out.println(comp);
    }
    
    @Test
    public void testConversion() throws IOException, URISyntaxException {

        File mzidFile = new File(TestMzIDToTsv.class.getClassLoader().getResource("test.mzid").toURI());

        File tsvFile = File.createTempFile("test", ".mgf");
        String[] argv = {"-i", mzidFile.getPath(), "-o", tsvFile.getPath(), "-showFormula", "1"};
        MzIDToTsv.main(argv);
    }
    
    @Test
    public void unimodReaderTest()
    {
        String compStr = "H(53) B C(37) N(6) O(6) F(2) S";
        System.out.println(UnimodComposition.getMass(compStr));
    }
    
    @Test
    @Ignore

    public void testMolecularFormula()
    {
        File tsvFile = new File("\\\\protoapps\\UserData\\Sangtae\\Examples\\EmpFormulaExample.tsv");
        TSVParser parser = new TSVParser();
        parser.parse(tsvFile.getPath());
        List<String> pepList = parser.getList("Peptide");
        List<String> formulaList = parser.getList("Formula");
        assertTrue(pepList.size() == formulaList.size());
        AminoAcidSet stdAASet = AminoAcidSet.getStandardAminoAcidSet();
        for(int i=0; i<pepList.size(); i++)
        {
            String pepStr = pepList.get(i);
            String formula = formulaList.get(i);
            Peptide peptide = new Peptide(pepStr, stdAASet);
            double mass1 = peptide.getAccurateMass() + Composition.H2O;
            double mass2 = UnimodComposition.getMass(formula);
            if(Math.abs(mass1-mass2) > 0.01f)
            {
                System.out.println(pepStr+"\t"+mass1);
                System.out.println(formula+"\t"+mass2);
                System.exit(-1);
            }
        }
    }
    
    @Test
    public void testReadingUnimod()
    {
        Unimod.getUnimod();
    }
}
