package msgfplus;

import static org.junit.Assert.*;

import java.io.File;
import java.net.URISyntaxException;

import org.junit.Test;

import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;

public class TestPercolator {

    @Test
    public void testAddFeatures() throws URISyntaxException {

        File specFile = new File(TestPercolator.class.getClassLoader().getResource("iprg-2013/F13.mgf").toURI());
        File dbFile = new File(TestPercolator.class.getClassLoader().getResource("iprg-2013/Homo_sapiens_non-redundant.GRCh37.68.pep.all_FPKM-cRAP.fasta").toURI());
        File modFile = new File(TestPercolator.class.getClassLoader().getResource("iprg-2013/Mods.txt").toURI());
        String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), "-addFeatures", "1", "-m", "3"};
        
        ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
        paramManager.addMSGFPlusParams();
        
        String msg = paramManager.parseParams(argv);
        assertTrue(msg == null);
        
        assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
    }

}
