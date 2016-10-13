package msgfplus;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;

public class TestPercolator {
    @Test
    public void testAddFeatures()
    {
        File dir = new File("/Users/kims336/Research/Data/QCShew/Percolator");

//        File specFile = new File(dir.getPath()+File.separator+"testCID.mgf");
//        File dbFile = new File(dir.getPath()+File.separator+"testCID.fasta");
//        File modFile = new File(dir.getPath()+File.separator+"testCIDMods.txt");
//        String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), "-addFeatures", "1", "-mod", modFile.getPath()};

        File specFile = new File(dir.getPath()+File.separator+"test.mgf");
        File dbFile = new File(dir.getPath()+File.separator+"test.fasta");
        String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), "-addFeatures", "1", "-m", "3"};
        
        ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
        paramManager.addMSGFPlusParams();
        
        String msg = paramManager.parseParams(argv);
        assertTrue(msg == null);
        
        assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
    }

}
