package msgfplus;

import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.Ignore;
import org.junit.Test;

import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;

public class TestIPRG {

    @Test
    @Ignore
    public void countProteins()
    {
        String[] accessions = { "P62894", "P00924", "P00330", "P02769"};
        
        File dir = new File("D:\\Research\\Data\\IPRG2014\\20ppm_TI3_NTT2");

        File specFile = new File(dir.getPath()+File.separator+"QC_Shew_12_02_2_1Aug12_Cougar_12-06-11_dta.txt");
        File dbFile = new File(dir.getPath()+File.separator+"ID_003456_9B916A8B.fasta");
        File modFile = new File(dir.getPath()+File.separator+"Mods.txt");
//        File outputFile = new File(dir.getPath()+File.separator+"Test"+"2013-07-26"+".txt");
        String versionString = MSGFPlus.VERSION.split("\\s+")[1];
        versionString = versionString.substring(versionString.indexOf('(')+1, versionString.lastIndexOf(')'));
        String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
                "-mod", modFile.getPath(), "-t", "10ppm", "-tda", "1", "-m", "1", "-ti", "0,1", "-ntt", "1",
                "-o", dir.getPath()+File.separator+"Test_"+versionString+".mzid"
                }; 

        ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
        paramManager.addMSGFPlusParams();
        
        String msg = paramManager.parseParams(argv);
        if(msg != null)
            System.err.println("Error: " + msg);
        assertTrue(msg == null);
        
        assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
    }
}
