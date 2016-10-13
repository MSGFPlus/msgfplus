package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.misc.CountPSMs;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.IonType;

public class TestMSUtils {
    @Test
    public void getKnownIonTypes()
    {
        for(IonType ionType : IonType.getAllKnownIonTypes(3, true, false, true, true))
        {
            if(ionType.getName().contains("y") && Math.round(ionType.getOffset()) == -227)
                System.out.println(ionType);
        }
    }
    
    @Test
    public void testParsingModFile()
    {
        File modFile = new File("H:\\Research\\Charles_TEDDY\\MSGFDB_Mods.txt");
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
        aaSet.printAASet();
    }
}
