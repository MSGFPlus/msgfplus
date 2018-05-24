package msgfplus;

import java.io.File;
import java.net.URISyntaxException;
import org.junit.Test;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.IonType;

public class TestMSUtils {

    @Test
    public void getKnownIonTypes() {
        for(IonType ionType : IonType.getAllKnownIonTypes(3, true, false, true, true)) {
            if(ionType.getName().contains("y") && Math.round(ionType.getOffset()) == -227)
                System.out.println(ionType);
        }
    }
    
    @Test
    public void testParsingModFile() throws URISyntaxException {
        File modFile = new File(TestMSUtils.class.getClassLoader().getResource("Mods.txt").toURI());
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
        aaSet.printAASet();
    }
}
