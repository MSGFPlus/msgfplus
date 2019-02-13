package msgfplus;

import java.io.File;
import java.net.URISyntaxException;

import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
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
        ParamManager paramManager = getParamManager();
        File modFile = new File(TestMSUtils.class.getClassLoader().getResource("Mods.txt").toURI());
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath(), paramManager);
        aaSet.printAASet();
    }

    private ParamManager getParamManager() {
        return new ParamManager("MS-GF+ Test", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "n/a");
    }

}
