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
		File modFile = new File("C:\\cygwin\\home\\kims336\\Data\\Debug\\MSGF_Plus_invalid_files\\Mods_second_step_toxo.txt");
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
		aaSet.printAASet();
	}
}
