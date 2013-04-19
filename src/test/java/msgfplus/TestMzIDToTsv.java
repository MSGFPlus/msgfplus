package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.mzid.Unimod;
import edu.ucsd.msjava.mzid.UnimodComposition;
import edu.ucsd.msjava.ui.MzIDToTsv;

public class TestMzIDToTsv {
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
	public void testConversion()
	{
		File dir = new File("C:\\cygwin\\home\\kims336\\Data\\QCShew");
		
		File mzidFile = new File(dir.getPath()+File.separator+"10ppm_TI_1_2_Len50_NumPeak5.mzid");
		File tsvFile = new File(dir.getPath()+File.separator+"test.tsv");
		
		String[] argv = {"-i", mzidFile.getPath(), "-o", tsvFile.getPath(), "-showFormula", "1"};
		MzIDToTsv.main(argv);
		
	}
}
