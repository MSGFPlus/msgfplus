package msgfplus;

import org.junit.Test;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.mzid.Unimod;
import edu.ucsd.msjava.mzid.UnimodComposition;

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
		
	}
}
