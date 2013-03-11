package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.fdr.ComputeFDR;

public class TestFDR {
	@Test
	public void testMultipleNTermMod()
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/Heejung/FDRTest");
		File inputFile = new File(dir.getPath()+File.separator+"NoQWithDecoy.tsv");
		File outputFile = new File(dir.getPath()+File.separator+"Test2.tsv");;

		String[] argv = {"-f", inputFile.getPath(), "10", "XXX", "-i", "0", "-n", "2", "-p", "9", "-s", "13", "0", "-o", outputFile.getPath(), "-decoy", "1"};
		
		try {
			ComputeFDR.main(argv);
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Done");		
	}
}
