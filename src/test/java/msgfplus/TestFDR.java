package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.fdr.ComputeFDR;
import edu.ucsd.msjava.fdr.ComputeQValue;

public class TestFDR {
	@Test
	public void testComputeQValue()
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/QCShew");
		File inputFile = new File(dir.getPath()+File.separator+"TestComputeQValue.tsv");
		File outputFile = new File(dir.getPath()+File.separator+"TestComputeQValueWithQValue.tsv");;

		String[] argv = {"-f", inputFile.getPath(), "-o", outputFile.getPath()};
		
		try {
			ComputeQValue.main(argv);
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Done");		
	}
	
	@Test
	public void testPepFDR()
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/Heejung/FDRTest");
		File inputFile = new File(dir.getPath()+File.separator+"NoQWithDecoy.tsv");
		File outputFile = new File(dir.getPath()+File.separator+"Test2NoDecoy.tsv");;

		String[] argv = {"-f", inputFile.getPath(), "10", "XXX", "-i", "0", "-n", "2", "-p", "9", "-s", "13", "0", "-o", outputFile.getPath(), "-decoy", "0"};
		
		try {
			ComputeFDR.main(argv);
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Done");		
	}
	
}
