package ims;

import java.io.*;
import org.junit.Test;
import edu.ucsd.msjava.msdbsearch.PeptideEnumerator;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.WindowFilter;

public class IMSMiscTest {

	@Test
	public void generateBSAPeptides() throws Exception {

		File bsaFastaFile = new File(IMSMiscTest.class.getClassLoader().getResource("BSA.fasta").toURI());
		File outputFile = File.createTempFile("BSAPeptides", "tsv");
		PeptideEnumerator.enumerate(bsaFastaFile, outputFile);
		outputFile.deleteOnExit();
		System.out.println("Done");
	}
	
	@Test
	public void centroiding() throws Exception
	{

		File outputFile = File.createTempFile("testCentroided", "mgf");
		PrintStream filteredSpecOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		
		File specFile = new File(IMSMiscTest.class.getClassLoader().getResource("test.mgf").toURI());
		SpectraAccessor specAcc = new SpectraAccessor(specFile);
		WindowFilter filter = new WindowFilter(1, 0.5f);
		while(specAcc.getSpecItr().hasNext())
		{
			Spectrum spec = specAcc.getSpecItr().next();
			Spectrum filteredSpec = filter.apply(spec);
			filteredSpec.outputMgf(filteredSpecOut);	
		}
		
		filteredSpecOut.close();
		outputFile.deleteOnExit();
		System.out.println("Done");
	}
	
}
