package ims;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import org.junit.Test;

import edu.ucsd.msjava.msdbsearch.PeptideEnumerator;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.SpectraContainer;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.WindowFilter;

public class Misc {
	@Test
	public void generateBSAPeptides()
	{
		File bsaFastaFile = new File("C:\\cygwin\\home\\kims336\\Data\\IMS_BSA\\BSA.fasta");
		File outputFile = new File("C:\\cygwin\\home\\kims336\\Data\\SuffixArray\\MSJavaBSAPeptides.tsv");
		try {
			PeptideEnumerator.enumerate(bsaFastaFile, outputFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Done");
	}
	
	@Test
	public void centroiding() throws Exception
	{
		File outputFile = new File("C:\\cygwin\\home\\kims336\\Data\\IMS_Sarc\\testCentroided.mgf");
		PrintStream filteredSpecOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		
		File specFile = new File("C:\\cygwin\\home\\kims336\\Data\\IMS_Sarc\\test.mgf");
		SpectraAccessor specAcc = new SpectraAccessor(specFile);
		WindowFilter filter = new WindowFilter(1, 0.5f);
		while(specAcc.getSpecItr().hasNext())
		{
			Spectrum spec = specAcc.getSpecItr().next();
			Spectrum filteredSpec = filter.apply(spec);
			filteredSpec.outputMgf(filteredSpecOut);	
		}
		
		filteredSpecOut.close();
		System.out.println("Done");
	}
	
}
