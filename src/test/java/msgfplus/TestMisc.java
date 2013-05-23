package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.misc.ConvertToMgf;
import edu.ucsd.msjava.misc.VennDiagram;

public class TestMisc {
	@Test
	public void testMisc()
	{
		String title = "Scan:25485 RT:62.983 PrecursorScan:25482 nMSN:19700 PrecursorMonoisoMZ:1134.2547 PEPMASS:Monoiso PrecursorMZ:1134.5891 PrecursorCharge:3 PrecursorScanFTMS:1 FTResolution:17500 IBP:3750861.85 ITot:126255817.49 max2med:29.41 InjTime:31.98 HCD=54.0063972473145eV IsolationMZ:1134.5900 PrecursorAb:8797869.00 MPY:1.00 ms1PrecursorTotAb:47857048671.30 ms1PrecursorInjTime:0.26 ms1PrecursorMZ:1134.5891 ms1PrecursorMzAvg:1134.8744 ms1PrecursorMzRMS:0.3998 ms1PrecursorIntens:8797869.00 ms1PrecursorRT:62.974 ms2IsolationWidth:2.50 ms1SelMZ:1134.2067-1135.5900 ms1SelAvgMZ:1134.7794 ms1SelRmsMZ:0.0636 PrecursorHasMax:1 ms1PrecursorAb:110982284.68 ms1PrecursorMax:111444926.50 numOCMF:270,270,0,7 PrecursorMaxMZ:1134.9262 PrecursorMaxAb:33904613.19 PrecursorMaxRT:62.988 PrecursorWayMMF:0.71 PrecursorMaxMMF:0.71 mzRmsMax:1.33 mzRmsMs2:1.32 maxDelMz:0.3310,5-0,99,99 ms2DelMz:0.3309,5-0,99,99 FilterMzPeakExists(25482):1 PCFD2,2;1,0,1134.2562,3,5,0.3331,0.0025,99,1,14,1.04,0.18,10,10,1.8,3.2,7.5,961;1,0,1134.5877,1,2,0.9994,0.0000,41,1,25,1.25,0.42,16,16,1.8,3.2,7.5,915 Precursor1HasMax:1 Precursor1MaxInjTime:0.26 Precursor1MaxTotAb:47857049600.00 Precursor1MaxAb:108699575.31 Precursor1MaxRT:62.974 Precursor1MaxWidth:0.982 Precursor1MaxWid50:0.174 Precursor1MaxRatio:1.0117 Precursor1MaxBkg:0.00 Precursor1AbuBkg:0.00 Precursor1MaxHW:0.50 Precursor1MaxSkew:-0.00 Precursor2HasMax:1 Precursor2MaxInjTime:0.26 Precursor2MaxTotAb:47857049600.00 Precursor2MaxAb:108699575.31 Precursor2MaxRT:62.974 Precursor2MaxWidth:0.960 Precursor2MaxWid50:0.174 Precursor2MaxRatio:1.0117 Precursor2MaxBkg:0.00 Precursor2AbuBkg:0.00 Precursor2MaxHW:0.50 Precursor2MaxSkew:-0.00 PrecursorMaxNoise:2.28 PrecursorRTStep:0.012 ConvVer:20120705a NumPeaks:472 Filter:FTMS + p NSI d Full ms2 1134.59@hcd28.00<mailto:1134.59@hcd28.00> [100.00-3495.00]";
		System.out.println(title.matches("^Scan:\\d+\\s.+"));
		System.out.println(title.matches("^Scan:\\d+\\sRT:\\d+\\.\\d+\\s.+"));
		System.out.println(title.matches("^Scan:\\d+\\sRT:\\d+\\.\\d+\\sPrecursorScan:\\d+\\??\\s.+"));
		String[] token = title.split("\\s+");
		int scanNum = Integer.parseInt(token[0].substring("Scan:".length()));
		System.out.println(scanNum);
	}
	
	@Test
	public void testVennDiagram()
	{
		File result1 = new File("/Users/kims336/Research/Data/Tao/Global/MSGFPlus_10ppm_TI1/CPTAC_OvC_JB5427_iTRAQ_01_9Apr12_Cougar_12-03-21_dta.tsv");
		File result2 = new File("/Users/kims336/Research/Data/Tao/Global/MSGFPlus_20ppm_TI2/CPTAC_OvC_JB5427_iTRAQ_01_9Apr12_Cougar_12-03-21_dta.tsv");
		
		try {
			VennDiagram.vennDiagram(result1, result2, 0.01f);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	@Test
	public void testMzMLParser()
	{
		File dir = new File("/Users/kims336/Research/Data/ASMS2013");
		File mzMLFile = new File(dir.getPath()+File.separator+"mzMLNoRefinement/QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.mzML");
		File mgfFile = new File(dir.getPath()+File.separator+"IPA/testIPA.mgf");
		try {
			ConvertToMgf.convert(mzMLFile, mgfFile, false, null, null, 10020, -1, true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
