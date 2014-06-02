package edu.ucsd.msjava.mzid;

import java.io.File;
import java.util.HashSet;

import edu.ucsd.msjava.parser.BufferedLineReader;

public class MzIDTest {
	
	private final File mzidFile;
	
	public MzIDTest(File mzidFile)
	{
		this.mzidFile = mzidFile;
	}
	
	public boolean isValid() throws Exception
	{
		BufferedLineReader in = new BufferedLineReader(mzidFile.getPath());
		
		boolean isValid = true;
		String s;
		HashSet<String> idSet = new HashSet<String>();
		
		String sir = null;
		boolean hasSII = false;
		while((s=in.readLine()) != null)
		{
			s = s.trim();
			
			// duplicate PepEvID
			if(s.startsWith("<PeptideEvidence isDecoy"))
			{
				String id = s.substring(s.lastIndexOf("id=")+4, s.lastIndexOf('"'));
				if(idSet.contains(id))
				{
					System.out.println("Duplicate id: " + id);
					isValid = false;
				}
				else
					idSet.add(id);
			}
			
			// cvRef check
			if(s.startsWith("<cvParam"))
			{
				if(!s.contains("cvRef="))
				{
					System.out.println("No cvRef: " + s);
					isValid = false;
				}
			}
			
			// SIR check
			if(s.startsWith("<SpectrumIdentificationResult"))
			{
				sir = s;
				hasSII = false;
			}
			
			if(sir != null && s.startsWith("<SpectrumIdentificationItem"))
			{
				hasSII = true;
			}
			
			if(s.startsWith("</SpectrumIdentificationResult"))
			{
				if(hasSII == false)
				{
					isValid = false;
					System.out.println("SIR doesn't have SII: " + sir);
				}
				sir = null;
			}
		}
		
		in.close();
		
		return isValid;
	}
	
	public static void main(String argv[]) throws Exception
	{
		File mzidFile = new File("C:\\cygwin\\home\\kims336\\Data\\Debug\\MSGF_Plus_invalid_files_non_unique_key_errors\\testAll.mzid");
		//mzidFile = new File("/Users/kims336/Research/Data/CompRef/Phospho/H20120518_JQ_CPTAC2_COMPREF4_IMAC_02.mzid");
		MzIDTest test = new MzIDTest(mzidFile);
		boolean isValid = test.isValid();
		System.out.println("IsValid? " + isValid);
	}
	
}
