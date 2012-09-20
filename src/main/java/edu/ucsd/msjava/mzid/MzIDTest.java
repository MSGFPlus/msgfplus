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
	
	public boolean pepEvDupCheck() throws Exception
	{
		BufferedLineReader in = new BufferedLineReader(mzidFile.getPath());
		
		boolean dupID = false;
		String s;
		HashSet<String> idSet = new HashSet<String>();
		while((s=in.readLine()) != null)
		{
			s = s.trim();
			if(s.startsWith("<PeptideEvidence isDecoy"))
			{
				String id = s.substring(s.lastIndexOf("id=")+4, s.lastIndexOf('"'));
				if(idSet.contains(id))
				{
					System.out.println("Duplicate id: " + id);
					dupID = true;
				}
				else
					idSet.add(id);
			}
		}
		
		in.close();
		
		return dupID;
	}
	
	public static void main(String argv[]) throws Exception
	{
		File mzidFile = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\QCShew\\QC_Shew_MSGFPlus_0_1.mzid");
		MzIDTest test = new MzIDTest(mzidFile);
		System.out.println("Is there duplicated PepEvID? " + test.pepEvDupCheck());
	}
	
}
