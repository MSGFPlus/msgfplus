package ipa;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.ucsd.msjava.parser.BufferedLineReader;

public class MSGFPlusResultMap {
	private List<PSM> psmList;
	
	public MSGFPlusResultMap(File msgfPlusResultFile)
	{
		parse(msgfPlusResultFile);
	}
	
	public List<PSM> getPSMList() { return psmList; }
	
	private void parse(File msgfPlusResultFile)
	{
		psmList = new ArrayList<PSM>();
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(msgfPlusResultFile.getPath());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String s;
		
		in.readLine();	// header
		while((s=in.readLine()) != null)
		{
			PSM psm = new PSM(s);
			psmList.add(psm);
		}
		
		try {
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
