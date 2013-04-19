package edu.ucsd.msjava.mzid;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

public class Unimod {
	
	public static Unimod getUnimod()
	{
		return unimod;
	}
	
	public String getRecordID(String name)
	{
		return recordIDMap.get(name);
	}
	
	public String getDeltaComposition(String id)
	{
		return idToDeltaCompositionMap.get(id);
	}

	private Map<String, String> recordIDMap;	// name -> record id
	private Map<String, String> idToDeltaCompositionMap;	// id -> delta_composition

	private Unimod()
	{
		readUnimodOBOFile();
	}
	
	private void readUnimodOBOFile()
	{
		InputStream is = Unimod.class.getClassLoader().getResourceAsStream(Constants.UNIMOD_RESOURCE_PATH);
		if(is == null)
		{
			System.err.println("Unable to access \"unimod.obo\".");
			System.exit(-1);
		}      
		BufferedReader in = new BufferedReader(new InputStreamReader(is));
		
		recordIDMap = new HashMap<String, String>();
		idToDeltaCompositionMap = new HashMap<String, String>();
		String s;
		String curID = null;
		try {
			while((s=in.readLine()) != null)
			{
				if(s.startsWith("id:"))
				{
					String id = s.split("\\s+")[1].trim();
					String nameLine = in.readLine();
					assert(nameLine.startsWith("name:"));
					String name = nameLine.split("\\s+")[1].trim();
					recordIDMap.put(name, id);
					curID = id;
				}
				if(s.startsWith("xref: delta_composition"))
				{
					String deltaComposition = s.substring(s.indexOf('"')+1, s.lastIndexOf('"'));
					idToDeltaCompositionMap.put(curID, deltaComposition);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}	
	
	private static Unimod unimod;
	
	static {
		unimod = new Unimod();
	}
}
