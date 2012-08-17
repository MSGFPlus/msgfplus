package edu.ucsd.msjava.mzid;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

public class Unimod {
	private Map<String, String> recordIDMap;	// name -> record id
	public Unimod()
	{
		readUnimodOBOFile();
	}
	
	public String getRecordID(String name)
	{
//		return "UNIMOD";
		return recordIDMap.get(name);
	}
	
	private void readUnimodOBOFile()
	{
//		InputStream is = ClassLoader.getSystemResourceAsStream(UNIMOD_FILE_NAME);
		InputStream is = Unimod.class.getClassLoader().getResourceAsStream(Constants.UNIMOD_RESOURCE_PATH);
		if(is == null)
		{
			System.err.println("Unable to access \"unimod.obo\".");
			System.exit(-1);
		}      
		BufferedReader in = new BufferedReader(new InputStreamReader(is));
		
		recordIDMap = new HashMap<String, String>();
		String s;
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
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
