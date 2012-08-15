package edu.ucsd.msjava.psi;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

public class Unimod {
	private Map<String, String> recordIDMap;	// name -> record id
	
	public Unimod(InputStream is)
	{
		BufferedReader in = new BufferedReader(new InputStreamReader(is));
		readUnimodOBOFile(in);
	}
	
	public String getRecordID(String name)
	{
		return recordIDMap.get(name);
	}
	
	private void readUnimodOBOFile(BufferedReader in)
	{
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
