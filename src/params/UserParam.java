package params;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map.Entry;

import parser.BufferedLineReader;

import msutil.ActivationMethod;
import msutil.Enzyme;

public class UserParam {
	private ArrayList<ActivationMethod> actMethodList = new ArrayList<ActivationMethod>();
	private ArrayList<Enzyme> enzList = new ArrayList<Enzyme>();
	
	public void parse() {
		// activation methods
		File actMethodFile = new File("params/activationMethod.txt");
		if(actMethodFile.exists())
		{
			ArrayList<String> paramStrs = parseFromFile(actMethodFile.getPath(), 2);
			for(String paramStr : paramStrs)
			{
				String[] token = paramStr.split(",");
				String shortName = token[0];
				String fullName = token[1];
				actMethodList.add(ActivationMethod.addOrChange(shortName, fullName));
			}
		}
	}
	
	public static ArrayList<String> parseFromFile(String fileName, int tokenLength)
	{
		ArrayList<String> paramStrs = new ArrayList<String>();
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#") || s.length() == 0)
				continue;
			String[] token = s.split(",");
			if(token.length != tokenLength)
				continue;
			else
				paramStrs.add(s);
		}
		return paramStrs;
	}
	
}
