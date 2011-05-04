package ui;

import java.io.FileNotFoundException;
import java.util.Hashtable;

import parser.BufferedLineReader;

/**
 * This class is for parsing parameter files used in MS-GF, MS-Dictionary and MS-Profile.
 * @author sangtaekim
 *
 */
public class ParameterParser {
	public static class Parameters extends Hashtable<String, String> {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public String getParameter(String name)
		{
			return get(name);
		}
		
		public Integer getIntParameter(String name)
		{
			String param = get(name);
			if(param == null)
				return null;
			else return Integer.parseInt(param);
		}

		public Float getFloatParameter(String name)
		{
			String param = get(name);
			if(param == null)
				return null;
			else return Float.parseFloat(param);
		}
	}
	
	/**
	 * Parses the specified parameter file.
	 * @param fileName the name of the parameter file.
	 * @return A table of parameters.
	 */
	public static Parameters parse(String fileName)
	{
		Parameters params = new Parameters();
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
			String[] token = s.split("=");
			if(token.length != 2)
				continue;
			else
				params.put(token[0].trim(), token[1].trim());
		}
		return params;
	}
}
