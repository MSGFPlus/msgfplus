package params;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

public class ParamManager {
	private LinkedHashMap<String,Parameter> params;
	
	public ParamManager()
	{
		params = new LinkedHashMap<String,Parameter>();
	}

	public boolean addParameter(Parameter param)
	{
		if(params.containsKey(param.getKey()))
		{
			System.err.println("ParamManager: duplicate key (" + param.getKey() + ")");
			System.exit(-1);
		}
		params.put(param.getKey(), param);
		return true;
	}
	
	public Parameter getParameter(String key)
	{
		return params.get(key);
	}
	
	public boolean isValid()
	{
		Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Parameter> entry = itr.next();
			if(!entry.getValue().isValid())
				return false;
		}
		return true;
	}
	
	public void printUsageInfo(String prefix)
	{
		Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Parameter> entry = itr.next();
			System.out.println(prefix+entry.getValue());
		}
	}
	
	public void printValues()
	{
		Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Parameter> entry = itr.next();
			if(entry.getValue() != null && entry.getValue().isValueAssigned())
				System.out.println(entry.getKey() + "\t" + entry.getValue().getValueAsString());
		}
	}
	
	public boolean parseParams(String argv[])
	{
		if(argv.length < 2 || argv.length % 2 != 0)
		{
			System.err.println("The number of parameters must be even.");
			return false;
		}
		
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length || argv[i].length() <= 1)
			{
				System.err.println("Illegal parameters!");
				return false;
			}
			else
			{
				String key = argv[i].substring(1);
				Parameter param = params.get(key);
				if(param == null)
				{
					System.err.println("Invalid parameter: " + argv[i]);
					return false;
				}
				else 
				{
					if(param.parse(argv[i+1]) == false)
					{
						System.err.println("Invalid value for parameter " + argv[i] + ": " + argv[i+1]);
						return false;
					}
					param.setValueAssigned();
				}
			}
		}
		
		return true;
	}
	
	public static void main(String argv[])
	{
		ParamManager paramManager = new ParamManager();
		
		FileParameter specFileParam = new FileParameter("s", "SpectrumFile", "*.mzXML, *.mzML, *.mgf, *.ms2, *.pkl or *_dta.txt", false);
		specFileParam.addExtension("mzXML");
		specFileParam.addExtension("mzML");
		specFileParam.addExtension("mgf");
		specFileParam.addExtension("ms2");
		specFileParam.addExtension("pkl");
		specFileParam.mustExist();
		paramManager.addParameter(specFileParam);
		
		FileParameter dbFileParam = new FileParameter("d", "DatabaseFile", "*.fasta or *.fa", false);
		dbFileParam.addExtension("fa");
		dbFileParam.addExtension("fasta");
		dbFileParam.mustExist();
		paramManager.addParameter(dbFileParam);
		
		ToleranceParameter pmTolParam = new ToleranceParameter("t", "ParentMassTolerance", "e.g. 2.5Da, 30ppm or 0.5Da,2.5Da", false);
		paramManager.addParameter(pmTolParam);

		FileParameter outputParam = new FileParameter("o", "DatabaseFile", "*.fasta or *.fa", true);
		outputParam.mustNotExist();
		paramManager.addParameter(outputParam);
		
		paramManager.printUsageInfo("\t");
		
		if(paramManager.parseParams(argv))
			paramManager.printValues();
		System.out.println("IsValid: " + paramManager.isValid());
		
	}
	
}
