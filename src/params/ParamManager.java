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
	
//	public boolean isValid()
//	{
//		Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
//		while(itr.hasNext())
//		{
//			Entry<String, Parameter> entry = itr.next();
//			if(!entry.getValue().isValid())
//				return false;
//		}
//		return true;
//	}

	public boolean isValid()
	{
		Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Parameter> entry = itr.next();
			Parameter param = entry.getValue();
			if(!param.isValid())
			{
				System.err.println("Parameter -" + param.getKey() + " (" + param.getName() + ") is missing!");
				return false;
			}
		}
		return true;
	}
	
	public void printUsageInfo(String prefix)
	{
		Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Parameter> entry = itr.next();
			if(!entry.getValue().isHidden())
				System.out.println(prefix+entry.getValue());
		}
	}
	
	public void printValues()
	{
		Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Parameter> entry = itr.next();
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
		
		FileParameter specFileParam = new FileParameter("s", "SpectrumFile", "*.mzXML, *.mzML, *.mgf, *.ms2, *.pkl or *_dta.txt");
		specFileParam.addExtension("mzXML");
		specFileParam.addExtension("mzML");
		specFileParam.addExtension("mgf");
		specFileParam.addExtension("ms2");
		specFileParam.addExtension("pkl");
		specFileParam.mustExist();
		paramManager.addParameter(specFileParam);
		
		FileParameter dbFileParam = new FileParameter("d", "DatabaseFile", "*.fasta or *.fa");
		dbFileParam.addExtension("fa");
		dbFileParam.addExtension("fasta");
		dbFileParam.mustExist();
		paramManager.addParameter(dbFileParam);
		
		ToleranceParameter pmTolParam = new ToleranceParameter("t", "ParentMassTolerance", "e.g. 2.5Da, 30ppm or 0.5Da,2.5Da");
		paramManager.addParameter(pmTolParam);

		FileParameter outputParam = new FileParameter("o", "DatabaseFile", "*.fasta or *.fa");
		outputParam.setAsOptional();
		outputParam.mustNotExist();
		paramManager.addParameter(outputParam);
		
		IntParameter numThreadParam = new IntParameter("thread", "NumThreads", "Number of concurrent threads to be executed, Default: Number of available cores");
		numThreadParam.defaultValue(Runtime.getRuntime().availableProcessors());
		paramManager.addParameter(numThreadParam);
		
		EnumParameter tdaParam = new EnumParameter("tda");
		tdaParam.registerEntry("don't search decoy database").setDefault();
		tdaParam.registerEntry("search decoy database to compute FDR");
		paramManager.addParameter(tdaParam);
		
		EnumParameter fragParam = new EnumParameter("m", "FragmentMethodID");
		fragParam.registerEntry("as written in the spectrum or CID if no info").setDefault();
		fragParam.registerEntry("CID");
		fragParam.registerEntry("ETD");
		fragParam.registerEntry("HCD");
		fragParam.registerEntry("Merge spectra from the same precursor");
		paramManager.addParameter(fragParam);
		
		EnumParameter instParam = new EnumParameter("inst", "InstrumentID");
		instParam.registerEntry("Low-res LCQ/LTQ").setDefault();
		instParam.registerEntry("TOF");
		instParam.registerEntry("High-res LTQ (Default for HCD)");
		paramManager.addParameter(instParam);

		EnumParameter enzParam = new EnumParameter("e", "EnzymeID");
		enzParam.registerEntry("No enzyme");
		enzParam.registerEntry("Trypsin").setDefault();
		enzParam.registerEntry("Chymotrypsin");
		enzParam.registerEntry("Lys-C");
		enzParam.registerEntry("Lys-N");
		enzParam.registerEntry("Glu-C");
		enzParam.registerEntry("Arg-C");
		enzParam.registerEntry("Asp-N");
		enzParam.registerEntry("alphaLP");
		enzParam.registerEntry("endogenous peptides");
		paramManager.addParameter(enzParam);
		
		EnumParameter c13Param = new EnumParameter("c13");
		c13Param.registerEntry("Consider only peptides matching precursor mass");
		c13Param.registerEntry("Consider peptides having one 13C").setDefault();
		c13Param.registerEntry("Consider peptides having up to two 13C");
		paramManager.addParameter(c13Param);
		
		EnumParameter nnetParam = new EnumParameter("nnet", null, "Number of allowed non-enzymatic termini");
		nnetParam.registerEntry("");
		nnetParam.registerEntry("").setDefault();
		nnetParam.registerEntry("");
		paramManager.addParameter(nnetParam);
		
		FileParameter modParam = new FileParameter("mod", "ModificationFileName", "Modification file, Default: standard amino acids with fixed C+57");
		modParam.setAsOptional();
		modParam.mustExist();
		paramManager.addParameter(modParam);
		
		IntParameter minLenParam = new IntParameter("minLength", "MinPepLength", "Minimum peptide length to consider, Default: 6");
		minLenParam.minValue(1);
		minLenParam.defaultValue(6);
		paramManager.addParameter(minLenParam);

		IntParameter maxLenParam = new IntParameter("maxLength", "MaxPepLength", "Maximum peptide length to consider, Default: 40");
		maxLenParam.minValue(1);
		maxLenParam.defaultValue(40);
		paramManager.addParameter(maxLenParam);

		IntParameter minCharge = new IntParameter("minCharge", "MinCharge", "Minimum precursor charge to consider if charges are not specified in the spectrum file, Default: 2");
		minCharge.minValue(1);
		minCharge.defaultValue(2);
		paramManager.addParameter(minCharge);

		IntParameter maxCharge = new IntParameter("maxCharge", "MaxCharge", "Maximum precursor charge to consider if charges are not specified in the spectrum file, Default: 3");
		maxCharge.minValue(1);
		maxCharge.defaultValue(3);
		paramManager.addParameter(maxCharge);
		
		IntParameter numMatchesParam = new IntParameter("n", "NumMatchesPerSpec", "Number of matches per spectrum to be reported, Default: 1");
		numMatchesParam.minValue(1);
		numMatchesParam.defaultValue(1);
		paramManager.addParameter(numMatchesParam);
		
		EnumParameter uniformAAProb = new EnumParameter("uniformAAProb");
		uniformAAProb.registerEntry("use amino acid probabilities computed from the input database").setDefault();
		uniformAAProb.registerEntry("use probability 0.05 for all amino acids");
		paramManager.addParameter(uniformAAProb);
		
		paramManager.printUsageInfo("\t");
		
		
		if(paramManager.parseParams(argv))
		{
			paramManager.printValues();
			System.out.println("IsValid: " + paramManager.isValid());
		}
	}
}
