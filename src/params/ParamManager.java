package params;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import ui.MSGFDB;

public class ParamManager {
	private LinkedHashMap<String,Parameter> params;
	private String toolName;
	private String version;
	private String date;
	private String command;
	private ArrayList<String> examples = new ArrayList<String>();
	
	public ParamManager(String toolName, String version, String date, String command)
	{
		this.toolName = toolName;
		this.version = version;
		this.date = date;
		this.command = command;
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

	public void addExample(String example)
	{
		this.examples.add(example);
	}
	
	public Parameter getParameter(String key)
	{
		return params.get(key);
	}
	
	public String isValid()
	{
		Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Parameter> entry = itr.next();
			Parameter param = entry.getValue();
			if(!param.isValid())
			{
				return "Parameter -" + param.getKey() + " (" + param.getName() + ") is missing";
			}
		}
		return null;
	}
	
	public void printUsageInfo()
	{
		System.out.println(this.toolName + " v" + this.version + " (" + this.date + ")");
		System.out.println("Usage: " + this.command);
		Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Parameter> entry = itr.next();
			Parameter param = entry.getValue();
			if(!param.isHidden())
			{
				System.out.println("\t"+param);
				if(param.getAdditionalDescription() != null)
					System.out.println("\t   "+param.getAdditionalDescription());
			}
		}
		for(String example : examples)
			System.out.println(example);
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
	
	public String parseParams(String argv[])
	{
		if(argv.length == 0)
		{
			return "No parameter specified.";
		}
		
		if(argv.length < 2 || argv.length % 2 != 0)
		{
			return "The number of parameters must be even.";
		}
		
		for(int i=0; i<argv.length; i+=2)
		{
			if(!argv[i].startsWith("-") || i+1 >= argv.length || argv[i].length() <= 1)
			{
				return "Syntax error.";
			}
			else
			{
				String key = argv[i].substring(1);
				Parameter param = params.get(key);
				if(param == null)
				{
					return "Invalid parameter: " + argv[i] + ".";
				}
				else 
				{
					String error = param.parse(argv[i+1]);
					if(error != null)
					{
						String err = "Invalid value for parameter " + argv[i] + ": " + argv[i+1];
						err += " (" + error + ")"; 
						return err;
					}
					param.setValueAssigned();
				}
			}
		}
		
		String error;
		if((error = isValid()) != null)
			return error;
			
		return null;
	}
	
	public void addSpecFileParam()
	{
		FileParameter specFileParam = new FileParameter("s", "SpectrumFile", "*.mzXML, *.mzML, *.mgf, *.ms2, *.pkl or *_dta.txt");
		specFileParam.addExtension("mzXML");
		specFileParam.addExtension("mzML");
		specFileParam.addExtension("mgf");
		specFileParam.addExtension("ms2");
		specFileParam.addExtension("pkl");
		specFileParam.fileMustExist();
		specFileParam.mustBeAFile();
		addParameter(specFileParam);
	}

	public void addDBFileParam()
	{
		FileParameter dbFileParam = new FileParameter("d", "DatabaseFile", "*.fasta or *.fa");
		dbFileParam.addExtension("fa");
		dbFileParam.addExtension("fasta");
		dbFileParam.fileMustExist();
		dbFileParam.mustBeAFile();
		addParameter(dbFileParam);
	}
	
	public void addPMTolParam()
	{
		ToleranceParameter pmTolParam = new ToleranceParameter("t", "ParentMassTolerance", "e.g. 2.5Da, 30ppm or 0.5Da,2.5Da");
		pmTolParam.setAdditionalDescription("Use comma to set asymmetric values. E.g. \"-t 0.5Da,2.5Da\" will set 0.5Da to the left (expMass<theoMass) and 2.5Da to the right (expMass>theoMass).");
		addParameter(pmTolParam);
	}

	public void addOutputFileParam() 
	{
		FileParameter outputParam = new FileParameter("o", "OutputFile", "*.fasta or *.fa");
		outputParam.setAsOptional();
		outputParam.fileMustNotExist();
		addParameter(outputParam);
	}
	
	public void addFragMethodParam()
	{
		EnumParameter fragParam = new EnumParameter("m", "FragmentMethodID");
		fragParam.registerEntry("as written in the spectrum or CID if no info").setDefault();
		fragParam.registerEntry("CID");
		fragParam.registerEntry("ETD");
		fragParam.registerEntry("HCD");
		fragParam.registerEntry("Merge spectra from the same precursor");
		addParameter(fragParam);
	}
	
	public void addInstTypeParam()
	{
		EnumParameter instParam = new EnumParameter("inst", "InstrumentID");
		instParam.registerEntry("Low-res LCQ/LTQ").setDefault();
		instParam.registerEntry("TOF");
		instParam.registerEntry("High-res LTQ (Default for HCD)");
		addParameter(instParam);
	}
	
	public void addEnzymeParam()
	{
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
		addParameter(enzParam);
	}
	
	public void addModFileParam()
	{
		FileParameter modParam = new FileParameter("mod", "ModificationFileName", "Modification file, Default: standard amino acids with fixed C+57");
		modParam.setAsOptional();
		modParam.fileMustExist();
		addParameter(modParam);
	}
	
	public static void main(String argv[])
	{
		ParamManager paramManager = new ParamManager("MSGFDB", MSGFDB.VERSION, MSGFDB.RELEASE_DATE, "java -Xmx2000M -jar MSGFDB.jar");
		
		paramManager.addSpecFileParam();
		paramManager.addDBFileParam();
		paramManager.addPMTolParam();
		paramManager.addOutputFileParam();
		
		IntParameter numThreadParam = new IntParameter("thread", "NumThreads", "Number of concurrent threads to be executed, Default: Number of available cores");
		numThreadParam.defaultValue(Runtime.getRuntime().availableProcessors());
		numThreadParam.minValue(1);
		paramManager.addParameter(numThreadParam);
		
		EnumParameter tdaParam = new EnumParameter("tda");
		tdaParam.registerEntry("don't search decoy database").setDefault();
		tdaParam.registerEntry("search decoy database to compute FDR");
		paramManager.addParameter(tdaParam);
		
		paramManager.addFragMethodParam();
		paramManager.addInstTypeParam();
		paramManager.addEnzymeParam();
		
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
		
		paramManager.addModFileParam();
		
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
		
		paramManager.addExample("Example (high-precision): java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 30ppm -c13 1 -nnet 0 -tda 1 -o testMSGFDB.tsv");
		paramManager.addExample("Example (low-precision): java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 0.5Da,2.5Da -nnet 0 -tda 1 -o testMSGFDB.tsv");
		
		// Hidden parameters
		FileParameter dbIndexDirParam = new FileParameter("dd", "DBIndexDir", "Path to the directory containing database index files");
		dbIndexDirParam.fileMustExist();
		dbIndexDirParam.mustBeADirectory();
		dbIndexDirParam.setAsOptional();
		dbIndexDirParam.setHidden();
		paramManager.addParameter(dbIndexDirParam);
		
		EnumParameter unitParam = new EnumParameter("u");
		unitParam.registerEntry("Da");
		unitParam.registerEntry("ppm");
		unitParam.registerEntry("Don't care").setDefault();
		unitParam.setHidden();
		paramManager.addParameter(unitParam);
		
		IntRangeParameter specIndexParam = new IntRangeParameter("index", "SpecIndex", "Range of spectrum index to be considered");
		specIndexParam.minValue(1);
		specIndexParam.defaultValue("1,"+(Integer.MAX_VALUE-1));
		specIndexParam.setHidden();
		paramManager.addParameter(specIndexParam);
		
		EnumParameter showFDRParam = new EnumParameter("showFDR");
		showFDRParam.registerEntry("do not show FDRs");
		showFDRParam.registerEntry("show FDRs").setDefault();
		showFDRParam.setHidden();
		paramManager.addParameter(showFDRParam);
		
		EnumParameter replicateMergedResParam = new EnumParameter("replicate");
		replicateMergedResParam.registerEntry("show merged spectra").setDefault();
		replicateMergedResParam.registerEntry("show individual spectra");
		replicateMergedResParam.setHidden();
		paramManager.addParameter(replicateMergedResParam);

		EnumParameter edgeScoreParam = new EnumParameter("edgeScore");
		edgeScoreParam.registerEntry("use edge scoring").setDefault();
		edgeScoreParam.registerEntry("do not use edge scoring");
		edgeScoreParam.setHidden();
		paramManager.addParameter(edgeScoreParam);
		
		String errMessage = paramManager.parseParams(argv); 
		if(errMessage == null)
		{
			paramManager.printValues();
		}
		else
		{
			System.err.println("[Error] " + errMessage);
			System.out.println();
			paramManager.printUsageInfo();
		}
	}
}
