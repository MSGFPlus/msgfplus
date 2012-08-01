package params;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import msutil.ActivationMethod;
import msutil.DBFileFormat;
import msutil.Enzyme;
import msutil.FileFormat;
import msutil.InstrumentType;
import msutil.Protocol;
import msutil.SpecFileFormat;

import ui.MSGF;

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
	
	public void printToolInfo()
	{
		System.out.println(this.toolName + " v" + this.version + " (" + this.date + ")");
	}
	
	public void printUsageInfo()
	{
		System.out.println(this.toolName + " v" + this.version + " (" + this.date + ")");
		System.out.println("Usage: " + this.command);
		
		ArrayList<Parameter> optParams = new ArrayList<Parameter>();
		Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Parameter> entry = itr.next();
			Parameter param = entry.getValue();
			if(!param.isHidden())
			{
				if(!param.isOptional())
				{
					System.out.println("\t"+param);
					if(param.getAdditionalDescription() != null)
						System.out.println("\t   "+param.getAdditionalDescription());
				}
				else
					optParams.add(param);
			}
		}
		
		for(Parameter param : optParams)
		{
			System.out.println("\t"+param);
			if(param.getAdditionalDescription() != null)
				System.out.println("\t   "+param.getAdditionalDescription());
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
			Parameter param = entry.getValue();
			System.out.println(param.getKey() + "\t" + param.getValueAsString());
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
		FileParameter specFileParam = new FileParameter("s", "SpectrumFile", "*.mzML, *.mzXML, *.mgf, *.ms2, *.pkl or *_dta.txt");
		specFileParam.addFileFormat(SpecFileFormat.MZML);
		specFileParam.addFileFormat(SpecFileFormat.MZXML);
		specFileParam.addFileFormat(SpecFileFormat.MGF);
		specFileParam.addFileFormat(SpecFileFormat.MS2);
		specFileParam.addFileFormat(SpecFileFormat.PKL);
		specFileParam.addFileFormat(SpecFileFormat.DTA_TXT);
		specFileParam.addFileFormat(FileFormat.DIRECTORY);
		specFileParam.fileMustExist();
		addParameter(specFileParam);
	}

	public void addDBFileParam()
	{
		addDBFileParam("d", "*.fasta or *.fa", false);
	}
	
	public void addDBFileParam(String key, String description, boolean isOptional)
	{
		FileParameter dbFileParam = new FileParameter(key, "DatabaseFile", description);
		if(isOptional)
			dbFileParam.setAsOptional();
		dbFileParam.addFileFormat(DBFileFormat.FASTA);
		dbFileParam.fileMustExist();
		dbFileParam.mustBeAFile();
		addParameter(dbFileParam);
	}
	
	public void addPMTolParam()
	{
		ToleranceParameter pmTolParam = new ToleranceParameter("t", "ParentMassTolerance", "e.g. 2.5Da, 30ppm or 0.5Da,2.5Da");
		pmTolParam.setAdditionalDescription("Use comma to set asymmetric values. E.g. \"-t 0.5Da,2.5Da\" will set 0.5Da to the left (expMass<theoMass) and 2.5Da to the right (expMass>theoMass)");
		addParameter(pmTolParam);
	}

	public void addMzIdOutputFileParam() 
	{
		FileParameter outputParam = new FileParameter("o", "OutputFile (*.mzid)", "Default: SpectrumFileName.mzid");
		outputParam.addFileFormat(new FileFormat(".mzid").setCaseSensitive());
		outputParam.setAsOptional();
		outputParam.fileMustNotExist();
		addParameter(outputParam);
	}
	
	public void addOutputFileParam() 
	{
		FileParameter outputParam = new FileParameter("o", "OutputFile", "Default: stdout");
		outputParam.setAsOptional();
		outputParam.fileMustNotExist();
		addParameter(outputParam);
	}
	
	public void addFragMethodParam()
	{
		addFragMethodParam(ActivationMethod.ASWRITTEN, false);
	}
	
	public void addFragMethodParam(ActivationMethod defaultMethod, boolean doNotAddMergeMode)
	{
		ObjectEnumParameter<ActivationMethod> fragParam = new ObjectEnumParameter<ActivationMethod>("m", "FragmentMethodID");
		ActivationMethod[] methods = ActivationMethod.getAllRegisteredActivationMethods();
		for(ActivationMethod m : methods)
		{
			if(doNotAddMergeMode && m == ActivationMethod.FUSION)
				continue;
			fragParam.registerObject(m);
			if(m == defaultMethod)
				fragParam.setDefault();
		}
		addParameter(fragParam);
	}

	public void addInstTypeParam()
	{
		addInstTypeParam(InstrumentType.LOW_RESOLUTION_LTQ);
	}
	
	public void addInstTypeParam(InstrumentType defaultInst)
	{
		ObjectEnumParameter<InstrumentType> instParam = new ObjectEnumParameter<InstrumentType>("inst", "InstrumentID");
		InstrumentType[] allInstTypes = InstrumentType.getAllRegisteredInstrumentTypes();
		for(InstrumentType inst : allInstTypes)
		{
			instParam.registerObject(inst);
			if(inst == defaultInst)
				instParam.setDefault();
		}
		addParameter(instParam);
	}
	
	public void addEnzymeParam()
	{
		addEnzymeParam(Enzyme.TRYPSIN);
	}
	
	public void addEnzymeParam(Enzyme defaulEnzyme)
	{
		ObjectEnumParameter<Enzyme> enzParam = new ObjectEnumParameter<Enzyme>("e", "EnzymeID");
		Enzyme[] allEnzymes = Enzyme.getAllRegisteredEnzymes();
		for(Enzyme e : allEnzymes)
		{
			enzParam.registerObject(e);
			if(e == defaulEnzyme)
				enzParam.setDefault();
		}
		addParameter(enzParam);
	}
	
	public void addProtocolParam()
	{
		addProtocolParam(Protocol.NOPROTOCOL);
	}
	
	public void addProtocolParam(Protocol defaultProtocol)
	{
		ObjectEnumParameter<Protocol> protocolParam = new ObjectEnumParameter<Protocol>("protocol", "ProtocolID");
		Protocol[] protocols = Protocol.getAllRegisteredProtocols();
		for(Protocol protocol : protocols)
		{
			protocolParam.registerObject(protocol);
			if(protocol == defaultProtocol)
				protocolParam.setDefault();
		}
		addParameter(protocolParam);
	}
	
	public void addModFileParam()
	{
		FileParameter modParam = new FileParameter("mod", "ModificationFileName", "Modification file, Default: standard amino acids with fixed C+57");
		modParam.setAsOptional();
		modParam.fileMustExist();
		addParameter(modParam);
	}
	
	public void addMSGFPlusParams()
	{
		addSpecFileParam();
		addDBFileParam();
		addMzIdOutputFileParam();

		ToleranceParameter pmTolParam = new ToleranceParameter("t", "PrecursorMassTolerance", "e.g. 2.5Da, 30ppm or 0.5Da,2.5Da, Default: 20ppm");
		pmTolParam.defaultValue("20ppm");
		pmTolParam.setAdditionalDescription("Use comma to set asymmetric values. E.g. \"-t 0.5Da,2.5Da\" will set 0.5Da to the minus (expMass<theoMass) and 2.5Da to plus (expMass>theoMass)");
		addParameter(pmTolParam);

		IntRangeParameter c13Range = new IntRangeParameter("c13", "Range13CError", "Range of allowed (13C-12C) mass error, Default:-1,2");
		c13Range.setAdditionalDescription("The combination of -t and -13c determins the precursor mass tolerance.\n" +
				"\t   E.g. \"-t 20ppm -nt -1,2\" tests abs(exp-calc-n*1.00335Da)<20ppm for n=-1, 0, 1, 2.");
		c13Range.defaultValue("-1,2");
		addParameter(c13Range);
		
		IntParameter numThreadParam = new IntParameter("thread", "NumThreads", "Number of concurrent threads to be executed, Default: Number of available cores");
		numThreadParam.defaultValue(Runtime.getRuntime().availableProcessors());
		numThreadParam.minValue(1);
		addParameter(numThreadParam);
		
		EnumParameter tdaParam = new EnumParameter("tda");
		tdaParam.registerEntry("don't search decoy database").setDefault();
		tdaParam.registerEntry("search decoy database");
		addParameter(tdaParam);
		
		addFragMethodParam();
		addInstTypeParam();
		addEnzymeParam();
		addProtocolParam();
		
		EnumParameter nttParam = new EnumParameter("ntt", null, "Number of tolerable termini");
		nttParam.setAdditionalDescription("E.g. For trypsin, 0: fully-tryptic peptides only, 1: semi-tryptic, 2: non-tryptic.");
		nttParam.registerEntry("");
		nttParam.registerEntry("").setDefault();
		nttParam.registerEntry("");
		addParameter(nttParam);
		
		addModFileParam();
		
		IntParameter minLenParam = new IntParameter("minLength", "MinPepLength", "Minimum peptide length to consider, Default: 6");
		minLenParam.minValue(1);
		minLenParam.defaultValue(6);
		addParameter(minLenParam);

		IntParameter maxLenParam = new IntParameter("maxLength", "MaxPepLength", "Maximum peptide length to consider, Default: 40");
		maxLenParam.minValue(1);
		maxLenParam.defaultValue(40);
		addParameter(maxLenParam);

		IntParameter minCharge = new IntParameter("minCharge", "MinCharge", "Minimum precursor charge to consider if charges are not specified in the spectrum file, Default: 2");
		minCharge.minValue(1);
		minCharge.defaultValue(2);
		addParameter(minCharge);

		IntParameter maxCharge = new IntParameter("maxCharge", "MaxCharge", "Maximum precursor charge to consider if charges are not specified in the spectrum file, Default: 3");
		maxCharge.minValue(1);
		maxCharge.defaultValue(3);
		addParameter(maxCharge);
		
		IntParameter numMatchesParam = new IntParameter("n", "NumMatchesPerSpec", "Number of matches per spectrum to be reported, Default: 1");
		numMatchesParam.minValue(1);
		numMatchesParam.defaultValue(1);
		addParameter(numMatchesParam);
		
		addExample("Example (high-precision): java -Xmx2000M -jar MSGFPlus.jar -s test.mzXML -d IPI_human_3.79.fasta -t 20ppm -c13 -1,2 -ntt 0 -tda 1 -o testMSGFDB.mzid");
		addExample("Example (low-precision): java -Xmx2000M -jar MSGFPlus.jar -s test.mzXML -d IPI_human_3.79.fasta -t 0.5Da,2.5Da -ntt 0 -tda 1 -o testMSGFDB.mzid");
		
		// Hidden parameters
		FileParameter dbIndexDirParam = new FileParameter("dd", "DBIndexDir", "Path to the directory containing database index files");
		dbIndexDirParam.fileMustExist();
		dbIndexDirParam.mustBeADirectory();
		dbIndexDirParam.setAsOptional();
		dbIndexDirParam.setHidden();
		addParameter(dbIndexDirParam);
		
		EnumParameter unitParam = new EnumParameter("u");
		unitParam.registerEntry("Da");
		unitParam.registerEntry("ppm");
		unitParam.registerEntry("Don't care").setDefault();
		unitParam.setHidden();
		addParameter(unitParam);
		
		IntRangeParameter specIndexParam = new IntRangeParameter("index", "SpecIndex", "Range of spectrum index to be considered");
		specIndexParam.minValue(1);
		specIndexParam.setMaxInclusive();
		specIndexParam.defaultValue("1,"+(Integer.MAX_VALUE-1));
		specIndexParam.setHidden();
		addParameter(specIndexParam);
		
		EnumParameter showFDRParam = new EnumParameter("showFDR");
		showFDRParam.registerEntry("do not show FDRs");
		showFDRParam.registerEntry("show FDRs").setDefault();
		showFDRParam.setHidden();
		addParameter(showFDRParam);

		EnumParameter showDecoyParam = new EnumParameter("showDecoy");
		showDecoyParam.registerEntry("do not show decoy PSMs").setDefault();
		showDecoyParam.registerEntry("show decoy PSMs");
		showDecoyParam.setHidden();
		addParameter(showDecoyParam);
		
		EnumParameter replicateMergedResParam = new EnumParameter("replicate");
		replicateMergedResParam.registerEntry("show merged spectra").setDefault();
		replicateMergedResParam.registerEntry("show individual spectra");
		replicateMergedResParam.setHidden();
		addParameter(replicateMergedResParam);

		EnumParameter edgeScoreParam = new EnumParameter("edgeScore");
		edgeScoreParam.registerEntry("use edge scoring").setDefault();
		edgeScoreParam.registerEntry("do not use edge scoring");
		edgeScoreParam.setHidden();
		addParameter(edgeScoreParam);
	}
	
	public void addMSGFDBParams()
	{
		addSpecFileParam();
		addDBFileParam();
		addPMTolParam();
		addOutputFileParam();
		
		IntParameter numThreadParam = new IntParameter("thread", "NumThreads", "Number of concurrent threads to be executed, Default: Number of available cores");
		numThreadParam.defaultValue(Runtime.getRuntime().availableProcessors());
		numThreadParam.minValue(1);
		addParameter(numThreadParam);
		
		EnumParameter tdaParam = new EnumParameter("tda");
		tdaParam.registerEntry("don't search decoy database").setDefault();
		tdaParam.registerEntry("search decoy database to compute FDR");
		addParameter(tdaParam);
		
		addFragMethodParam();
		addInstTypeParam();
		addEnzymeParam();
		addProtocolParam();
		
		EnumParameter c13Param = new EnumParameter("c13");
		c13Param.registerEntry("Consider only peptides matching precursor mass");
		c13Param.registerEntry("Consider peptides having one 13C").setDefault();
		c13Param.registerEntry("Consider peptides having up to two 13C");
		addParameter(c13Param);
		
		EnumParameter nnetParam = new EnumParameter("nnet", null, "Number of allowed non-enzymatic termini");
		nnetParam.registerEntry("");
		nnetParam.registerEntry("").setDefault();
		nnetParam.registerEntry("");
		addParameter(nnetParam);
		
		addModFileParam();
		
//		FloatRangeParameter itrqParam = new FloatRangeParameter("itraq", "minMass,maxMass", "Remove MS/MS peaks in the mass range between minMass and maxMass (for iTRAQ analysis).");
//		itrqParam.minValue(0f);
//		itrqParam.setMaxInclusive();
//		itrqParam.defaultValue("0,0");
//		itrqParam.setHidden();
//		addParameter(itrqParam);
		
		IntParameter minLenParam = new IntParameter("minLength", "MinPepLength", "Minimum peptide length to consider, Default: 6");
		minLenParam.minValue(1);
		minLenParam.defaultValue(6);
		addParameter(minLenParam);

		IntParameter maxLenParam = new IntParameter("maxLength", "MaxPepLength", "Maximum peptide length to consider, Default: 40");
		maxLenParam.minValue(1);
		maxLenParam.defaultValue(40);
		addParameter(maxLenParam);

		IntParameter minCharge = new IntParameter("minCharge", "MinCharge", "Minimum precursor charge to consider if charges are not specified in the spectrum file, Default: 2");
		minCharge.minValue(1);
		minCharge.defaultValue(2);
		addParameter(minCharge);

		IntParameter maxCharge = new IntParameter("maxCharge", "MaxCharge", "Maximum precursor charge to consider if charges are not specified in the spectrum file, Default: 3");
		maxCharge.minValue(1);
		maxCharge.defaultValue(3);
		addParameter(maxCharge);
		
		IntParameter numMatchesParam = new IntParameter("n", "NumMatchesPerSpec", "Number of matches per spectrum to be reported, Default: 1");
		numMatchesParam.minValue(1);
		numMatchesParam.defaultValue(1);
		addParameter(numMatchesParam);
		
		EnumParameter uniformAAProb = new EnumParameter("uniformAAProb");
		uniformAAProb.registerEntry("use amino acid probabilities computed from the input database").setDefault();
		uniformAAProb.registerEntry("use probability 0.05 for all amino acids");
		addParameter(uniformAAProb);
		
		
		addExample("Example (high-precision): java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 30ppm -c13 1 -nnet 0 -tda 1 -o testMSGFDB.tsv");
		addExample("Example (low-precision): java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 0.5Da,2.5Da -nnet 0 -tda 1 -o testMSGFDB.tsv");
		
		// Hidden parameters
		FileParameter dbIndexDirParam = new FileParameter("dd", "DBIndexDir", "Path to the directory containing database index files");
		dbIndexDirParam.fileMustExist();
		dbIndexDirParam.mustBeADirectory();
		dbIndexDirParam.setAsOptional();
		dbIndexDirParam.setHidden();
		addParameter(dbIndexDirParam);
		
		EnumParameter unitParam = new EnumParameter("u");
		unitParam.registerEntry("Da");
		unitParam.registerEntry("ppm");
		unitParam.registerEntry("Don't care").setDefault();
		unitParam.setHidden();
		addParameter(unitParam);
		
		IntRangeParameter specIndexParam = new IntRangeParameter("index", "SpecIndex", "Range of spectrum index to be considered");
		specIndexParam.minValue(1);
		specIndexParam.setMaxInclusive();
		specIndexParam.defaultValue("1,"+(Integer.MAX_VALUE-1));
		specIndexParam.setHidden();
		addParameter(specIndexParam);
		
		EnumParameter showFDRParam = new EnumParameter("showFDR");
		showFDRParam.registerEntry("do not show FDRs");
		showFDRParam.registerEntry("show FDRs").setDefault();
		showFDRParam.setHidden();
		addParameter(showFDRParam);

		EnumParameter showDecoyParam = new EnumParameter("showDecoy");
		showDecoyParam.registerEntry("do not show decoy PSMs").setDefault();
		showDecoyParam.registerEntry("show decoy PSMs");
		showDecoyParam.setHidden();
		addParameter(showDecoyParam);
		
		EnumParameter replicateMergedResParam = new EnumParameter("replicate");
		replicateMergedResParam.registerEntry("show merged spectra").setDefault();
		replicateMergedResParam.registerEntry("show individual spectra");
		replicateMergedResParam.setHidden();
		addParameter(replicateMergedResParam);

		EnumParameter edgeScoreParam = new EnumParameter("edgeScore");
		edgeScoreParam.registerEntry("use edge scoring").setDefault();
		edgeScoreParam.registerEntry("do not use edge scoring");
		edgeScoreParam.setHidden();
		addParameter(edgeScoreParam);
		
		EnumParameter precolatorParam = new EnumParameter("percolator");
		edgeScoreParam.registerEntry("normal").setDefault();
		edgeScoreParam.registerEntry("for MS-GF+Percolator");
		edgeScoreParam.setHidden();
		addParameter(precolatorParam);		
		
	}

	public void addMSGFParams()
	{
		// SpectrumFile
		FileParameter resFileParam = new FileParameter("i", "ResultFile", "ResultFile");
		resFileParam.fileMustExist();
		addParameter(resFileParam);
		
		// SpecDir
		FileParameter specDirParam = new FileParameter("d", "SpecDir", "Path to directory containing spectrum files");
		specDirParam.mustBeADirectory();
		specDirParam.fileMustExist();
		addParameter(specDirParam);

		// OutputFileName
		addOutputFileParam();
		
		// DBFile
		addDBFileParam("db", "To get AA frequencies, if not specified, 1/20 is used for all AAs", true);

		// Fragmentation method
		addFragMethodParam(ActivationMethod.ASWRITTEN, true);
		
		// Instrument type
		addInstTypeParam();
		
		// Enzyme
		addEnzymeParam();
		
		// FixedMod
		EnumParameter fixModParam = new EnumParameter("fixMod");
		fixModParam.registerEntry("NoCysteineProtection");
		fixModParam.registerEntry("Carbamidomethyl-C").setDefault();
		fixModParam.registerEntry("Carboxymethyl-C");
		addParameter(fixModParam);
		
		// -x
		EnumParameter numSpecParam = new EnumParameter("x");
		numSpecParam.registerEntry("All").setDefault();
		numSpecParam.registerEntry("OnePerSpec");
		addParameter(numSpecParam);
		
		// -p
		FloatParameter spThParam = new FloatParameter("p", "SpecProbThreshold", "Spectral probability threshold (Default: 1)");
		spThParam.minValue(0f).setMinExclusive();
		spThParam.maxValue(1f).setMaxInclusive();
		spThParam.defaultValue(1f);
		addParameter(spThParam);
		
		// -addScore
		EnumParameter addScoreParam = new EnumParameter("addScore");
		addScoreParam.registerEntry("Don't add MSGFScore").setDefault();
		addScoreParam.registerEntry("Add MSGFScore");
		addParameter(addScoreParam);
	}
	
	public void addMSGFLibParams()
	{
		addSpecFileParam();
		
		// Add library file param
		FileParameter libFileParam = new FileParameter("d", "LibraryFile", "*.sptxt");
		libFileParam.addFileFormat(new FileFormat(".sptxt"));
		libFileParam.fileMustExist();
		libFileParam.mustBeAFile();
		addParameter(libFileParam);

		addPMTolParam();
		addOutputFileParam();
		
		IntParameter numThreadParam = new IntParameter("thread", "NumThreads", "Number of concurrent threads to be executed, Default: Number of available cores");
		numThreadParam.defaultValue(Runtime.getRuntime().availableProcessors());
		numThreadParam.minValue(1);
		addParameter(numThreadParam);
		
		addFragMethodParam();
		addInstTypeParam();
		addEnzymeParam();
		addProtocolParam();
		
		EnumParameter c13Param = new EnumParameter("c13");
		c13Param.registerEntry("Consider only peptides matching precursor mass");
		c13Param.registerEntry("Consider peptides having one 13C").setDefault();
		c13Param.registerEntry("Consider peptides having up to two 13C");
		addParameter(c13Param);
		
		IntParameter numMatchesParam = new IntParameter("n", "NumMatchesPerSpec", "Number of matches per spectrum to be reported, Default: 1");
		numMatchesParam.minValue(1);
		numMatchesParam.defaultValue(1);
		addParameter(numMatchesParam);
		
		addExample("Example: java -Xmx2000M -jar MSGFLib.jar -s test.mzXML -d IPI_human_3.79.fasta -t 30ppm -c13 1 -nnet 0 -o testMSGFDB.tsv");
	}
	
	public FileParameter getSpecFileParam()
	{
		return ((FileParameter)getParameter("s"));
	}

	public FileParameter getDBFileParam()
	{
		return ((FileParameter)getParameter("d"));
	}
	
	public ToleranceParameter getPMTolParam()
	{
		return ((ToleranceParameter)getParameter("t"));
	}
	
	public FileParameter getOutputFileParam()
	{
		return ((FileParameter)getParameter("o"));
	}
	
	public ActivationMethod getActivationMethod()
	{
		return (ActivationMethod)((ObjectEnumParameter<?>)getParameter("m")).getObject();
	}

	public InstrumentType getInstType()
	{
		return (InstrumentType)((ObjectEnumParameter<?>)getParameter("inst")).getObject();
	}

	public Enzyme getEnzyme()
	{
		return (Enzyme)((ObjectEnumParameter<?>)getParameter("e")).getObject();
	}

	public Protocol getProtocol()
	{
		return (Protocol)((ObjectEnumParameter<?>)getParameter("protocol")).getObject();
	}
	
	public FileParameter getModFileParam()
	{
		return ((FileParameter)getParameter("mod"));
	}
	
	public int getIntValue(String key)
	{
		Parameter param = this.getParameter(key);
		if(param instanceof IntParameter)
			return ((IntParameter)param).getValue();
		else
		{
			System.err.println("[Error] in ParamManager.getIntValue: " + key + " is not an instance of IntParameter.");
			System.exit(-1);
		}
		return -1;
	}

	public float getFloatValue(String key)
	{
		Parameter param = this.getParameter(key);
		if(param instanceof FloatParameter)
			return ((FloatParameter)param).getValue();
		else
		{
			System.err.println("[Error] in ParamManager.getFloatValue: " + key + " is not an instance of FloatParameter.");
			System.exit(-1);
		}
		return -1;
	}
	
	public File getFile(String key)
	{
		Parameter param = this.getParameter(key);
		if(param instanceof FileParameter)
			return ((FileParameter)param).getFile();
		else
		{
			System.err.println("[Error] in ParamManager.getFile: " + key + " is not an instance of FileParameter.");
			System.exit(-1);
		}
		return null;
	}

	public File[] getFiles(String key)
	{
		Parameter param = this.getParameter(key);
		if(param instanceof FileListParameter)
			return ((FileListParameter)param).getFiles();
		else
		{
			System.err.println("[Error] in ParamManager.getFile: " + key + " is not an instance of FileListParameter.");
			System.exit(-1);
		}
		return null;
	}
	
	public static void main(String argv[])
	{
		ParamManager paramManager = new ParamManager("MSGF", MSGF.VERSION, MSGF.RELEASE_DATE, "java -Xmx2000M -jar MSGFDB.jar");
		paramManager.addMSGFDBParams();
		
//		FileListParameter testParam = new FileListParameter("test", "test", "test");
//		testParam.addFileFormat(SpecFileFormat.MGF);
//		testParam.addFileFormat(SpecFileFormat.MZXML);
//		paramManager.addParameter(testParam);
		
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
