package edu.ucsd.msjava.params;

import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.sequences.Constants;
import edu.ucsd.msjava.ui.MSGF;
import edu.ucsd.msjava.ui.MSGFPlus;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map.Entry;


public class ParamManager {
    private CaseInsensitiveLinkedHashMapParam params;
    private String toolName;
    private String version;
    private String date;
    private String command;
    private ArrayList<String> examples = new ArrayList<String>();


    public enum ParamNameEnum{

        SPECTRUM_FILE("s", "SpectrumFile", "*.mzML, *.mzXML, *.mgf, *.ms2, *.pkl or *_dta.txt",
                "Spectra should be centroided (see below for MSConvert example). Profile spectra will be ignored."),

        DB_FILE("d", "DatabaseFile", "*.fasta or *.fa or *.faa", null),

        DECOY_PREFIX("decoy", "DecoyPrefix",
                "Prefix for decoy protein names; default is " + MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX, null),

        PARENT_MASS_TOLERANCE("t", "ParentMassTolerance", "e.g. 2.5Da, 30ppm or 0.5Da,2.5Da",
                "Use a comma to set asymmetric values. E.g. \"-t 0.5Da,2.5Da\" will set 0.5Da to the left (expMass<theoMass) and 2.5Da to the right (expMass>theoMass)"),

        MZID_OUTPUT_FILE("o", "OutputFile (*.mzid)", "Default: [SpectrumFileName].mzid", null),
        OUTPUT_FILE("o", "OutputFile", "Default: stdout", null),
        FRAG_METHOD("m", "FragmentMethodID", null, null),
        INTRUMENT_TYPE("inst", "InstrumentID", null, null),
        ENZYME_ID("e", "EnzymeID", null, null),
        PROTOCOL_ID("protocol", "ProtocolID", null, null),
        MOD_FILE("mod", "ModificationFileName", "Modification file, Default: standard amino acids with fixed C+57; only if -mod is not specified", null),
        CONFIGURATION_FILE("conf", "ConfigurationFileName", "Configuration file, Default: all parameters are provided by the commandline interface", null),

        NUM_THREADS("thread", "NumThreads", "Number of concurrent threads to be executed, Default: Number of available cores", "This is best set to the number of physical cores in a single NUMA node.\n" +
                "\t   Generally a single NUMA node is 1 physical processor.\n" +
                "\t   The default will try to use hyperthreading cores, which can increase the amount of time this process will take.\n" +
                "\t   This is because the part of Scoring param generation that is multithreaded is also I/O intensive."),

        NUM_TASKS("tasks", "NumTasks", "Override the number of tasks to use on the threads, Default: (internally calculated based on inputs)", "More tasks than threads will reduce the memory requirements of the search, but will be slower (how much depends on the inputs).\n" +
                "\t   1 <= tasks <= numThreads: will create one task per thread, which is the original behavior.\n" +
                "\t   tasks = 0: use default calculation - minimum of: (threads*3) and (numSpectra/250).\n" +
                "\t   tasks < 0: multiply number of threads by abs(tasks) to determine number of tasks (i.e., -2 means \"2 * numThreads\" tasks).\n" +
                "\t   One task per thread will use the most memory, but will usually finish the fastest.\n" +
                "\t   2-3 tasks per thread will use comparably less memory, but may cause the search to take 1.5 to 2 times as long."),

        PRECURSOR_MASS_TOLERANCE("t", "PrecursorMassTolerance", "e.g. 2.5Da, 20ppm or 0.5Da,2.5Da, Default: 20ppm", "Use a comma to set asymmetric values. E.g. \"-t 0.5Da,2.5Da\" will set 0.5Da to the left (ObsMass < TheoMass) and 2.5Da to the right (ObsMass > TheoMass)"),

        ISOTOPE_ERROR("ti", "IsotopeErrorRange", "Range of allowed isotope peak errors, Default:0,1", "Takes into account the error introduced by choosing a non-monoisotopic peak for fragmentation.\n" +
                "\t   The combination of -t and -ti determines the precursor mass tolerance.\n" +
                "\t   E.g. \"-t 20ppm -ti -1,2\" tests abs(ObservedPepMass - TheoreticalPepMass - n * 1.00335Da) < 20ppm for n = -1, 0, 1, 2."),

        ENZYME_SPECIFICITY("ntt", "NTT", "Number of Tolerable Termini", "E.g. For trypsin, 0: non-tryptic, 1: semi-tryptic, 2: fully-tryptic peptides only."),

        MIN_PEPTIDE_LENGTH("minLength", "MinPepLength", "Minimum peptide length to consider, Default: 6", null),
        MAX_PEPTIDE_LENGTH("maxLength", "MaxPepLength", "Maximum peptide length to consider, Default: 40", null),
        MIN_CHARGE("minCharge", "MinCharge", "Minimum precursor charge to consider if charges are not specified in the spectrum file, Default: 2", null),
        MAX_CHARGE("maxCharge", "MaxCharge", "Maximum precursor charge to consider if charges are not specified in the spectrum file, Default: 3", null),
        NUM_MATCHES_SPEC("n", "NumMatchesPerSpec", "Number of matches per spectrum to be reported, Default: 1", null),
        CHARGE_CARRIER_MASSES("ccm", "ChargeCarrierMass", "Mass of charge carrier, Default: mass of proton (1.00727649)", null),
        MIN_NUM_PEAKS("minNumPeaks", "MinNumPeaksPerSpectrum", "Minimum number of peaks per spectrum, Default: " + Constants.MIN_NUM_PEAKS_PER_SPECTRUM, null),
        NUM_ISOFORMS("iso", "NumIsoforms", "Number of isoforms to consider per peptide, Default: 128" + Constants.NUM_VARIANTS_PER_PEPTIDE, null),

        MIN_DE_NOVO_SCORE("minDeNovoScore", "MinDeNovoScore", "Minimum de Novo score, Default: " + Constants.MIN_DE_NOVO_SCORE, null),
        SPEC_INDEX("index", "SpecIndex", "Range of spectrum index to be considered", null),
        MAX_MISSED_CLEAVAGES("maxMissedCleavages", "Count", "Exclude peptides with more than this number of missed cleavages from the search, Default: -1 (no limit)", null),
        TDA_STRATEGY("tda", "TDA", "Target decoy strategy", null),
        ADD_FEATURES("addFeatures", "AddFeatures", "Add features in the output", null),
        DD_DIRECTORY("dd", "DBIndexDir", "Path to the directory containing database index files", null),
        UNIFORM_AA_PROBABILITY("uniformAAProb", "UniformAAProb", null, null),
        MAX_NUM_MODS("numMods", "NumMods", "Maximum number of modifications", null),
        STATIC_MODIFICATION("staticMod", "StaticMod", "Static/Fixed modification", null),
        DYNAMIC_MODIFICATION("dynamicMod", "DynamicMod", "Dynamic/Variable modification", null),
        CUSTOM_AA("customAA", "CustomAA", "Custom amino acid", null),


        VERBOSE("verbose", null, null,null);

        private String name;
        private String commandlineName;
        private String description;
        private String additionalDescription;

        ParamNameEnum(String commandlineName, String name, String description, String additionalDescription){
            this.name = name;
            this.commandlineName = commandlineName;
            this.description = description;
            this.additionalDescription = additionalDescription;
        }

        public String getName() {
            return name;
        }

        public String getCommandlineName() {
            return commandlineName;
        }

        public String getDescription() {
            return description;
        }

        public String getAdditionalDescription() {
            return additionalDescription;
        }

        /**
         * Check is the parameter line contains the ParamValue
         * @param line Param Line
         * @return is the line.
         */
        public boolean isLine(String line) {
            return ((getName()!= null && line.contains(getName().toLowerCase())));
        }
    }

    public ParamManager(String toolName, String version, String date, String command) {
        this.toolName = toolName;
        this.version = version;
        this.date = date;
        this.command = command;
        params = new CaseInsensitiveLinkedHashMapParam();
    }

    public boolean addParameter(Parameter param) {
        if (params.containsKey(param.getKey())) {
            System.err.println("ParamManager: duplicate key (" + param.getKey() + ")");
            System.exit(-1);
        }
        params.put(param.getKey(), param);
        return true;
    }

    public void addExample(String example) {
        this.examples.add(example);
    }

    public Parameter getParameter(String key) {
        return params.get(key);
    }

    /**
     * Validates the parameters
     * @return Error message if an error, otherwise null
     */
    public String isValid() {
        Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
        while (itr.hasNext()) {
            Entry<String, Parameter> entry = itr.next();
            Parameter param = entry.getValue();
            if (!param.isValid()) {
                return "Parameter -" + param.getKey() + " (" + param.getName() + ") is missing";
            }
        }
        return null;
    }

    public void printToolInfo() {
        System.out.println(this.toolName + " " + this.version + " (" + this.date + ")");
    }

    public void printUsageInfo() {
        System.out.println();
        System.out.println(this.toolName + " " + this.version + " (" + this.date + ")");
        System.out.println();
        System.out.println("Usage: " + this.command);

        ArrayList<Parameter> optParams = new ArrayList<Parameter>();
        Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
        while (itr.hasNext()) {
            Entry<String, Parameter> entry = itr.next();
            Parameter param = entry.getValue();
            if (!param.isHidden()) {
                if (!param.isOptional()) {
                    System.out.println("\t" + param);
                    if (param.getAdditionalDescription() != null)
                        System.out.println("\t   " + param.getAdditionalDescription());
                } else
                    optParams.add(param);
            }
        }

        for (Parameter param : optParams) {
            System.out.println("\t" + param);
            if (param.getAdditionalDescription() != null)
                System.out.println("\t   " + param.getAdditionalDescription());
        }

        System.out.println();
        for (String example : examples)
            System.out.println(example);

        System.out.println();
        System.out.println("For Thermo .raw files, obtain a centroided .mzML file using MSConvert, which is part of ProteoWizard (http://proteowizard.sourceforge.net/)");
        System.out.println("  MSConvert.exe DatasetName.raw --filter \"peakPicking true 1-\" --mzML --32");
        System.out.println();
        System.out.println("To add or override the enzyme definitions, create a file named enzymes.txt in a directory named params below the working directory.");
        System.out.println("For example, create file C:\\Work\\params\\enzymes.txt when the working directory is C:\\Work");
        System.out.println("Example enzymes.txt file: https://github.com/MSGFPlus/msgfplus/blob/master/docs/examples/enzymes.txt");
        System.out.println();
        System.out.println("Documentation: https://msgfplus.github.io/msgfplus/");
        System.out.println("Releases:      https://github.com/MSGFPlus/msgfplus/releases");
    }

    public void printValues() {
        Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
        while (itr.hasNext()) {
            Entry<String, Parameter> entry = itr.next();
            Parameter param = entry.getValue();
            System.out.println(param.getKey() + "\t" + param.getValueAsString());
        }
    }

    public String parseParams(String argv[]) {
        if (argv.length == 0) {
            return "No parameter specified.";
        }

        if (argv.length < 2 || argv.length % 2 != 0) {
            return "The number of parameters must be even. If a file path has a space, surround it with double quotes.";
        }

        for (int i = 0; i < argv.length; i += 2) {
            if (!argv[i].startsWith("-") || i + 1 >= argv.length || argv[i].length() <= 1) {
                return "Syntax error; parameter names must start with a dash: " + argv[i];
            } else {
                String key = argv[i].substring(1);
                Parameter param = params.get(key);
                if (param == null) {
                    return "Invalid parameter: " + argv[i] + ".";
                } else {
                    String error = param.parse(argv[i + 1]);
                    if (error != null) {
                        String err = "Invalid value for parameter " + argv[i] + ": " + argv[i + 1];
                        err += "\n        (" + error + ")";
                        return err;
                    }
                    param.setValueAssigned();
                }
            }
        }

        String error = isValid();
        if (error != null)
            return error;

        return null;
    }

    public void addSpecFileParam() {
        FileParameter specFileParam = new FileParameter(ParamNameEnum.SPECTRUM_FILE.commandlineName,
                ParamNameEnum.SPECTRUM_FILE.name, ParamNameEnum.SPECTRUM_FILE.description);
        specFileParam.addFileFormat(SpecFileFormat.MZML);
        specFileParam.addFileFormat(SpecFileFormat.MZXML);
        specFileParam.addFileFormat(SpecFileFormat.MGF);
        specFileParam.addFileFormat(SpecFileFormat.MS2);
        specFileParam.addFileFormat(SpecFileFormat.PKL);
        specFileParam.addFileFormat(SpecFileFormat.DTA_TXT);
        specFileParam.addFileFormat(FileFormat.DIRECTORY);
        specFileParam.fileMustExist();
        specFileParam.setAdditionalDescription(ParamNameEnum.SPECTRUM_FILE.additionalDescription);
        addParameter(specFileParam);
    }

    public void addDBFileParam() {
        addDBFileParam(ParamNameEnum.DB_FILE.commandlineName, ParamNameEnum.DB_FILE.description, false);
    }

    public void addDBFileParam(String key, String description, boolean isOptional) {
        FileParameter dbFileParam = new FileParameter(key, ParamNameEnum.DB_FILE.name, description);
        if (isOptional)
            dbFileParam.setAsOptional();
        dbFileParam.addFileFormat(DBFileFormat.FASTA);
        dbFileParam.fileMustExist();
        dbFileParam.mustBeAFile();
        addParameter(dbFileParam);
    }

    public void addDecoyPrefixParam() {
        addDecoyPrefixParam(MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX);
    }

    public void addDecoyPrefixParam(String defaultDecoyPrefix) {
        StringParameter decoyPrefixParam = new StringParameter(ParamNameEnum.DECOY_PREFIX.commandlineName, ParamNameEnum.DECOY_PREFIX.name, ParamNameEnum.DECOY_PREFIX.description);
        // Note that defining a default value auto-sets isOptional to True
        decoyPrefixParam.defaultValue(defaultDecoyPrefix);
        addParameter(decoyPrefixParam);
    }

    public void addPMTolParam() {
        ToleranceParameter pmTolParam = new ToleranceParameter(ParamNameEnum.PARENT_MASS_TOLERANCE.commandlineName, ParamNameEnum.PARENT_MASS_TOLERANCE.name, ParamNameEnum.PARENT_MASS_TOLERANCE.description);
        pmTolParam.setAdditionalDescription(ParamNameEnum.PARENT_MASS_TOLERANCE.additionalDescription);
        addParameter(pmTolParam);
    }

    public void addMzIdOutputFileParam() {
        FileParameter outputParam = new FileParameter(ParamNameEnum.MZID_OUTPUT_FILE.commandlineName, ParamNameEnum.MZID_OUTPUT_FILE.name, ParamNameEnum.MZID_OUTPUT_FILE.description);
        outputParam.addFileFormat(new FileFormat(".mzid").setCaseSensitive());
        outputParam.setAsOptional();
        addParameter(outputParam);
    }

    public void addOutputFileParam() {
        FileParameter outputParam = new FileParameter(ParamNameEnum.OUTPUT_FILE.commandlineName, ParamNameEnum.OUTPUT_FILE.name, ParamNameEnum.OUTPUT_FILE.description);
        outputParam.setAsOptional();
        outputParam.fileMustNotExist();
        addParameter(outputParam);
    }

    public void addFragMethodParam() {
        addFragMethodParam(ActivationMethod.ASWRITTEN, false);
    }

    /**
     * Used by both MS-GFDB and MS-GF+
     * MS-GF+ passes True for doNotAddMergeMode, thus ignoring ActivationMethod.FUSION
     *
     * @param defaultMethod
     * @param doNotAddMergeMode
     */
    public void addFragMethodParam(ActivationMethod defaultMethod, boolean doNotAddMergeMode) {
        ObjectEnumParameter<ActivationMethod> fragParam = new ObjectEnumParameter<ActivationMethod>(ParamNameEnum.FRAG_METHOD.commandlineName, ParamNameEnum.FRAG_METHOD.name);
        ActivationMethod[] methods = ActivationMethod.getAllRegisteredActivationMethods();
        for (ActivationMethod m : methods) {
            if (doNotAddMergeMode && m == ActivationMethod.FUSION)
                continue;
            fragParam.registerObject(m);
            if (m == defaultMethod)
                fragParam.setDefault();
        }
        addParameter(fragParam);
    }

    public void addInstTypeParam() {
        addInstTypeParam(InstrumentType.LOW_RESOLUTION_LTQ);
    }

    public void addInstTypeParam(InstrumentType defaultInst) {
        ObjectEnumParameter<InstrumentType> instParam = new ObjectEnumParameter<InstrumentType>(ParamNameEnum.INTRUMENT_TYPE.commandlineName, ParamNameEnum.INTRUMENT_TYPE.name);
        InstrumentType[] allInstTypes = InstrumentType.getAllRegisteredInstrumentTypes();
        for (InstrumentType inst : allInstTypes) {
            instParam.registerObject(inst);
            if (inst == defaultInst)
                instParam.setDefault();
        }
        addParameter(instParam);
    }

    public void addEnzymeParam() {
        addEnzymeParam(Enzyme.TRYPSIN);
    }

    public void addEnzymeParam(Enzyme enzymeId) {
        ObjectEnumParameter<Enzyme> enzParam = new ObjectEnumParameter<Enzyme>(ParamNameEnum.ENZYME_ID.commandlineName, ParamNameEnum.ENZYME_ID.name);
        Enzyme[] allEnzymes = Enzyme.getAllRegisteredEnzymes();
        for (Enzyme e : allEnzymes) {
            enzParam.registerObject(e);
            if (e == enzymeId)
                enzParam.setDefault();
        }
        addParameter(enzParam);
    }

    public void addProtocolParam() {
        addProtocolParam(Protocol.AUTOMATIC);
    }

    public void addProtocolParam(Protocol defaultProtocol) {
        ObjectEnumParameter<Protocol> protocolParam = new ObjectEnumParameter<Protocol>(ParamNameEnum.PROTOCOL_ID.commandlineName, ParamNameEnum.PROTOCOL_ID.name);
        Protocol[] protocols = Protocol.getAllRegisteredProtocols();
        for (Protocol protocol : protocols) {
            protocolParam.registerObject(protocol);
            if (protocol == defaultProtocol)
                protocolParam.setDefault();
        }
        addParameter(protocolParam);
    }

    public void addModFileParam() {
        FileParameter modParam = new FileParameter(ParamNameEnum.MOD_FILE);
        modParam.setAsOptional();
        modParam.fileMustExist();
        addParameter(modParam);
    }

    public void addConfigFileParam(){
        FileParameter configFile = new FileParameter(ParamNameEnum.CONFIGURATION_FILE);
        configFile.setAsOptional();
        configFile.fileMustExist();
        addParameter(configFile);
    }

    /**
     * Add parameters for MS-GF+
     */
    public void addMSGFPlusParams() {
        // -s SpectrumFile (*.mzML, *.mzXML, *.mgf, *.ms2, *.pkl or *_dta.txt)
        addSpecFileParam();

        // -d DatabaseFile (*.fasta or *.fa or *.faa)
        addDBFileParam();

        // decoy DecoyPrefix
        addDecoyPrefixParam();

        // [-o OutputFile (*.mzid)] (Default: [SpectrumFileName].mzid)
        addMzIdOutputFileParam();

        ToleranceParameter pmTolParam = new ToleranceParameter(ParamNameEnum.PRECURSOR_MASS_TOLERANCE);
        pmTolParam.defaultValue("20ppm");
        pmTolParam.setAdditionalDescription(ParamNameEnum.PRECURSOR_MASS_TOLERANCE.additionalDescription);
        addParameter(pmTolParam);

        IntRangeParameter isotopeRange = new IntRangeParameter(ParamNameEnum.ISOTOPE_ERROR);
        isotopeRange.setAdditionalDescription(ParamNameEnum.ISOTOPE_ERROR.additionalDescription);
        isotopeRange.setMaxInclusive();
        isotopeRange.defaultValue("0,1");
        addParameter(isotopeRange);

        IntParameter numThreadParam = new IntParameter(ParamNameEnum.NUM_THREADS);
        numThreadParam.defaultValue(Runtime.getRuntime().availableProcessors());
        numThreadParam.minValue(1);
        addParameter(numThreadParam);

        IntParameter numTasksParam = new IntParameter(ParamNameEnum.NUM_TASKS);
        numTasksParam.setAdditionalDescription(ParamNameEnum.NUM_TASKS.additionalDescription);
        numTasksParam.defaultValue(0);
        numTasksParam.minValue(-10);
        addParameter(numTasksParam);

        EnumParameter verboseOutputParam = new EnumParameter(ParamNameEnum.VERBOSE.commandlineName);
        verboseOutputParam.registerEntry("Report total progress only").setDefault();
        verboseOutputParam.registerEntry("Report total and per-thread progress/status");
        addParameter(verboseOutputParam);

        EnumParameter tdaParam = new EnumParameter(ParamNameEnum.TDA_STRATEGY);
        tdaParam.registerEntry("Don't search decoy database").setDefault();
        tdaParam.registerEntry("Search decoy database");
        addParameter(tdaParam);

        addFragMethodParam(ActivationMethod.ASWRITTEN, true);
        addInstTypeParam();
        addEnzymeParam();
        addProtocolParam();

        EnumParameter nttParam = new EnumParameter(ParamNameEnum.ENZYME_SPECIFICITY.commandlineName, null, ParamNameEnum.ENZYME_SPECIFICITY.description);
        nttParam.setAdditionalDescription(ParamNameEnum.ENZYME_SPECIFICITY.additionalDescription);
        nttParam.registerEntry("");
        nttParam.registerEntry("");
        nttParam.registerEntry("").setDefault();
        addParameter(nttParam);

        addModFileParam();

        addConfigFileParam();

        IntParameter minLenParam = new IntParameter(ParamNameEnum.MIN_PEPTIDE_LENGTH);
        minLenParam.minValue(1);
        minLenParam.defaultValue(6);
        addParameter(minLenParam);

        IntParameter maxLenParam = new IntParameter(ParamNameEnum.MAX_PEPTIDE_LENGTH);
        maxLenParam.minValue(1);
        maxLenParam.defaultValue(40);
        addParameter(maxLenParam);

        IntParameter minCharge = new IntParameter(ParamNameEnum.MIN_CHARGE);
        minCharge.minValue(1);
        minCharge.defaultValue(2);
        addParameter(minCharge);

        IntParameter maxCharge = new IntParameter(ParamNameEnum.MAX_CHARGE);
        maxCharge.minValue(1);
        maxCharge.defaultValue(3);
        addParameter(maxCharge);

        IntParameter numMatchesParam = new IntParameter(ParamNameEnum.NUM_MATCHES_SPEC);
        numMatchesParam.minValue(1);
        numMatchesParam.defaultValue(1);
        addParameter(numMatchesParam);

        EnumParameter addFeatureParam = new EnumParameter(ParamNameEnum.ADD_FEATURES);
        addFeatureParam.registerEntry("Output basic scores only").setDefault();
        addFeatureParam.registerEntry("Output additional features");
        addParameter(addFeatureParam);

        DoubleParameter chargeCarrierMassParam = new DoubleParameter(ParamNameEnum.CHARGE_CARRIER_MASSES);
        chargeCarrierMassParam.minValue(0.1);
        chargeCarrierMassParam.setMaxInclusive();
        chargeCarrierMassParam.defaultValue(Composition.PROTON);
        addParameter(chargeCarrierMassParam);

        /* Maximum number of missed cleavages to allow on searched peptides */
        IntParameter maxMissedCleavages = new IntParameter(ParamNameEnum.MAX_MISSED_CLEAVAGES);
        maxMissedCleavages.minValue(-1);
        maxMissedCleavages.defaultValue(-1);
        addParameter(maxMissedCleavages);

        addExample("Example (high-precision): java -Xmx3500M -jar MSGFPlus.jar -s test.mzML -d IPI_human_3.79.fasta -inst 1 -t 20ppm -ti -1,2 -ntt 2 -tda 1 -o testMSGFPlus.mzid -mod Mods.txt");
        addExample("Example (low-precision):  java -Xmx3500M -jar MSGFPlus.jar -s test.mzML -d IPI_human_3.79.fasta -inst 0 -t 0.5Da,2.5Da    -ntt 2 -tda 1 -o testMSGFPlus.mzid -mod Mods.txt");

        // Hidden parameters
        FileParameter dbIndexDirParam = new FileParameter(ParamNameEnum.DD_DIRECTORY);
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

        IntRangeParameter specIndexParam = new IntRangeParameter(ParamNameEnum.SPEC_INDEX);
        specIndexParam.minValue(1);
        specIndexParam.setMaxInclusive();
        specIndexParam.defaultValue("1," + (Integer.MAX_VALUE - 1));
        specIndexParam.setHidden();
        addParameter(specIndexParam);

//		EnumParameter showFDRParam = new EnumParameter("showQValue");
//		showFDRParam.registerEntry("do not show Q-values");
//		showFDRParam.registerEntry("show Q-values").setDefault();
//		showFDRParam.setHidden();
//		addParameter(showFDRParam);

//		EnumParameter showDecoyParam = new EnumParameter("showDecoy");
//		showDecoyParam.registerEntry("do not show decoy PSMs").setDefault();
//		showDecoyParam.registerEntry("show decoy PSMs");
//		addParameter(showDecoyParam);

        EnumParameter edgeScoreParam = new EnumParameter("edgeScore");
        edgeScoreParam.registerEntry("Use edge scoring").setDefault();
        edgeScoreParam.registerEntry("Do not use edge scoring");
        edgeScoreParam.setHidden();
        addParameter(edgeScoreParam);

        IntParameter minNumPeaksParam = new IntParameter(ParamNameEnum.MIN_NUM_PEAKS);
        minNumPeaksParam.defaultValue(Constants.MIN_NUM_PEAKS_PER_SPECTRUM);
        minNumPeaksParam.setHidden();
        addParameter(minNumPeaksParam);

        IntParameter isoParam = new IntParameter(ParamNameEnum.NUM_ISOFORMS);
        isoParam.defaultValue(Constants.NUM_VARIANTS_PER_PEPTIDE);
        isoParam.setHidden();
        addParameter(isoParam);

        EnumParameter metCleavageParam = new EnumParameter("ignoreMetCleavage");
        metCleavageParam.registerEntry("Consider protein N-term Met cleavage").setDefault();
        metCleavageParam.registerEntry("Ignore protein N-term Met cleavage");
        metCleavageParam.setHidden();
        addParameter(metCleavageParam);

        IntParameter minDeNovoScoreParam = new IntParameter(ParamNameEnum.MIN_DE_NOVO_SCORE);
        minDeNovoScoreParam.minValue(Integer.MIN_VALUE);
        minDeNovoScoreParam.defaultValue(Constants.MIN_DE_NOVO_SCORE);
        minDeNovoScoreParam.setHidden();
        addParameter(minDeNovoScoreParam);

    }

    /**
     * Used by MS-GF+
     */
    public void addScoringParamGenParams() {
        FileListParameter resFileParam = new FileListParameter("i", "ResultPath", "MSGFDBResultFile (*.mzid) or MSGFDBResultDir");
        resFileParam.addFileFormat(new FileFormat(".mzid"));
        resFileParam.addFileFormat(new FileFormat(".tsv"));
        resFileParam.setAdditionalDescription("mzid files are converted to tsv using default settings before use.\n" +
                "\t   If you are going to run ScoringParamGen multiple times on the same data (with different parameters),\n" +
                "\t   convert any mzid files to tsv prior to running ScoringParamGen.");
        addParameter(resFileParam);

        FileParameter specDirParam = new FileParameter("d", "SpecDir", "Path to directory containing spectrum files");
        specDirParam.mustBeADirectory();
        specDirParam.fileMustExist();
        addParameter(specDirParam);

        // ActivationMethod
        ObjectEnumParameter<ActivationMethod> fragParam = new ObjectEnumParameter<ActivationMethod>("m", "FragmentMethodID");
        ActivationMethod[] methods = ActivationMethod.getAllRegisteredActivationMethods();
        for (int i = 1; i < methods.length; i++) {
            ActivationMethod m = methods[i];
            if (m != ActivationMethod.FUSION)
                fragParam.registerObject(m);
        }
        addParameter(fragParam);

        // Instrument type
        addInstTypeParam(null);

        // Enzyme
        ObjectEnumParameter<Enzyme> enzParam = new ObjectEnumParameter<Enzyme>("e", "EnzymeID");
        Enzyme[] allEnzymes = Enzyme.getAllRegisteredEnzymes();
        for (int i = 1; i < allEnzymes.length; i++) {
            Enzyme e = allEnzymes[i];
            enzParam.registerObject(e);
        }
        addParameter(enzParam);

        // Protocol
        addProtocolParam();

        IntParameter numThreadParam = new IntParameter(ParamNameEnum.NUM_THREADS);
        numThreadParam.defaultValue(Runtime.getRuntime().availableProcessors() / 2);
        numThreadParam.setAdditionalDescription(ParamNameEnum.NUM_THREADS.additionalDescription);
        numThreadParam.minValue(1);
        addParameter(numThreadParam);

        EnumParameter dropErrors = new EnumParameter("dropErrors");
        dropErrors.registerEntry("Fail on first dataset with errors").setDefault();
        dropErrors.registerEntry("Drop results from datasets with errors");
        numThreadParam.setAdditionalDescription("If 0, the first dataset encountered.\n" +
                "\t   Generally a single NUMA node is 1 physical processor.\n" +
                "\t   The default will try to use hyperthreading cores, which can increase the amount of time this process will take.\n" +
                "\t   This is because the part of Scoring param generation that is multithreaded is also I/O intensive.");
        addParameter(dropErrors);

        addExample("Example (high-precision): java -Xmx4G -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen -i resultsFolder -d spectraFolder -m 2 -e 1 -protocol 5 -thread 4 -dropErrors 1");

        EnumParameter mgfParam = new EnumParameter("mgf");
        mgfParam.registerEntry("Do not create annotated mgf").setDefault();
        mgfParam.registerEntry("Create annotated mgf");
        mgfParam.setHidden();
        addParameter(mgfParam);


//		paramManager.addModFileParam();
//		StringParameter nlParam = new StringParameter("nl", "NeutralLosses", "Comma separated neutral losses to consider. Specify compositions or masses");
//		nlParam.setAdditionalDescription("E.g. '-nl H3PO4', '-nl 97.995,64.064'");
//		nlParam.defaultValue(null);
//		paramManager.addParameter(nlParam);

    }

    @Deprecated
    public void addMSGFDBParams() {
        addSpecFileParam();
        addDBFileParam();
        addPMTolParam();
        addOutputFileParam();

        IntParameter numThreadParam = new IntParameter("thread", "NumThreads", "Number of concurrent threads to be executed, Default: Number of available cores");
        numThreadParam.defaultValue(Runtime.getRuntime().availableProcessors());
        numThreadParam.minValue(1);
        addParameter(numThreadParam);

        EnumParameter tdaParam = new EnumParameter("tda");
        tdaParam.registerEntry("Don't search decoy database").setDefault();
        tdaParam.registerEntry("Search decoy database to compute FDR");
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

//		FloatRangeParameter itraqParam = new FloatRangeParameter("itraq", "minMass,maxMass", "Remove MS/MS peaks in the mass range between minMass and maxMass (for iTRAQ analysis).");
//		itraqParam.minValue(0f);
//		itraqParam.setMaxInclusive();
//		itraqParam.defaultValue("0,0");
//		itraqParam.setHidden();
//		addParameter(itraqParam);

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

        EnumParameter uniformAAProb = new EnumParameter(ParamNameEnum.UNIFORM_AA_PROBABILITY.commandlineName, ParamNameEnum.UNIFORM_AA_PROBABILITY.name);
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
        specIndexParam.defaultValue("1," + (Integer.MAX_VALUE - 1));
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

//		EnumParameter percolatorParam = new EnumParameter("percolator");
//		edgeScoreParam.registerEntry("normal").setDefault();
//		edgeScoreParam.registerEntry("for MS-GF+Percolator");
//		edgeScoreParam.setHidden();
//		addParameter(percolatorParam);
    }

    public void addMSGFParams() {
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

    public void addMSGFLibParams() {
        addSpecFileParam();

        // Add library file param
        FileParameter libFileParam = new FileParameter("d", "LibraryFile", "*.sptxt");
        libFileParam.addFileFormat(new FileFormat(".sptxt"));
        libFileParam.fileMustExist();
        libFileParam.mustBeAFile();
        addParameter(libFileParam);

        addPMTolParam();
        addOutputFileParam();

        IntParameter numThreadParam = new IntParameter(ParamNameEnum.NUM_THREADS.commandlineName, ParamNameEnum.NUM_THREADS.name, ParamNameEnum.NUM_THREADS.description);

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

    public FileParameter getSpecFileParam() {
        return ((FileParameter) getParameter("s"));
    }

    public FileParameter getDBFileParam() {
        return ((FileParameter) getParameter("d"));
    }

    public String getDecoyProteinPrefix() {
        StringParameter decoyProteinPrefixParam = (StringParameter) getParameter("decoy");
        return (decoyProteinPrefixParam.value);
    }

    public ToleranceParameter getPMTolParam() {
        return ((ToleranceParameter) getParameter("t"));
    }

    public FileParameter getOutputFileParam() {
        return ((FileParameter) getParameter("o"));
    }

    public ActivationMethod getActivationMethod() {
        return (ActivationMethod) ((ObjectEnumParameter<?>) getParameter("m")).getObject();
    }

    public InstrumentType getInstType() {
        return (InstrumentType) ((ObjectEnumParameter<?>) getParameter("inst")).getObject();
    }

    public Enzyme getEnzyme() {
        return (Enzyme) ((ObjectEnumParameter<?>) getParameter("e")).getObject();
    }

    public Protocol getProtocol() {
        return (Protocol) ((ObjectEnumParameter<?>) getParameter("protocol")).getObject();
    }

    public FileParameter getModFileParam() {
        return ((FileParameter) getParameter("mod"));
    }

    public FileParameter getConfigFileParam(){
        return ((FileParameter) getParameter("conf"));
    }

    public int getIntValue(String key) {
        Parameter param = this.getParameter(key);
        if (param instanceof IntParameter)
            return ((IntParameter) param).getValue();
        else {
            System.err.println("[Error] in ParamManager.getIntValue: " + key + " is not an instance of IntParameter.");
            System.exit(-1);
        }
        return -1;
    }

    public float getFloatValue(String key) {
        Parameter param = this.getParameter(key);
        if (param instanceof FloatParameter)
            return ((FloatParameter) param).getValue();
        else {
            System.err.println("[Error] in ParamManager.getFloatValue: " + key + " is not an instance of FloatParameter.");
            System.exit(-1);
        }
        return -1;
    }

    public double getDoubleValue(String key) {
        Parameter param = this.getParameter(key);
        if (param instanceof DoubleParameter)
            return ((DoubleParameter) param).getValue();
        else {
            System.err.println("[Error] in ParamManager.getDoubleValue: " + key + " is not an instance of DoubleParameter.");
            System.exit(-1);
        }
        return -1;
    }

    public File getFile(String key) {
        Parameter param = this.getParameter(key);
        if (param instanceof FileParameter)
            return ((FileParameter) param).getFile();
        else {
            System.err.println("[Error] in ParamManager.getFile: " + key + " is not an instance of FileParameter.");
            System.exit(-1);
        }
        return null;
    }

    public File[] getFiles(String key) {
        Parameter param = this.getParameter(key);
        if (param instanceof FileListParameter)
            return ((FileListParameter) param).getFiles();
        else {
            System.err.println("[Error] in ParamManager.getFile: " + key + " is not an instance of FileListParameter.");
            System.exit(-1);
        }
        return null;
    }

    public static void main(String argv[]) {
        ParamManager paramManager = new ParamManager("MSGF", MSGF.VERSION, MSGF.RELEASE_DATE, "java -Xmx2000M -jar MSGFDB.jar");
        paramManager.addMSGFDBParams();

//		FileListParameter testParam = new FileListParameter("test", "test", "test");
//		testParam.addFileFormat(SpecFileFormat.MGF);
//		testParam.addFileFormat(SpecFileFormat.MZXML);
//		paramManager.addParameter(testParam);

        String errMessage = paramManager.parseParams(argv);
        if (errMessage == null) {
            paramManager.printValues();
        } else {
            System.out.println();
            System.err.println("[Error] " + errMessage);
            System.out.println();
            paramManager.printUsageInfo();
        }
    }
}
