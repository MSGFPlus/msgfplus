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

    /**
     * Keys in this HashMap are the parameter key (typically the command line names), values are the parameter definition
     */
    private CaseInsensitiveLinkedHashMapParam params;

    private String toolName;
    private String version;
    private String date;
    private String command;
    private ArrayList<String> examples = new ArrayList<>();

    public enum ParamNameEnum {

        CONFIGURATION_FILE("conf", "ConfigurationFile",
                "Configuration file path; options specified at the command line will override settings in the config file",
                "Example parameter file is at https://github.com/MSGFPlus/msgfplus/blob/master/docs/examples/MSGFPlus_Params.txt"),

        SPECTRUM_FILE("s", "SpectrumFile", "*.mzML, *.mzXML, *.mgf, *.ms2, *.pkl or *_dta.txt",
                "Spectra should be centroided (see below for MSConvert example). Profile spectra will be ignored."),

        DB_FILE("d", "DatabaseFile", "*.fasta or *.fa or *.faa", null),

        DECOY_PREFIX("decoy", "DecoyPrefix",
                "Prefix for decoy protein names; Default: " + MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX, null),

        // Used by MS-GF+
        MZID_OUTPUT_FILE("o", "OutputFile (*.mzid)", "Default: [SpectrumFileName].mzid", null),

        // Used by MSGF and MS-GFDB
        OUTPUT_FILE("o", "OutputFile", "Default: stdout", null),

        //  MS-GF+, MSGF, and MS-GFDB
        PRECURSOR_MASS_TOLERANCE("t", "PrecursorMassTolerance", "e.g. 2.5Da, 20ppm or 0.5Da,2.5Da; Default: 20ppm",
                "Use a comma to define asymmetric values. E.g. \"-t 0.5Da,2.5Da\" will set 0.5Da to the left (ObsMass < TheoMass) and 2.5Da to the right (ObsMass > TheoMass)"),

        PRECURSOR_MASS_TOLERANCE_UNITS("u", "PrecursorMassToleranceUnits", "Units for the precursor mass tolerance; only useful if you do not include units in the PrecursorMassTolerance specification",
                "0 means Ds\n" +
                        "\t   1 means ppm\n" +
                        "\t   2 means use units specified by the PrecursorMassTolerance (Default)"),

        // aka Activation method
        FRAG_METHOD("m", "FragmentationMethodID", "Fragmentation Method",
                "0 means as written in the spectrum or CID if no info (Default)\n" +
                        "\t   1 means CID\n" +
                        "\t   2 means ETD\n" +
                        "\t   3 means HCD"),

        INSTRUMENT_TYPE("inst", "InstrumentID", null, null),

        ENZYME_ID("e", "EnzymeID", null, null),

        PROTOCOL_ID("protocol", "ProtocolID", null, null),

        MOD_FILE("mod", "ModificationFileName", "Modification file; Default: standard amino acids with fixed C+57; only if -mod is not specified", null),

        NUM_THREADS("thread", "NumThreads", "Number of concurrent threads to be executed; Default: Number of available cores",
                "This is best set to the number of physical cores in a single NUMA node.\n" +
                "\t   Generally a single NUMA node is 1 physical processor.\n" +
                "\t   The default will try to use hyperthreading cores, which can increase the amount of time this process will take.\n" +
                "\t   This is because the part of Scoring param generation that is multithreaded is also I/O intensive."),

        NUM_TASKS("tasks", "NumTasks", "Override the number of tasks to use on the threads; Default: (internally calculated based on inputs)",
                "More tasks than threads will reduce the memory requirements of the search, but will be slower (how much depends on the inputs).\n" +
                "\t   1 <= tasks <= numThreads: will create one task per thread, which is the original behavior.\n" +
                "\t   tasks = 0: use default calculation - minimum of: (threads*3) and (numSpectra/250).\n" +
                "\t   tasks < 0: multiply number of threads by abs(tasks) to determine number of tasks (i.e., -2 means \"2 * numThreads\" tasks).\n" +
                "\t   One task per thread will use the most memory, but will usually finish the fastest.\n" +
                "\t   2-3 tasks per thread will use comparably less memory, but may cause the search to take 1.5 to 2 times as long."),

        // Used by MS-GF+
        ISOTOPE_ERROR("ti", "IsotopeErrorRange", "Range of allowed isotope peak errors; Default: 0,1",
                "Takes into account the error introduced by choosing a non-monoisotopic peak for fragmentation.\n" +
                "\t   The combination of -t and -ti determines the precursor mass tolerance.\n" +
                "\t   E.g. \"-t 20ppm -ti -1,2\" tests abs(ObservedPepMass - TheoreticalPepMass - n * 1.00335Da) < 20ppm for n = -1, 0, 1, 2."),

        ENZYME_SPECIFICITY("ntt", "NTT", "Number of Tolerable Termini",
                "E.g. For trypsin, 0: non-tryptic, 1: semi-tryptic, 2: fully-tryptic peptides only."),

        // Used by MS-GFDB
        C13("c13", null, "Precursor isotope peak error",
                "0 means consider only peptides matching precursor mass\n" +
                        "\t   1 means Consider peptides having one 13C (Default)\n" +
                        "\t   2 means Consider peptides having up to two 13C"),

        // Used by MS-GFDB
        NNET("nnet", null, "Number of allowed non-enzymatic termini", null),

        MIN_PEPTIDE_LENGTH("minLength", "MinPepLength", "Minimum peptide length to consider; Default: 6", null),
        MAX_PEPTIDE_LENGTH("maxLength", "MaxPepLength", "Maximum peptide length to consider; Default: 40", null),

        MIN_CHARGE("minCharge", "MinCharge", "Minimum precursor charge to consider if charges are not specified in the spectrum file; Default: 2", null),
        MAX_CHARGE("maxCharge", "MaxCharge", "Maximum precursor charge to consider if charges are not specified in the spectrum file; Default: 3", null),

        NUM_MATCHES_SPEC("n", "NumMatchesPerSpec", "Number of matches per spectrum to be reported; Default: 1", null),

        CHARGE_CARRIER_MASSES("ccm", "ChargeCarrierMass", "Mass of charge carrier; Default: mass of proton (1.00727649)", null),

        MIN_NUM_PEAKS("minNumPeaks", "MinNumPeaksPerSpectrum", "Minimum number of peaks per spectrum; Default: " + Constants.MIN_NUM_PEAKS_PER_SPECTRUM, null),

        NUM_ISOFORMS("iso", "NumIsoforms", "Number of isoforms to consider per peptide; Default: 128" + Constants.NUM_VARIANTS_PER_PEPTIDE, null),

        IGNORE_MET_CLEAVAGE("ignoreMetCleavage", "IgnoreMetCleavage", "When 1, ignore N-terminal methionine cleavage",
                "0 means Consider protein N-term Met cleavage (Default)\n" +
                "\t   1 means Ignore protein N-term Met cleavage"),

        MIN_DE_NOVO_SCORE("minDeNovoScore", "MinDeNovoScore", "Minimum de Novo score; Default: " + Constants.MIN_DE_NOVO_SCORE, null),

        SPEC_INDEX("index", "SpecIndex", "Range of spectrum indices to be considered",
                "For example, to analyze the first 1000 spectra use -index 1,1000"),

        MAX_MISSED_CLEAVAGES("maxMissedCleavages", "MaxMissedCleavages", "Exclude peptides with more than this number of missed cleavages from the search; Default: -1 (no limit)", null),

        TDA_STRATEGY("tda", "TDA", "Target decoy strategy",
                "0 means Don't search decoy database (Default)\n" +
                "\t   1 means search the decoy database (forward + reverse proteins)"),

        ADD_FEATURES("addFeatures", "AddFeatures", "Include additional features in the output (enable this to post-process results with Percolator)",
                "0 means Output basic scores only (Default)\n" +
                        "\t   1 means Output additional features"),

        DD_DIRECTORY("dd", "DBIndexDir", "Path to the directory containing database index files", null),

        EDGE_SCORE("edgeScore", "EdgeScore", "Toggle edge scoring",
                "0 means Use Edge Scoring (Default)\n" +
                "\t   1 means Do not use edge scoring"),

        // Only used by MS-GFDB
        @Deprecated
        UNIFORM_AA_PROBABILITY("uniformAAProb", "UniformAAProb", null, null),

        MAX_NUM_MODS("numMods", "NumMods", "Maximum number of dynamic (variable) modifications per peptide; Default: 3", null),

        // Note that static and dynamic modifications cannot be specified at the command line
        // Use -mod or -conf
        STATIC_MODIFICATION("staticMod", "StaticMod", "Static/Fixed modification", null),

        DYNAMIC_MODIFICATION("dynamicMod", "DynamicMod", "Dynamic/Variable modification", null),

        CUSTOM_AA("customAA", "CustomAA", "Custom amino acid", null),

        VERBOSE("verbose", null, "Console output message verbosity",
                "0 means Report total progress only\n" +
                        "\t   1 means Report total and per-thread progress/status");

        private String key;
        private String name;
        private String description;
        private String additionalDescription;

        ParamNameEnum(String key, String name, String description, String additionalDescription) {
            this.key = key;
            this.name = name;
            this.description = description;
            this.additionalDescription = additionalDescription;
        }

        /**
         * Parameter key; defines the command line argument for this parameter
         * @return
         */
        public String getKey() {
            return key;
        }

        /**
         * Parameter name when used in a configuration file
         * @return
         */
        public String getName() {
            return name;
        }

        /**
         * Parameter description
         * @return
         */
        public String getDescription() {
            return description;
        }

        /**
         * Additional description
         * @return
         */
        public String getAdditionalDescription() {
            return additionalDescription;
        }

        /**
         * Check whether the parameter line matches this parameter's name
         * @param paramName Parameter name from the config file
         * @return True if it matches the parameter name of this class (more specifically, of a class that inherits from this class)
         */
        public boolean isThisParam(String paramName) {
            return ((getName() != null && paramName.equalsIgnoreCase(getName())));
        }

        public static String getParamNameFromLine(String lineSetting) {
            String[] lineParts = lineSetting.split("=");
            if (lineParts.length < 2)
                return "";

            String paramName = lineParts[0].trim();

            // Auto-update some names to change from abbreviations / alternate names to the standard name
            if (paramName.equalsIgnoreCase("IsotopeError")) {
                paramName = "IsotopeErrorRange";
            } else if (paramName.equalsIgnoreCase("TargetDecoyAnalysis")) {
                paramName = "TDA";
            } else if (paramName.equalsIgnoreCase("FragmentationMethod")) {
                paramName = "FragmentationMethodID";
            } else if (paramName.equalsIgnoreCase("Instrument")) {
                paramName = "InstrumentID";
            } else if (paramName.equalsIgnoreCase("Enzyme")) {
                paramName = "EnzymeID";
            } else if (paramName.equalsIgnoreCase("Protocol")) {
                paramName = "ProtocolID";
            } else if (paramName.equalsIgnoreCase("NumTolerableTermini")) {
                paramName = "NTT";
            } else if (paramName.equalsIgnoreCase("MinNumPeaks")) {
                paramName = "MinNumPeaksPerSpectrum";
            } else if (paramName.equalsIgnoreCase("MaxNumMods") || paramName.equalsIgnoreCase("MaxNumModsPerPeptide")) {
                paramName = "NumMods";
            } else if (paramName.equalsIgnoreCase("minLength") || paramName.equalsIgnoreCase("MinPeptideLength")) {
                paramName = "MinPepLength";
            } else if (paramName.equalsIgnoreCase("maxLength") || paramName.equalsIgnoreCase("MaxPeptideLength")) {
                paramName = "MaxPepLength";
            } else if (paramName.equalsIgnoreCase("PMTolerance") || paramName.equalsIgnoreCase("ParentMassTolerance")) {
                paramName = "PrecursorMassTolerance";
            }

            return paramName;
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

    private void addExample(String example) {
        this.examples.add(example);
    }

    public Parameter getParameter(String key) {
        return params.get(key);
    }

    /**
     * Validates that required parameters are defined
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
    
    public void printJVMInfo() {
        System.out.println("Java " + System.getProperty("java.version") + " (" + System.getProperty("java.vendor") + ")");
        System.out.println(System.getProperty("os.name") + " (" + System.getProperty("os.arch") + ", version " + System.getProperty("os.version") + ")");
    }

    public void printUsageInfo() {
        System.out.println();
        System.out.println(this.toolName + " " + this.version + " (" + this.date + ")");
        System.out.println();
        System.out.println("Usage: " + this.command);

        ArrayList<Parameter> optParams = new ArrayList<>();
        Iterator<Entry<String, Parameter>> itr = params.entrySet().iterator();
        while (itr.hasNext()) {
            Entry<String, Parameter> entry = itr.next();
            Parameter param = entry.getValue();
            if (!param.isHidden()) {
                if (!param.isOptional()) {
                    System.out.println("\t" + param);
                    if (param.getAdditionalDescription() != null)
                        System.out.println("\t   " + param.getAdditionalDescription());
                } else {
                    optParams.add(param);
                }
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

    public void addSpecFileParam(boolean isOptional) {
        FileParameter specFileParam = new FileParameter(ParamNameEnum.SPECTRUM_FILE);
        if (isOptional) {
            specFileParam.setAsOptional();
        }
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

    private void addDBFileParam(boolean isOptional) {
        addDBFileParam(ParamNameEnum.DB_FILE, isOptional);
    }

    private void addDBFileParam(ParamNameEnum paramInfo, boolean isOptional) {
        FileParameter dbFileParam = new FileParameter(paramInfo);
        if (isOptional) {
            dbFileParam.setAsOptional();
        }
        dbFileParam.addFileFormat(DBFileFormat.FASTA);
        dbFileParam.fileMustExist();
        dbFileParam.mustBeAFile();
        addParameter(dbFileParam);
    }

    private void addDBFileParam(String key, String description, boolean isOptional) {
        FileParameter dbFileParam = new FileParameter(key, ParamNameEnum.DB_FILE.name, description);
        if (isOptional)
            dbFileParam.setAsOptional();
        dbFileParam.addFileFormat(DBFileFormat.FASTA);
        dbFileParam.fileMustExist();
        dbFileParam.mustBeAFile();
        addParameter(dbFileParam);
    }

    private void addDecoyPrefixParam() {
        addDecoyPrefixParam(MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX);
    }

    private void addDecoyPrefixParam(String defaultDecoyPrefix) {
        StringParameter decoyPrefixParam = new StringParameter(ParamNameEnum.DECOY_PREFIX);
        // Note that defining a default value auto-sets isOptional to True
        decoyPrefixParam.defaultValue(defaultDecoyPrefix);
        addParameter(decoyPrefixParam);
    }

    private void addPrecursorMassToleranceParam() {
        ToleranceParameter pmTolParam = new ToleranceParameter(ParamNameEnum.PRECURSOR_MASS_TOLERANCE);
        pmTolParam.defaultValue("20ppm");
        addParameter(pmTolParam);
    }

    /**
     * -o for MS-GF+
     */
    private void addMzIdOutputFileParam() {
        FileParameter outputParam = new FileParameter(ParamNameEnum.MZID_OUTPUT_FILE);
        outputParam.addFileFormat(new FileFormat(".mzid").setCaseSensitive());
        outputParam.setAsOptional();
        addParameter(outputParam);
    }

    /**
     * -o for MSGF and MS-GFDB
     */
    private void addOutputFileParam() {
        FileParameter outputParam = new FileParameter(ParamNameEnum.OUTPUT_FILE);
        outputParam.setAsOptional();
        outputParam.fileMustNotExist();
        addParameter(outputParam);
    }

    /**
     * Used by both MS-GF+ and MS-GFDB
     * MS-GF+ passes True for doNotAddMergeMode, thus ignoring ActivationMethod.FUSION
     *
     * @param defaultMethod
     * @param doNotAddMergeMode
     */
    private void addFragMethodParam(ActivationMethod defaultMethod, boolean doNotAddMergeMode) {
        ObjectEnumParameter<ActivationMethod> fragParam = new ObjectEnumParameter<>(ParamNameEnum.FRAG_METHOD);
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

    private void addInstTypeParam() {
        addInstTypeParam(InstrumentType.LOW_RESOLUTION_LTQ);
    }

    private void addInstTypeParam(InstrumentType defaultInst) {
        ObjectEnumParameter<InstrumentType> instParam = new ObjectEnumParameter<InstrumentType>(ParamNameEnum.INSTRUMENT_TYPE);
        InstrumentType[] allInstTypes = InstrumentType.getAllRegisteredInstrumentTypes();
        for (InstrumentType inst : allInstTypes) {
            instParam.registerObject(inst);
            if (inst == defaultInst)
                instParam.setDefault();
        }
        addParameter(instParam);
    }

    private void addEnzymeParam() {
        addEnzymeParam(Enzyme.TRYPSIN);
    }

    private void addEnzymeParam(Enzyme enzymeId) {
        ObjectEnumParameter<Enzyme> enzParam = new ObjectEnumParameter<>(ParamNameEnum.ENZYME_ID);
        Enzyme[] allEnzymes = Enzyme.getAllRegisteredEnzymes();
        for (Enzyme e : allEnzymes) {
            enzParam.registerObject(e);
            if (e == enzymeId)
                enzParam.setDefault();
        }
        addParameter(enzParam);
    }

    private void addProtocolParam() {
        addProtocolParam(Protocol.AUTOMATIC);
    }

    private void addProtocolParam(Protocol defaultProtocol) {
        ObjectEnumParameter<Protocol> protocolParam = new ObjectEnumParameter<Protocol>(ParamNameEnum.PROTOCOL_ID);
        Protocol[] protocols = Protocol.getAllRegisteredProtocols();
        for (Protocol protocol : protocols) {
            protocolParam.registerObject(protocol);
            if (protocol == defaultProtocol)
                protocolParam.setDefault();
        }
        addParameter(protocolParam);
    }

    private void addEnzymeSpecificityParam() {
        EnumParameter nttParam = new EnumParameter(ParamNameEnum.ENZYME_SPECIFICITY);
        nttParam.registerEntry("");
        nttParam.registerEntry("");
        nttParam.registerEntry("").setDefault();
        addParameter(nttParam);
    }

    private void addModFileParam() {
        FileParameter modParam = new FileParameter(ParamNameEnum.MOD_FILE);
        modParam.setAsOptional();
        modParam.fileMustExist();
        addParameter(modParam);
    }

    private void addConfigFileParam() {
        FileParameter configFile = new FileParameter(ParamNameEnum.CONFIGURATION_FILE);
        configFile.setAsOptional();
        configFile.fileMustExist();
        addParameter(configFile);
    }

    private void addIsotopeRangeParam() {
        IntRangeParameter isotopeRange = new IntRangeParameter(ParamNameEnum.ISOTOPE_ERROR);
        isotopeRange.setMaxInclusive();
        isotopeRange.defaultValue("0,1");
        addParameter(isotopeRange);
    }

    private IntParameter addNumThreadsParam() {
        IntParameter numThreadsParam = new IntParameter(ParamNameEnum.NUM_THREADS);
        numThreadsParam.defaultValue(Runtime.getRuntime().availableProcessors());
        numThreadsParam.minValue(1);
        addParameter(numThreadsParam);
        return numThreadsParam;
    }

    private void addVerboseModeParam() {
        EnumParameter verboseOutputParam = new EnumParameter(ParamNameEnum.VERBOSE);
        verboseOutputParam.registerEntry("Report total progress only").setDefault();
        verboseOutputParam.registerEntry("Report total and per-thread progress/status");
        addParameter(verboseOutputParam);
    }

    private void addNumTasksParam() {
        IntParameter numTasksParam = new IntParameter(ParamNameEnum.NUM_TASKS);
        numTasksParam.defaultValue(0);
        numTasksParam.minValue(-10);
        addParameter(numTasksParam);
    }

    private void addTdaParam() {
        EnumParameter tdaParam = new EnumParameter(ParamNameEnum.TDA_STRATEGY);
        tdaParam.registerEntry("Don't search decoy database").setDefault();
        tdaParam.registerEntry("Search decoy database");
        addParameter(tdaParam);
    }

    private void addMinPeptideLengthParam() {
        IntParameter minLenParam = new IntParameter(ParamNameEnum.MIN_PEPTIDE_LENGTH);
        minLenParam.minValue(1);
        minLenParam.defaultValue(6);
        addParameter(minLenParam);
    }

    private void addMaxPeptideLengthParam() {
        IntParameter maxLenParam = new IntParameter(ParamNameEnum.MAX_PEPTIDE_LENGTH);
        maxLenParam.minValue(1);
        maxLenParam.defaultValue(40);
        addParameter(maxLenParam);
    }

    private void addMinChargeParam() {
        IntParameter minCharge = new IntParameter(ParamNameEnum.MIN_CHARGE);
        minCharge.minValue(1);
        minCharge.defaultValue(2);
        addParameter(minCharge);
    }

    private void addMaxChargeParam() {
        IntParameter maxCharge = new IntParameter(ParamNameEnum.MAX_CHARGE);
        maxCharge.minValue(1);
        maxCharge.defaultValue(3);
        addParameter(maxCharge);
    }

    private void addNumMatchesPerSpecParam() {
        IntParameter numMatchesParam = new IntParameter(ParamNameEnum.NUM_MATCHES_SPEC);
        numMatchesParam.minValue(1);
        numMatchesParam.defaultValue(1);
        addParameter(numMatchesParam);
    }

    private void addAddFeaturesParam() {
        EnumParameter addFeatureParam = new EnumParameter(ParamNameEnum.ADD_FEATURES);
        addFeatureParam.registerEntry("Output basic scores only").setDefault();
        addFeatureParam.registerEntry("Output additional features");
        addParameter(addFeatureParam);
    }

    private void addChargeCarrierMassParam() {
        DoubleParameter chargeCarrierMassParam = new DoubleParameter(ParamNameEnum.CHARGE_CARRIER_MASSES);
        chargeCarrierMassParam.minValue(0.1);
        chargeCarrierMassParam.setMaxInclusive();
        chargeCarrierMassParam.defaultValue(Composition.PROTON);
        addParameter(chargeCarrierMassParam);
    }

    private void addMaxMissedCleavagesParam() {
        IntParameter maxMissedCleavages = new IntParameter(ParamNameEnum.MAX_MISSED_CLEAVAGES);
        maxMissedCleavages.minValue(-1);
        maxMissedCleavages.defaultValue(-1);
        addParameter(maxMissedCleavages);
    }

    private void addMaxNumModsParam() {
        IntParameter maxNumMods = new IntParameter(ParamNameEnum.MAX_NUM_MODS);
        maxNumMods.minValue(0);
        maxNumMods.defaultValue(3);
        addParameter(maxNumMods);
    }

    private void addDbIndexDirParam(boolean isHidden) {
        FileParameter dbIndexDirParam = new FileParameter(ParamNameEnum.DD_DIRECTORY);
        dbIndexDirParam.fileMustExist();
        dbIndexDirParam.mustBeADirectory();
        dbIndexDirParam.setAsOptional();
        if (isHidden) {
            dbIndexDirParam.setHidden();
        }
        addParameter(dbIndexDirParam);
    }

    private void addPrecursorMassToleranceUnitsParam(boolean isHidden) {
        EnumParameter unitParam = new EnumParameter(ParamNameEnum.PRECURSOR_MASS_TOLERANCE_UNITS);
        unitParam.registerEntry("Da");
        unitParam.registerEntry("ppm");
        unitParam.registerEntry("Don't care").setDefault();
        if (isHidden) {
            unitParam.setHidden();
        }
        addParameter(unitParam);
    }

    private void addSpecIndexRangeParam(boolean isHidden) {
        IntRangeParameter specIndexParam = new IntRangeParameter(ParamNameEnum.SPEC_INDEX);
        specIndexParam.minValue(1);
        specIndexParam.setMaxInclusive();
        specIndexParam.defaultValue("1," + (Integer.MAX_VALUE - 1));
        if (isHidden) {
            specIndexParam.setHidden();
        }
        addParameter(specIndexParam);
    }

    private void addEdgeScoreParam(boolean isHidden) {
        EnumParameter edgeScoreParam = new EnumParameter(ParamNameEnum.EDGE_SCORE.key);
        edgeScoreParam.registerEntry("Use edge scoring").setDefault();
        edgeScoreParam.registerEntry("Do not use edge scoring");
        if (isHidden) {
            edgeScoreParam.setHidden();
        }
        addParameter(edgeScoreParam);
    }

    private void addMinNumPeaksParam(boolean isHidden) {
        IntParameter minNumPeaksParam = new IntParameter(ParamNameEnum.MIN_NUM_PEAKS);
        minNumPeaksParam.defaultValue(Constants.MIN_NUM_PEAKS_PER_SPECTRUM);
        if (isHidden) {
            minNumPeaksParam.setHidden();
        }
        addParameter(minNumPeaksParam);
    }

    private void addNumIsoformsParam(boolean isHidden) {
        IntParameter isoParam = new IntParameter(ParamNameEnum.NUM_ISOFORMS);
        isoParam.defaultValue(Constants.NUM_VARIANTS_PER_PEPTIDE);
        if (isHidden) {
            isoParam.setHidden();
        }
        addParameter(isoParam);
    }

    private void addMetCleavageParamParam(boolean isHidden) {
        EnumParameter metCleavageParam = new EnumParameter(ParamNameEnum.IGNORE_MET_CLEAVAGE);
        metCleavageParam.registerEntry("Consider protein N-term Met cleavage").setDefault();
        metCleavageParam.registerEntry("Ignore protein N-term Met cleavage");
        if (isHidden) {
            metCleavageParam.setHidden();
        }
        addParameter(metCleavageParam);
    }

    private void addMinDeNovoScoreParam(boolean isHidden) {
        IntParameter minDeNovoScoreParam = new IntParameter(ParamNameEnum.MIN_DE_NOVO_SCORE);
        minDeNovoScoreParam.minValue(Integer.MIN_VALUE);
        minDeNovoScoreParam.defaultValue(Constants.MIN_DE_NOVO_SCORE);
        if (isHidden) {
            minDeNovoScoreParam.setHidden();
        }
        addParameter(minDeNovoScoreParam);
    }

    /**
     * Add parameters for MS-GF+
     */
    public void addMSGFPlusParams() {

        // -conf ConfigurationFileName
        addConfigFileParam();

        // -s SpectrumFile (*.mzML, *.mzXML, *.mgf, *.ms2, *.pkl or *_dta.txt)
        addSpecFileParam(true);

        // -d DatabaseFile (*.fasta or *.fa or *.faa)
        addDBFileParam(true);
        addDecoyPrefixParam();

        // [-o OutputFile (*.mzid)] (Default: [SpectrumFileName].mzid)
        addMzIdOutputFileParam();

        addPrecursorMassToleranceParam();
        addPrecursorMassToleranceUnitsParam(true);

        addIsotopeRangeParam();

        addNumThreadsParam();
        addNumTasksParam();
        addVerboseModeParam();

        addTdaParam();

        addFragMethodParam(ActivationMethod.ASWRITTEN, true);
        addInstTypeParam();
        addEnzymeParam();
        addProtocolParam();
        addEnzymeSpecificityParam();

        addModFileParam();

        addMinPeptideLengthParam();
        addMaxPeptideLengthParam();
        addMinChargeParam();
        addMaxChargeParam();

        addNumMatchesPerSpecParam();
        addAddFeaturesParam();
        addChargeCarrierMassParam();
        addMaxMissedCleavagesParam();
        addMaxNumModsParam();

        addExample("Example (high-precision): java -Xmx3500M -jar MSGFPlus.jar -s test.mzML -d IPI_human_3.79.fasta -inst 1 -t 20ppm -ti -1,2 -ntt 2 -tda 1 -o testMSGFPlus.mzid -mod Mods.txt");
        addExample("Example (low-precision):  java -Xmx3500M -jar MSGFPlus.jar -s test.mzML -d IPI_human_3.79.fasta -inst 0 -t 0.5Da,2.5Da    -ntt 2 -tda 1 -o testMSGFPlus.mzid -mod Mods.txt");

        // Hidden parameters
        addDbIndexDirParam(true);
        addSpecIndexRangeParam(true);
        addEdgeScoreParam(true);
        addMinNumPeaksParam(true);
        addNumIsoformsParam(true);
        addMetCleavageParamParam(true);
        addMinDeNovoScoreParam(true);

    } // MSGFPlusParams

    /**
     * Used by class edu.ucsd.msjava.ui.ScoringParamGen
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
        ObjectEnumParameter<ActivationMethod> fragParam = new ObjectEnumParameter<>(ParamNameEnum.FRAG_METHOD);
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
        addEnzymeParam();

        // Protocol
        addProtocolParam();

        IntParameter numThreadsParam = addNumThreadsParam();
        numThreadsParam.defaultValue(Runtime.getRuntime().availableProcessors() / 2);

        EnumParameter dropErrors = new EnumParameter("dropErrors");
        dropErrors.setAdditionalDescription("If 0, stop processing if an error occurs; if 1, discard results from datasets with errors.");
        dropErrors.registerEntry("Fail on first dataset with errors").setDefault();
        dropErrors.registerEntry("Drop results from datasets with errors");
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
        addSpecFileParam(false);
        addDBFileParam(false);

        addPrecursorMassToleranceParam();
        addPrecursorMassToleranceUnitsParam(true);

        addOutputFileParam();

        addNumThreadsParam();

        addTdaParam();

        addFragMethodParam(ActivationMethod.ASWRITTEN, false);
        addInstTypeParam();
        addEnzymeParam();
        addProtocolParam();

        EnumParameter c13Param = new EnumParameter(ParamNameEnum.C13);
        c13Param.registerEntry("Consider only peptides matching precursor mass");
        c13Param.registerEntry("Consider peptides having one 13C").setDefault();
        c13Param.registerEntry("Consider peptides having up to two 13C");
        addParameter(c13Param);

        EnumParameter nnetParam = new EnumParameter(ParamNameEnum.NNET);
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

        addMinPeptideLengthParam();
        addMaxPeptideLengthParam();
        addMinChargeParam();
        addMaxChargeParam();

        addNumMatchesPerSpecParam();

        EnumParameter uniformAAProb = new EnumParameter(ParamNameEnum.UNIFORM_AA_PROBABILITY);
        uniformAAProb.registerEntry("Use amino acid probabilities computed from the input database").setDefault();
        uniformAAProb.registerEntry("Use probability 0.05 for all amino acids");
        addParameter(uniformAAProb);

        addExample("Example (high-precision): java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 30ppm -c13 1 -nnet 0 -tda 1 -o testMSGFDB.tsv");
        addExample("Example (low-precision):  java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 0.5Da,2.5Da  -nnet 0 -tda 1 -o testMSGFDB.tsv");

        // Hidden parameters
        addDbIndexDirParam(true);
        addSpecIndexRangeParam(true);

        EnumParameter showFDRParam = new EnumParameter("showFDR");
        showFDRParam.registerEntry("Do not show FDRs");
        showFDRParam.registerEntry("Show FDRs").setDefault();
        showFDRParam.setHidden();
        addParameter(showFDRParam);

        EnumParameter showDecoyParam = new EnumParameter("showDecoy");
        showDecoyParam.registerEntry("Do not show decoy PSMs").setDefault();
        showDecoyParam.registerEntry("Show decoy PSMs");
        showDecoyParam.setHidden();
        addParameter(showDecoyParam);

        EnumParameter replicateMergedResParam = new EnumParameter("replicate");
        replicateMergedResParam.registerEntry("Show merged spectra").setDefault();
        replicateMergedResParam.registerEntry("Show individual spectra");
        replicateMergedResParam.setHidden();
        addParameter(replicateMergedResParam);

        addEdgeScoreParam(true);

//		EnumParameter percolatorParam = new EnumParameter("percolator");
//		edgeScoreParam.registerEntry("normal").setDefault();
//		edgeScoreParam.registerEntry("for MS-GF+Percolator");
//		edgeScoreParam.setHidden();
//		addParameter(percolatorParam);

    } // MSGFDBParams

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
        addSpecFileParam(false);

        // Add library file param
        FileParameter libFileParam = new FileParameter("d", "LibraryFile", "*.sptxt");
        libFileParam.addFileFormat(new FileFormat(".sptxt"));
        libFileParam.fileMustExist();
        libFileParam.mustBeAFile();
        addParameter(libFileParam);

        addPrecursorMassToleranceParam();

        addOutputFileParam();

        addNumThreadsParam();

        addFragMethodParam(ActivationMethod.ASWRITTEN, false);
        addInstTypeParam();
        addEnzymeParam();
        addProtocolParam();

        EnumParameter c13Param = new EnumParameter(ParamNameEnum.C13);
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
        return ((FileParameter) getParameter(ParamNameEnum.SPECTRUM_FILE.key));
    }

    public FileParameter getDBFileParam() {
        return ((FileParameter) getParameter(ParamNameEnum.DB_FILE.key));
    }

    public String getDecoyProteinPrefix() {
        StringParameter decoyProteinPrefixParam = (StringParameter)getParameter(ParamNameEnum.DECOY_PREFIX.key);
        return (decoyProteinPrefixParam.value);
    }

    public double getChargeCarrierMass() {
        return getDoubleValue(ParamNameEnum.CHARGE_CARRIER_MASSES.key);
    }

    public ToleranceParameter getPrecursorMassToleranceParam() {
        return ((ToleranceParameter) getParameter(ParamNameEnum.PRECURSOR_MASS_TOLERANCE.key));
    }

    public int getToleranceUnit() {
        return getIntValue(ParamNameEnum.PRECURSOR_MASS_TOLERANCE_UNITS.key);
    }

    public IntRangeParameter getIsotopeRangeParameter() {
        return (IntRangeParameter) getParameter(ParamNameEnum.ISOTOPE_ERROR.key);
    }

    public FileParameter getOutputFileParam() {
        return ((FileParameter) getParameter(ParamNameEnum.OUTPUT_FILE.key));
    }

    public ActivationMethod getActivationMethod() {
        return (ActivationMethod) ((ObjectEnumParameter<?>) getParameter(ParamNameEnum.FRAG_METHOD.key)).getObject();
    }

    public InstrumentType getInstType() {
        return (InstrumentType) ((ObjectEnumParameter<?>) getParameter(ParamNameEnum.INSTRUMENT_TYPE.key)).getObject();
    }

    public Enzyme getEnzyme() {
        return (Enzyme) ((ObjectEnumParameter<?>) getParameter(ParamNameEnum.ENZYME_ID.key)).getObject();
    }

    public int getNumTolerableTermini() {
        return getIntValue(ParamNameEnum.ENZYME_SPECIFICITY.key);
    }

    public int getNumMatchesPerSpectrum() {
        return getIntValue(ParamNameEnum.NUM_MATCHES_SPEC.key);
    }

    public IntRangeParameter getSpecIndexParameter() {
        return ((IntRangeParameter) getParameter(ParamNameEnum.SPEC_INDEX.key));
    }

    public int getTDA() {
        return getIntValue(ParamNameEnum.TDA_STRATEGY.key);
    }

    public int getIgnoreMetCleavage() {
        return getIntValue(ParamNameEnum.IGNORE_MET_CLEAVAGE.key);
    }

    public int getOutputAdditionalFeatures() {
        return getIntValue(ParamNameEnum.ADD_FEATURES.key);
    }

    public int getMinPeptideLength() {
        return getIntValue(ParamNameEnum.MIN_PEPTIDE_LENGTH.key);
    }

    public int getMaxPeptideLength() {
        return getIntValue(ParamNameEnum.MAX_PEPTIDE_LENGTH.key);
    }

    public int getMaxNumVariantsPerPeptide() {
        return getIntValue(ParamNameEnum.NUM_ISOFORMS.key);
    }

    public int getMinCharge() {
        return getIntValue(ParamNameEnum.MIN_CHARGE.key);
    }

    public int getMaxCharge() {
        return getIntValue(ParamNameEnum.MAX_CHARGE.key);
    }

    public int getNumThreads() {
        return getIntValue(ParamNameEnum.NUM_THREADS.key);
    }

    public int getNumTasks() {
        return getIntValue(ParamNameEnum.NUM_TASKS.key);
    }

    public int getVerboseFlag() {
        return getIntValue(ParamNameEnum.VERBOSE.key);
    }

    public int getEdgeScoreFlag() {
        return getIntValue(ParamNameEnum.EDGE_SCORE.key);
    }

    // Used by MS-GF+
    public File getDatabaseIndexDir() {
        return getFile("dd");
    }

    public int getMinNumPeaksPerSpectrum() {
        return getIntValue(ParamNameEnum.MIN_NUM_PEAKS.key);
    }

    public int getMinDeNovoScore() {
        return getIntValue(ParamNameEnum.MIN_DE_NOVO_SCORE.key);
    }

    public int getMaxMissedCleavages() {
        return getIntValue(ParamNameEnum.MAX_MISSED_CLEAVAGES.key);
    }

    public int getMaxNumModsPerPeptide() {
        Parameter param = this.getParameter(ParamNameEnum.MAX_NUM_MODS.key);
        if (param == null) {
            this.addMaxNumModsParam();
        }
        return getIntValue(ParamNameEnum.MAX_NUM_MODS.key);
    }

    public Protocol getProtocol() {
        return (Protocol) ((ObjectEnumParameter<?>) getParameter(ParamNameEnum.PROTOCOL_ID.key)).getObject();
    }

    public FileParameter getModFileParam() {
        return ((FileParameter) getParameter(ParamNameEnum.MOD_FILE.key));
    }

    // Used by MS-GF+
    public FileParameter getConfigFileParam() {
        return ((FileParameter) getParameter(ParamNameEnum.CONFIGURATION_FILE.key));
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

    public void setMaxNumMods(int numMods) {
        Parameter numModsParam = getParameter(ParamManager.ParamNameEnum.MAX_NUM_MODS.getKey());
        numModsParam.parse(String.valueOf(numMods));
    }

    // This class is not typically instantiated
    @Deprecated
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
