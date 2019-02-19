package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.params.*;
import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

import static edu.ucsd.msjava.msutil.Composition.POTASSIUM_CHARGE_CARRIER_MASS;
import static edu.ucsd.msjava.msutil.Composition.PROTON;
import static edu.ucsd.msjava.msutil.Composition.SODIUM_CHARGE_CARRIER_MASS;

public class SearchParams {
    private List<DBSearchIOFiles> dbSearchIOList;
    private File databaseFile;
    private String decoyProteinPrefix;
    private Tolerance leftPrecursorMassTolerance;
    private Tolerance rightPrecursorMassTolerance;
    private int minIsotopeError;
    private int maxIsotopeError;
    private Enzyme enzyme;
    private int numTolerableTermini;
    private ActivationMethod activationMethod;
    private InstrumentType instType;
    private Protocol protocol;
    private AminoAcidSet aaSet;
    private int numMatchesPerSpec;
    private int startSpecIndex;
    private int endSpecIndex;
    private boolean useTDA;
    private boolean ignoreMetCleavage;
    private int minPeptideLength;
    private int maxPeptideLength;
    private int maxNumVariantsPerPeptide;
    private int minCharge;
    private int maxCharge;
    private int numThreads;
    private int numTasks;
    private boolean verbose;
    private boolean doNotUseEdgeScore;
    private File dbIndexDir;
    private boolean outputAdditionalFeatures;
    private int minNumPeaksPerSpectrum;
    private int minDeNovoScore;
    private double chargeCarrierMass;
    private int maxMissedCleavages;
    private int maxNumMods;

    public SearchParams() {
    }

    // Used by MS-GF+
    public List<DBSearchIOFiles> getDBSearchIOList() {
        return dbSearchIOList;
    }

    // Used by MS-GF+
    public File getDatabaseFile() {
        return databaseFile;
    }

    // Used by MS-GF+
    public String getDecoyProteinPrefix() { return decoyProteinPrefix; }

    // Used by MS-GF+
    public Tolerance getLeftPrecursorMassTolerance() {
        return leftPrecursorMassTolerance;
    }

    // Used by MS-GF+
    public Tolerance getRightPrecursorMassTolerance() {
        return rightPrecursorMassTolerance;
    }

    // Used by MS-GF+
    public int getMinIsotopeError() {
        return minIsotopeError;
    }

    // Used by MS-GF+
    public int getMaxIsotopeError() {
        return maxIsotopeError;
    }

    // Used by MS-GF+
    public Enzyme getEnzyme() {
        return enzyme;
    }

    public int getNumTolerableTermini() {
        return numTolerableTermini;
    }

    // Used by MS-GF+
    public ActivationMethod getActivationMethod() {
        return activationMethod;
    }

    // Used by MS-GF+
    public InstrumentType getInstType() {
        return instType;
    }

    // Used by MS-GF+
    public Protocol getProtocol() {
        return protocol;
    }

    // Used by MS-GF+
    public AminoAcidSet getAASet() {
        return aaSet;
    }

    // Used by MS-GF+
    public int getNumMatchesPerSpec() {
        return numMatchesPerSpec;
    }

    // Used by MS-GF+
    public int getStartSpecIndex() {
        return startSpecIndex;
    }

    // Used by MS-GF+
    public int getEndSpecIndex() {
        return endSpecIndex;
    }

    // Used by MS-GF+
    public boolean useTDA() {
        return useTDA;
    }

    // Used by MS-GF+
    public boolean ignoreMetCleavage() {
        return ignoreMetCleavage;
    }

    // Used by MS-GF+
    public int getMinPeptideLength() {
        return minPeptideLength;
    }

    // Used by MS-GF+
    public int getMaxPeptideLength() {
        return maxPeptideLength;
    }

    // Used by MS-GF+
    public int getMaxNumVariantsPerPeptide() {
        return maxNumVariantsPerPeptide;
    }

    // Used by MS-GF+
    public int getMinCharge() {
        return minCharge;
    }

    // Used by MS-GF+
    public int getMaxCharge() {
        return maxCharge;
    }

    // Used by MS-GF+
    public int getNumThreads() {
        return numThreads;
    }

    public int getNumTasks() {
        return numTasks;
    }

    public boolean getVerbose() {
        return verbose;
    }

    // Used by MS-GF+
    public boolean doNotUseEdgeScore() {
        return doNotUseEdgeScore;
    }

    // Used by MS-GF+
    public File getDBIndexDir() {
        return dbIndexDir;
    }

    public boolean outputAdditionalFeatures() {
        return outputAdditionalFeatures;
    }

    // Used by MS-GF+
    public int getMinNumPeaksPerSpectrum() {
        return minNumPeaksPerSpectrum;
    }

    // Used by MS-GF+
    public int getMinDeNovoScore() {
        return minDeNovoScore;
    }

    public double getChargeCarrierMass() {
        return chargeCarrierMass;
    }

    // Used by MS-GF+
    public int getMaxMissedCleavages() {
        return maxMissedCleavages;
    }

    // Used by MS-GF+
    public String parse(ParamManager paramManager) {
        AminoAcidSet configAASet = null;
        FileParameter configFileParam = paramManager.getConfigFileParam();

        if (configFileParam != null && configFileParam.getFile() != null) {
            configAASet = parseConfigParamFile(paramManager);
        }

        // Charge carrier mass
        chargeCarrierMass = paramManager.getChargeCarrierMass();
        Composition.setChargeCarrierMass(chargeCarrierMass);

        // Spectrum file
        FileParameter specParam = paramManager.getSpecFileParam();
        File specPath = specParam.getFile();

        dbSearchIOList = new ArrayList<>();

        if (!specPath.isDirectory()) {
            // Spectrum format
            SpecFileFormat specFormat = (SpecFileFormat) specParam.getFileFormat();
            // Output file
            File outputFile = paramManager.getOutputFileParam().getFile();
            if (outputFile == null) {
                String outputFilePath = specPath.getPath().substring(0, specPath.getPath().lastIndexOf('.')) + ".mzid";
                outputFile = new File(outputFilePath);
//				if (outputFile.exists())
//					return outputFile.getPath() + " already exists!";
            }

            dbSearchIOList = new ArrayList<>();
            dbSearchIOList.add(new DBSearchIOFiles(specPath, specFormat, outputFile));
        } else    // spectrum directory
        {
            dbSearchIOList = new ArrayList<>();
            for (File f : specPath.listFiles()) {
                SpecFileFormat specFormat = SpecFileFormat.getSpecFileFormat(f.getName());
                if (specParam.isSupported(specFormat)) {
                    String outputFileName = f.getName().substring(0, f.getName().lastIndexOf('.')) + ".mzid";
                    File outputFile = new File(outputFileName);
//					if (outputFile.exists())
//						return outputFile.getPath() + " already exists!";
                    dbSearchIOList.add(new DBSearchIOFiles(f, specFormat, outputFile));
                }
            }
        }

        // FASTA file
        databaseFile = paramManager.getDBFileParam().getFile();

        decoyProteinPrefix = paramManager.getDecoyProteinPrefix();

        // Precursor mass tolerance
        ToleranceParameter tol = paramManager.getPrecursorMassToleranceParam();
        leftPrecursorMassTolerance = tol.getLeftTolerance();
        rightPrecursorMassTolerance = tol.getRightTolerance();

        int toleranceUnit = paramManager.getToleranceUnit();
        if (toleranceUnit != 2) {
            boolean isTolerancePPM;
            isTolerancePPM = toleranceUnit != 0;
            leftPrecursorMassTolerance = new Tolerance(leftPrecursorMassTolerance.getValue(), isTolerancePPM);
            rightPrecursorMassTolerance = new Tolerance(rightPrecursorMassTolerance.getValue(), isTolerancePPM);
        }

        IntRangeParameter isotopeParam = paramManager.getIsotopeRangeParameter();
        this.minIsotopeError = isotopeParam.getMin();
        this.maxIsotopeError = isotopeParam.getMax();

        if (rightPrecursorMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f ||
            leftPrecursorMassTolerance.getToleranceAsDa(1000, 2)  >= 0.5f) {
            minIsotopeError = maxIsotopeError = 0;
        }

        enzyme = paramManager.getEnzyme();
        numTolerableTermini = paramManager.getNumTolerableTermini();
        activationMethod = paramManager.getActivationMethod();
        instType = paramManager.getInstType();
        if (activationMethod == ActivationMethod.HCD && instType != InstrumentType.HIGH_RESOLUTION_LTQ && instType != InstrumentType.QEXACTIVE)
            instType = InstrumentType.QEXACTIVE;    // by default use Q-Exactive model for HCD

        protocol = paramManager.getProtocol();

        aaSet = null;
        File modFile = paramManager.getModFileParam().getFile();
        if (modFile == null && configAASet == null)
            aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
        else {
            if (modFile != null) {
                String modFileName = modFile.getName();
                String ext = modFileName.substring(modFileName.lastIndexOf('.') + 1);
                if (ext.equalsIgnoreCase("xml"))
                    aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(modFile.getPath());
                else
                    aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath(), paramManager);
            } else {
                aaSet = configAASet;
            }

            if (protocol == Protocol.AUTOMATIC) {
                if (aaSet.containsITRAQ()) {
                    if (aaSet.containsPhosphorylation())
                        protocol = Protocol.ITRAQPHOSPHO;
                    else
                        protocol = Protocol.ITRAQ;
                } else if (aaSet.containsTMT()) {
                    protocol = Protocol.TMT;
                } else {
                    if (aaSet.containsPhosphorylation())
                        protocol = Protocol.PHOSPHORYLATION;
                    else
                        protocol = Protocol.STANDARD;
                }
            }
        }

        numMatchesPerSpec = paramManager.getNumMatchesPerSpectrum();

        IntRangeParameter specIndexParam = paramManager.getSpecIndexParameter();
        startSpecIndex = specIndexParam.getMin();
        endSpecIndex = specIndexParam.getMax();

        useTDA = paramManager.getTDA() == 1;
        ignoreMetCleavage = paramManager.getIgnoreMetCleavage() == 1;
        outputAdditionalFeatures = paramManager.getOutputAdditionalFeatures() == 1;

        minPeptideLength = paramManager.getMinPeptideLength();
        maxPeptideLength = paramManager.getMaxPeptideLength();

        // Number of isoforms to consider per peptide, Default: 128
        maxNumVariantsPerPeptide = paramManager.getMaxNumVariantsPerPeptide();

        if (minPeptideLength > maxPeptideLength) {
            return "MinPepLength must not be larger than MaxPepLength";
        }

        minCharge = paramManager.getMinCharge();
        maxCharge = paramManager.getMaxCharge();
        if (minCharge > maxCharge) {
            return "MinCharge must not be larger than MaxCharge";
        }

        numThreads = paramManager.getNumThreads();
        numTasks = paramManager.getNumTasks();
        verbose = paramManager.getVerboseFlag() == 1;
        doNotUseEdgeScore = paramManager.getEdgeScoreFlag() == 1;

        dbIndexDir = paramManager.getDatabaseIndexDir();

        minNumPeaksPerSpectrum = paramManager.getMinNumPeaksPerSpectrum();

        minDeNovoScore = paramManager.getMinDeNovoScore();

        /* Make sure max missed cleavages is a valid value and that it is not
         * being mixed with an unspecific or no-cleave enzyme
         */
        maxMissedCleavages = paramManager.getMaxMissedCleavages();
        if (maxMissedCleavages > -1 && enzyme.getName().equals("UnspecificCleavage")) {
            return "Cannot specify a MaxMissedCleavages when using unspecific cleavage enzyme";
        } else if (maxMissedCleavages > -1 && enzyme.getName().equals("NoCleavage")) {
            return "Cannot specify a MaxMissedCleavages when using no cleavage enzyme";
        }

        maxNumMods = paramManager.getMaxNumModsPerPeptide();
        int maxNumModsCompare = aaSet.getMaxNumberOfVariableModificationsPerPeptide();

        if (maxNumMods != maxNumModsCompare) {
            System.err.println("Error, code bug: " +
                    "MaxNumModsPerPeptide tracked by the ParamManager does not match the value tracked by the AminoAcidSet: " +
                    maxNumMods + " vs. " + maxNumModsCompare);
            System.exit(-1);
        }

        return null;
    }

    // Used by MS-GF+
    private AminoAcidSet parseConfigParamFile(ParamManager paramManager) {

        BufferedLineReader reader = null;

        File paramFile = paramManager.getConfigFileParam().getFile();

        try {
            reader = new BufferedLineReader(paramFile.getPath());
        } catch (IOException e) {
            System.err.println("Error opening parameter file " + paramFile.getPath());
            e.printStackTrace();
            System.exit(-1);
        }

        String dataLine;
        int lineNum = 0;

        // Keys in this table are line numbers
        // Values are the text from the config file, after the equals sign, defining a custom amino acid
        Hashtable<Integer, String> customAAByLine = new Hashtable<>();

        // Keys in this table are line numbers
        // Values are the text from the config file, after the equals sign, defining a static or dynamic mod
        Hashtable<Integer, String> modsByLine = new Hashtable<>();

        // Parse the settings

        int invalidParameterCount = 0;

        assert reader != null;
        while ((dataLine = reader.readLine()) != null) {
            lineNum++;

            // Check for a comment in the line (starts with a #)
            String[] tokenArray = dataLine.split("#");
            if (tokenArray.length == 0) {
                continue;
            }

            String lineSetting = tokenArray[0].trim();
            if (lineSetting.length() == 0) {
                continue;
            }

            String paramName = ParamManager.ParamNameEnum.getParamNameFromLine(lineSetting);
            if (paramName.isEmpty()) {
                continue;
            }

            if (ParamManager.ParamNameEnum.DYNAMIC_MODIFICATION.isThisParam(paramName) ||
                ParamManager.ParamNameEnum.STATIC_MODIFICATION.isThisParam(paramName) ||
                ParamManager.ParamNameEnum.CUSTOM_AA.isThisParam(paramName)) {

                String value = lineSetting.split("=")[1].trim();
                if (!value.equalsIgnoreCase("none")) {
                    // Store the text after the equals sign
                    if (ParamManager.ParamNameEnum.CUSTOM_AA.isThisParam(paramName))
                        customAAByLine.put(lineNum, value);
                    else
                        modsByLine.put(lineNum, value);
                }
                continue;
            }

            boolean validParameter = false;
            for (ParamManager.ParamNameEnum param: ParamManager.ParamNameEnum.values()) {
                if (param.isThisParam(paramName)) {
                    Parameter commandLineParam = paramManager.getParameter(param.getKey());
                    if (commandLineParam != null) {
                        validParameter = true;
                        if (!commandLineParam.isValueAssigned()) {
                            String value = lineSetting.split("=")[1].trim();
                            String parseError = commandLineParam.parse(value);
                            if (parseError == null || parseError.isEmpty()) {
                                commandLineParam.setValueAssigned();
                                continue;
                            }

                            if (commandLineParam.getKey().equals(ParamManager.ParamNameEnum.NUM_THREADS.getKey()) &&
                                    value.equalsIgnoreCase("all")) {
                                // Config file has: NumThreads=All
                                // This is acceptable
                                continue;
                            }

                            System.err.println("Error parsing '" + lineSetting + "' in config file " +
                                    paramFile.getAbsolutePath() + ": " + parseError);
                            System.exit(-1);
                        }
                    }
                }
            }

            if (!validParameter) {
                System.out.println("Warning, unrecognized parameter '" + lineSetting + "' in config file " + paramFile.getName());
                invalidParameterCount++;
            }

        }

        if (invalidParameterCount > 0) {
            System.out.println("Valid parameters are described in the example parameter file at " +
                    "https://github.com/MSGFPlus/msgfplus/blob/master/docs/examples/MSGFPlus_Params.txt");
        }

        return AminoAcidSet.getAminoAcidSetFromList(paramFile.getName(), customAAByLine, modsByLine, paramManager);
    }

    @Override
    public String toString() {
        StringBuffer buf = new StringBuffer();

//		buf.append("Spectrum File(s):\n");
//		for(DBSearchIOFiles ioFile : this.dbSearchIOList)
//		{
//			buf.append("\t"+ioFile.getSpecFile().getAbsolutePath()+"\n");
//		}
//		buf.append("Database File: " + this.databaseFile.getAbsolutePath() + "\n");

        buf.append("\tPrecursorMassTolerance: ");
        if (leftPrecursorMassTolerance.equals(rightPrecursorMassTolerance)) {
            buf.append(leftPrecursorMassTolerance);
        } else {
            buf.append("[" + leftPrecursorMassTolerance + "," + rightPrecursorMassTolerance + "]");
        }
        buf.append("\n");

        buf.append("\tIsotopeError: " + this.minIsotopeError + "," + this.maxIsotopeError + "\n");
        buf.append("\tTargetDecoyAnalysis: " + this.useTDA + "\n");
        buf.append("\tFragmentationMethod: " + this.activationMethod + "\n");
        buf.append("\tInstrument: " + (instType == null ? "null" : this.instType.getNameAndDescription()) + "\n");
        buf.append("\tEnzyme: " + (enzyme == null ? "null" : this.enzyme.getName()) + "\n");

        String customEnzymeFile = Enzyme.getCustomEnzymeFilePath();
        if (customEnzymeFile != null && !customEnzymeFile.isEmpty()) {
            buf.append("\tEnzyme file: " + customEnzymeFile + "\n");
        }

        ArrayList<String> customEnzymeMessages = Enzyme.getCustomEnzymeMessages();
        for (String message : customEnzymeMessages) {
            buf.append("\tEnzyme info: " + message + "\n");
        }

        buf.append("\tProtocol: " + (protocol == null ? "null" : this.protocol.getName()) + "\n");
        buf.append("\tNumTolerableTermini: " + this.numTolerableTermini + "\n");
        buf.append("\tMinPepLength: " + this.minPeptideLength + "\n");
        buf.append("\tMaxPepLength: " + this.maxPeptideLength + "\n");
        buf.append("\tMinCharge: " + this.minCharge + "\n");
        buf.append("\tMaxCharge: " + this.maxCharge + "\n");
        buf.append("\tNumMatchesPerSpec: " + this.numMatchesPerSpec + "\n");
        buf.append("\tMaxMissedCleavages: " + this.maxMissedCleavages + "\n");
        buf.append("\tMaxNumModsPerPeptide: " + this.maxNumMods + "\n");
        buf.append("\tChargeCarrierMass: " + this.chargeCarrierMass);

        if (Math.abs(this.chargeCarrierMass - PROTON) < 0.005) {
            buf.append(" (proton)\n");
        } else if (Math.abs(this.chargeCarrierMass - POTASSIUM_CHARGE_CARRIER_MASS) < 0.005) {
            buf.append(" (potassium)\n");
        } else if (Math.abs(this.chargeCarrierMass - SODIUM_CHARGE_CARRIER_MASS) < 0.005) {
            buf.append(" (sodium)\n");
        } else {
            buf.append(" (custom)\n");
        }

        buf.append("\tMinNumPeaksPerSpectrum: " + this.minNumPeaksPerSpectrum + "\n");
        buf.append("\tNumIsoforms: " + this.maxNumVariantsPerPeptide + "\n");

        ArrayList<String> modificationsInUse = aaSet.getModificationsInUse();

        if (modificationsInUse.size() == 0) {
            buf.append("No static or dynamic post translational modifications are defined.\n");
        } else {
            buf.append("Post translational modifications in use:\n");
            for (String modInfo : modificationsInUse )
                buf.append("\t" + modInfo + "\n");
        }

        return buf.toString();
    }
}
