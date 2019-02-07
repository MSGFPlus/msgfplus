package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.params.*;
import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static edu.ucsd.msjava.msutil.Composition.POTASSIUM_CHARGE_CARRIER_MASS;
import static edu.ucsd.msjava.msutil.Composition.PROTON;
import static edu.ucsd.msjava.msutil.Composition.SODIUM_CHARGE_CARRIER_MASS;

public class SearchParams {
    private List<DBSearchIOFiles> dbSearchIOList;
    private File databaseFile;
    private String decoyProteinPrefix;
    private Tolerance leftParentMassTolerance;
    private Tolerance rightParentMassTolerance;
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
    //	private boolean showFDR;
//	private boolean showDecoy;
    private int minPeptideLength;
    private int maxPeptideLength;
    private int maxNumVariatsPerPeptide;
    private int minCharge;
    private int maxCharge;
    private int numThreads;
    private int numTasks;
    private boolean verbose;
    private boolean replicateMergedResults;
    private boolean doNotUseEdgeScore;
    private File dbIndexDir;
    private boolean outputAdditionalFeatures;
    private int minNumPeaksPerSpectrum;
    private int minDeNovoScore;
    private double chargeCarrierMass;
    private int maxMissedCleavages;

    public SearchParams() {
    }

    public List<DBSearchIOFiles> getDBSearchIOList() {
        return dbSearchIOList;
    }

    public File getDatabaseFile() {
        return databaseFile;
    }

    public String getDecoyProteinPrefix() {return decoyProteinPrefix; }

    public Tolerance getLeftParentMassTolerance() {
        return leftParentMassTolerance;
    }

    public Tolerance getRightParentMassTolerance() {
        return rightParentMassTolerance;
    }

    public int getMinIsotopeError() {
        return minIsotopeError;
    }

    public int getMaxIsotopeError() {
        return maxIsotopeError;
    }

    public Enzyme getEnzyme() {
        return enzyme;
    }

    public int getNumTolerableTermini() {
        return numTolerableTermini;
    }

    public ActivationMethod getActivationMethod() {
        return activationMethod;
    }

    public InstrumentType getInstType() {
        return instType;
    }

    public Protocol getProtocol() {
        return protocol;
    }

    public AminoAcidSet getAASet() {
        return aaSet;
    }

    public int getNumMatchesPerSpec() {
        return numMatchesPerSpec;
    }

    public int getStartSpecIndex() {
        return startSpecIndex;
    }

    public int getEndSpecIndex() {
        return endSpecIndex;
    }

    public boolean useTDA() {
        return useTDA;
    }

    public boolean ignoreMetCleavage() {
        return ignoreMetCleavage;
    }

    public int getMinPeptideLength() {
        return minPeptideLength;
    }

    public int getMaxPeptideLength() {
        return maxPeptideLength;
    }

    public int getMaxNumVariatsPerPeptide() {
        return maxNumVariatsPerPeptide;
    }

    public int getMinCharge() {
        return minCharge;
    }

    public int getMaxCharge() {
        return maxCharge;
    }

    public int getNumThreads() {
        return numThreads;
    }

    public int getNumTasks() {
        return numTasks;
    }

    public boolean getVerbose() {
        return verbose;
    }

    public boolean replicateMergedResults() {
        return replicateMergedResults;
    }

    public boolean doNotUseEdgeScore() {
        return doNotUseEdgeScore;
    }

    public File getDBIndexDir() {
        return dbIndexDir;
    }

    public boolean outputAdditionalFeatures() {
        return outputAdditionalFeatures;
    }

    public int getMinNumPeaksPerSpectrum() {
        return minNumPeaksPerSpectrum;
    }

    public int getMinDeNovoScore() {
        return minDeNovoScore;
    }

    public double getChargeCarrierMass() {
        return chargeCarrierMass;
    }

    public int getMaxMissedCleavages() {
        return maxMissedCleavages;
    }

    public String parse(ParamManager paramManager) {
        // Charge carrier mass
        AminoAcidSet configAASet = null;
        if(paramManager.getConfigFileParam() != null){
            configAASet = parseConfigParamFile(paramManager);
        }
        chargeCarrierMass = paramManager.getDoubleValue("ccm");
        Composition.setChargeCarrierMass(chargeCarrierMass);

        // Spectrum file
        FileParameter specParam = paramManager.getSpecFileParam();
        File specPath = specParam.getFile();

        dbSearchIOList = new ArrayList<DBSearchIOFiles>();

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

            dbSearchIOList = new ArrayList<DBSearchIOFiles>();
            dbSearchIOList.add(new DBSearchIOFiles(specPath, specFormat, outputFile));
        } else    // spectrum directory
        {
            dbSearchIOList = new ArrayList<DBSearchIOFiles>();
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

        // DB file
        databaseFile = paramManager.getDBFileParam().getFile();

        decoyProteinPrefix = paramManager.getDecoyProteinPrefix();

        // PM tolerance
        ToleranceParameter tol = ((ToleranceParameter) paramManager.getParameter("t"));
        leftParentMassTolerance = tol.getLeftTolerance();
        rightParentMassTolerance = tol.getRightTolerance();

        int toleranceUnit = paramManager.getIntValue("u");
        if (toleranceUnit != 2) {
            boolean isTolerancePPM;
            isTolerancePPM = toleranceUnit != 0;
            leftParentMassTolerance = new Tolerance(leftParentMassTolerance.getValue(), isTolerancePPM);
            rightParentMassTolerance = new Tolerance(rightParentMassTolerance.getValue(), isTolerancePPM);
        }

        IntRangeParameter isotopeParam = (IntRangeParameter) paramManager.getParameter("ti");
        this.minIsotopeError = isotopeParam.getMin();
        this.maxIsotopeError = isotopeParam.getMax();

        if (rightParentMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f || leftParentMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f)
            minIsotopeError = maxIsotopeError = 0;

        enzyme = paramManager.getEnzyme();
        numTolerableTermini = paramManager.getIntValue("ntt");
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
            if(modFile != null){
                String modFileName = modFile.getName();
                String ext = modFileName.substring(modFileName.lastIndexOf('.') + 1);
                if (ext.equalsIgnoreCase("xml"))
                    aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(modFile.getPath());
                else
                    aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
            }else
                aaSet = configAASet;

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

        numMatchesPerSpec = paramManager.getIntValue("n");

        startSpecIndex = ((IntRangeParameter) paramManager.getParameter("index")).getMin();
        endSpecIndex = ((IntRangeParameter) paramManager.getParameter("index")).getMax();

        useTDA = paramManager.getIntValue("tda") == 1;
        ignoreMetCleavage = paramManager.getIntValue("ignoreMetCleavage") == 1;
        outputAdditionalFeatures = paramManager.getIntValue("addFeatures") == 1;

        minPeptideLength = paramManager.getIntValue("minLength");
        maxPeptideLength = paramManager.getIntValue("maxLength");

        maxNumVariatsPerPeptide = paramManager.getIntValue("iso");

        if (minPeptideLength > maxPeptideLength) {
            return "MinPepLength must not be larger than MaxPepLength";
        }

        minCharge = paramManager.getIntValue("minCharge");
        maxCharge = paramManager.getIntValue("maxCharge");
        if (minCharge > maxCharge) {
            return "MinCharge must not be larger than MaxCharge";
        }

        numThreads = paramManager.getIntValue("thread");
        numTasks = paramManager.getIntValue("tasks");
        verbose = paramManager.getIntValue("verbose") == 1;
        doNotUseEdgeScore = paramManager.getIntValue("edgeScore") == 1;

        dbIndexDir = paramManager.getFile("dd");

        minNumPeaksPerSpectrum = paramManager.getIntValue("minNumPeaks");

        minDeNovoScore = paramManager.getIntValue("minDeNovoScore");

        /* Make sure max missed cleavages is valid value and that it is not
         * being mixed with an unspecific or no-cleave enzyme
         *
         * String comparison to name is fragile here. It would be better if
         * there was a stable identifier to use for the comparision.
         */
        maxMissedCleavages = paramManager.getIntValue("maxMissedCleavages");
        if (maxMissedCleavages > -1 && enzyme.getName().equals("UnspecificCleavage")) {
            return "Cannot specify a MaxMissedCleavages when using unspecific cleavage enzyme";
        } else if (maxMissedCleavages > -1 && enzyme.getName().equals("NoCleavage")) {
            return "Cannot specify a MaxMissedCleavages when using no cleavage enzyme";
        }

        return null;
    }

    private AminoAcidSet parseConfigParamFile(ParamManager paramManager) {

        BufferedLineReader reader = null;

        File paramFile = paramManager.getConfigFileParam().getFile();

        try {
            reader = new BufferedLineReader(paramFile.getPath());
        } catch (IOException e) {
            e.printStackTrace();
        }
        int numMods = 2;

        // parse modifications
        String dataLine;
        List<String> mods = new ArrayList<>();

        assert reader != null;
        while ((dataLine = reader.readLine()) != null) {
            String[] tokenArray = dataLine.split("#");
            if (tokenArray.length == 0)
                continue;

            String lineSetting = tokenArray[0].trim();
            if (lineSetting.length() == 0) {
                continue;
            } else if(ParamManager.ParamNameEnum.DYNAMIC_MODIFICATION.isLine(lineSetting.toLowerCase()) || ParamManager.ParamNameEnum.STATIC_MODIFICATION.isLine(lineSetting.toLowerCase()) || ParamManager.ParamNameEnum.CUSTOM_AA.isLine(lineSetting.toLowerCase())){
                String value = lineSetting.split("=")[1].trim();
                mods.add(value);
            }else {
                for(ParamManager.ParamNameEnum param: ParamManager.ParamNameEnum.values()){
                    if (param.isLine(lineSetting.toLowerCase()) && (paramManager.getParameter(param.getCommandlineName()) != null && !paramManager.getParameter(param.getCommandlineName()).isValueAssigned())){
                        String value = lineSetting.split("=")[1].trim();
                        Parameter currentParam = paramManager.getParameter(param.getCommandlineName());
                        currentParam.parse(value);
                        currentParam.setValueAssigned();
                    }
                }
                if (ParamManager.ParamNameEnum.MAX_NUM_MODS.isLine(lineSetting.toLowerCase())){
                    String value = lineSetting.split("=")[1].trim();
                    try{
                        numMods = Integer.parseInt(value);
                    }catch (NumberFormatException e){

                    }
                }
            }
        }
        return AminoAcidSet.getAminoAcidSetFromList(paramFile.getName(), mods, numMods);
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
        if (leftParentMassTolerance.equals(rightParentMassTolerance))
            buf.append(leftParentMassTolerance);
        else
            buf.append("[" + leftParentMassTolerance + "," + rightParentMassTolerance + "]");
        buf.append("\n");

        buf.append("\tIsotopeError: " + this.minIsotopeError + "," + this.maxIsotopeError + "\n");
        buf.append("\tTargetDecoyAnalysis: " + this.useTDA + "\n");
        buf.append("\tFragmentationMethod: " + this.activationMethod + "\n");
        buf.append("\tInstrument: " + (instType == null ? "null" : this.instType.getNameAndDescription()) + "\n");
        buf.append("\tEnzyme: " + (enzyme == null ? "null" : this.enzyme.getName()) + "\n");
        buf.append("\tProtocol: " + (protocol == null ? "null" : this.protocol.getName()) + "\n");
        buf.append("\tNumTolerableTermini: " + this.numTolerableTermini + "\n");
        buf.append("\tMinPeptideLength: " + this.minPeptideLength + "\n");
        buf.append("\tMaxPeptideLength: " + this.maxPeptideLength + "\n");
        buf.append("\tNumMatchesPerSpec: " + this.numMatchesPerSpec + "\n");
        buf.append("\tMaxMissedCleavages: " + this.maxMissedCleavages + "\n");
        buf.append("\tChargeCarrierMass: " + this.chargeCarrierMass);

        if (Math.abs(this.chargeCarrierMass - PROTON) < 0.005) {
            buf.append(" (proton)");
        } else if (Math.abs(this.chargeCarrierMass - POTASSIUM_CHARGE_CARRIER_MASS) < 0.005) {
            buf.append(" (potassium)");
        } else if (Math.abs(this.chargeCarrierMass - SODIUM_CHARGE_CARRIER_MASS) < 0.005) {
            buf.append(" (sodium)");
        } else {
            buf.append(" (custom)");
        }

        return buf.toString();
    }
}
