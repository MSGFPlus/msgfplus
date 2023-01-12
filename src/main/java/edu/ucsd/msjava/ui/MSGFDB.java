package edu.ucsd.msjava.ui;

import edu.ucsd.msjava.msdbsearch.*;
import edu.ucsd.msjava.msgf.MSGFDBResultGenerator;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.IntRangeParameter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.params.ToleranceParameter;
import edu.ucsd.msjava.sequences.Constants;

import java.io.*;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * This class is deprecated
 * Instead, use MSGFPlus
 */
@Deprecated
public class MSGFDB {
    public static final String VERSION = "8091";
    public static final String RELEASE_DATE = "08/06/2012";

    public static final String DECOY_PROTEIN_PREFIX = "XXX";
    public static final String DECOY_DB_EXTENSION = ".revConcat.fasta";

    public static void main(String argv[]) {
        long time = System.currentTimeMillis();

        ParamManager paramManager = new ParamManager("MSGFDB", MSGFDB.VERSION, MSGFDB.RELEASE_DATE, "java -Xmx2000M -jar MSGFDB.jar");
        paramManager.addMSGFDBParams();

        if (argv.length == 0) {
            paramManager.printUsageInfo();
            return;
        }

        // Parse parameters
        String errMessage = paramManager.parseParams(argv);
        if (errMessage != null) {
            System.err.println("[Error] " + errMessage);
            System.out.println();
            paramManager.printUsageInfo();
            return;
        }

        // Running MS-GFDB (deprecated)
        paramManager.printToolInfo();
        paramManager.printJVMInfo();
        String errorMessage = null;
        try {
            errorMessage = runMSGFDB(paramManager);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }

        if (errorMessage != null) {
            System.err.println("[Error] " + errorMessage);
            System.out.println();
            System.exit(-1);
        } else
            System.out.format("MS-GFDB complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis() - time) / (float) 1000);
    }

    public static String runMSGFDB(ParamManager paramManager) {
        // Spectrum file
        FileParameter specParam = paramManager.getSpecFileParam();
        File specPath = specParam.getFile();

        if (!specPath.exists()) {
            return "Spectrum file not found: " + specPath.getPath();
        }

        if (!specPath.isDirectory()) {
            // Spectrum format
            SpecFileFormat specFormat = (SpecFileFormat) specParam.getFileFormat();

            // Output file
            File outputFile = paramManager.getOutputFileParam().getFile();

            return runMSGFDB(specPath, specFormat, outputFile, paramManager);
        } else    // spectrum directory
        {
            for (File f : specPath.listFiles()) {
                SpecFileFormat specFormat = SpecFileFormat.getSpecFileFormat(f.getName());
                if (specParam.isSupported(specFormat)) {
                    System.out.println("\nProcessing " + f.getAbsolutePath());
                    String outputFileName = f.getName().substring(0, f.getName().lastIndexOf('.')) + ".tsv";
                    File outputFile = new File(outputFileName);
                    if (outputFile.exists())
                        return outputFile.getAbsolutePath() + " already exists!";
                    System.out.println("Writing results to " + outputFile.getAbsolutePath());
                    String errMsg = runMSGFDB(f, specFormat, outputFile, paramManager);
                    if (errMsg != null)
                        return errMsg;
                }
            }
            return null;
        }
    }

    private static String runMSGFDB(File specFile, SpecFileFormat specFormat, File outputFile, ParamManager paramManager) {
        long time = System.currentTimeMillis();

        // Verify that the output directory exists and can be written to
        File outputDirectory = outputFile.getParentFile();
        if (outputDirectory != null) {

            if (!outputDirectory.exists()) {
                System.out.println("Creating directory " + outputDirectory.getPath());
                boolean success = outputDirectory.mkdirs();
                if (!success) {
                    return "Unable to create the missing directory: " + outputDirectory.getPath();
                }
            } else if (!outputDirectory.isDirectory()) {
                return "Invalid output file path (file path instead of directory path?): " + outputDirectory.getPath();
            }

            // An easy way to test for write access is outputDirectory.canWrite()
            // However, on Windows this is not always accurate
            // Thus, create a temporary file then delete it
            try {
                File testFile = File.createTempFile("MSGFPlus", ".tmp", outputDirectory);
                testFile.delete();
            } catch (java.io.IOException e) {
                return "Cannot create files in the output directory: " + e.getMessage();
            } catch (SecurityException e) {
                return "Cannot create files in the output directory; permission denied for: " + outputDirectory.getPath();
            }
        }

        // DB file
        File databaseFile = paramManager.getDBFileParam().getFile();

        // Precursor mass tolerance
        ToleranceParameter tol = ((ToleranceParameter) paramManager.getParameter(ParamManager.ParamNameEnum.PRECURSOR_MASS_TOLERANCE.getKey()));
        Tolerance leftPrecursorMassTolerance = tol.getLeftTolerance();
        Tolerance rightPrecursorMassTolerance = tol.getRightTolerance();

        int toleranceUnit = paramManager.getIntValue((ParamManager.ParamNameEnum.PRECURSOR_MASS_TOLERANCE_UNITS.getKey()));
        if (toleranceUnit != 2) {
            boolean isTolerancePPM;
            isTolerancePPM = toleranceUnit != 0;
            leftPrecursorMassTolerance = new Tolerance(leftPrecursorMassTolerance.getValue(), isTolerancePPM);
            rightPrecursorMassTolerance = new Tolerance(rightPrecursorMassTolerance.getValue(), isTolerancePPM);
        }

        int numAllowedC13 = paramManager.getIntValue(ParamManager.ParamNameEnum.C13.getKey());
        if (rightPrecursorMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f)
            numAllowedC13 = 0;

        Enzyme enzyme = paramManager.getEnzyme();
        int numAllowedNonEnzymaticTermini = paramManager.getIntValue(ParamManager.ParamNameEnum.NNET.getKey());
        ActivationMethod activationMethod = paramManager.getActivationMethod();
        InstrumentType instType = paramManager.getInstType();
        if (activationMethod == ActivationMethod.HCD)
            instType = InstrumentType.HIGH_RESOLUTION_LTQ;

        Protocol protocol = paramManager.getProtocol();

        AminoAcidSet aaSet = null;
        File modFile = paramManager.getModFileParam().getFile();
        if (modFile == null) {
            aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
        } else {
            String modFileName = modFile.getName();
            String ext = modFileName.substring(modFileName.lastIndexOf('.') + 1);
            if (ext.equalsIgnoreCase("xml"))
                aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(modFile.getAbsolutePath());
            else
                aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getAbsolutePath(), paramManager);
            if (aaSet.containsPhosphorylation()) {
                protocol = Protocol.PHOSPHORYLATION;
            }
        }

        int numMatchesPerSpec = paramManager.getNumMatchesPerSpectrum();

        IntRangeParameter specIndexParam = paramManager.getSpecIndexParameter();
        int startSpecIndex = specIndexParam.getMin();
        int endSpecIndex = specIndexParam.getMax();

        boolean useTDA = paramManager.getIntValue(ParamManager.ParamNameEnum.TDA_STRATEGY.getKey()) == 1;
        boolean showFDR = paramManager.getIntValue("showFDR") == 1;
        boolean showDecoy = paramManager.getIntValue("showDecoy") == 1;

        int minPeptideLength = paramManager.getIntValue(ParamManager.ParamNameEnum.MIN_PEPTIDE_LENGTH.getKey());
        int maxPeptideLength = paramManager.getIntValue(ParamManager.ParamNameEnum.MAX_PEPTIDE_LENGTH.getKey());
        if (minPeptideLength > maxPeptideLength) {
            return "MinPepLength must not be larger than MaxPepLength";
        }

        int minCharge = paramManager.getIntValue(ParamManager.ParamNameEnum.MIN_CHARGE.getKey());
        int maxCharge = paramManager.getIntValue(ParamManager.ParamNameEnum.MAX_CHARGE.getKey());
        if (minCharge > maxCharge) {
            return "MinCharge must not be larger than MaxCharge";
        }

        int numThreads = paramManager.getIntValue(ParamManager.ParamNameEnum.NUM_THREADS.getKey());
        boolean useUniformAAProb = paramManager.getIntValue(ParamManager.ParamNameEnum.UNIFORM_AA_PROBABILITY.getKey()) == 1;
        boolean replicateMergedResults = paramManager.getIntValue("replicate") == 1;
        boolean doNotDseEdgeScore = paramManager.getIntValue(ParamManager.ParamNameEnum.EDGE_SCORE.getKey()) == 1;
        boolean allowDenseCentroidedPeaks = paramManager.getIntValue(ParamManager.ParamNameEnum.ALLOW_DENSE_CENTROIDED_PEAKS.getKey()) == 1;

        System.out.println("Loading database files...");
        File dbIndexDir = paramManager.getFile(ParamManager.ParamNameEnum.DD_DIRECTORY.getKey());
        if (dbIndexDir != null) {

            File newDBFile = new File(Paths.get(dbIndexDir.getAbsolutePath(), databaseFile.getName()).toString());
            if (!useTDA) {
                if (!newDBFile.exists()) {
                    System.out.println("Creating " + newDBFile.getAbsolutePath() + ".");
                    ReverseDB.copyDB(databaseFile.getAbsolutePath(), newDBFile.getAbsolutePath());
                }
            }
            databaseFile = newDBFile;
        }

        if (useTDA) {
            String dbFileName = databaseFile.getName();
            String concatDBFileName = dbFileName.substring(0, dbFileName.lastIndexOf('.')) + DECOY_DB_EXTENSION;

            String concatDBFilePath = Paths.get(databaseFile.getAbsoluteFile().getParent(), concatDBFileName).toString();
            File concatTargetDecoyDBFile = new File(concatDBFilePath);

            if (!concatTargetDecoyDBFile.exists()) {
                System.out.println("Creating " + concatTargetDecoyDBFile.getAbsolutePath() + ".");
                if (ReverseDB.reverseDB(databaseFile.getAbsolutePath(), concatTargetDecoyDBFile.getAbsolutePath(), true, DECOY_PROTEIN_PREFIX) == false) {
                    return "Cannot create a decoy database file!";
                }
            }
            databaseFile = concatTargetDecoyDBFile;
        }

        if (!useUniformAAProb)
            DBScanner.setAminoAcidProbabilities(databaseFile.getAbsolutePath(), aaSet);

        aaSet.registerEnzyme(enzyme);

        CompactFastaSequence fastaSequence = new CompactFastaSequence(databaseFile.getAbsolutePath()).truncateAnnotation();
        if (useTDA) {
            float ratioUniqueProteins = fastaSequence.getRatioUniqueProteins();
            if (ratioUniqueProteins < 0.5f) {
                fastaSequence.printTooManyDuplicateSequencesMessage(databaseFile.getName(), "MS-GFDB");
                System.exit(-1);
            }
        }

        CompactSuffixArray sa = new CompactSuffixArray(fastaSequence, maxPeptideLength);
        System.out.print("Loading database finished ");
        System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - time) / 1000);

        System.out.println("Reading spectra...");

        // Show a message of the form "Opening mzML file QC_Mam_19_01_PNNL_10_06Jan21_Arwen_WBEH-20-12-01.mzML"
        System.out.printf("Opening %s %s\n", specFormat.getPSIName(), specFile.getName());

        SpectraAccessor specAcc = new SpectraAccessor(specFile, specFormat);

        if (specAcc.getSpecMap() == null || specAcc.getSpecItr() == null)
            return "Error while parsing spectrum file: " + specFile.getPath();


        if (enzyme == null)
            numAllowedNonEnzymaticTermini = 2;

        // determine the number of spectra to be scanned together
        long maxMemory = Runtime.getRuntime().maxMemory() - sa.getSize() - 1 << 28;

        int avgPeptideMass = 2000;
        int numBytesPerMass = 12;
        int numSpecScannedTogether = (int) ((float) maxMemory / avgPeptideMass / numBytesPerMass);
        ArrayList<SpecKey> specKeyList = SpecKey.getSpecKeyList(specAcc.getSpecItr(), startSpecIndex, endSpecIndex, minCharge, maxCharge, activationMethod, Constants.MIN_NUM_PEAKS_PER_SPECTRUM, allowDenseCentroidedPeaks);
        int specSize = specKeyList.size();

        System.out.print("Reading spectra finished ");
        System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - time) / 1000);

        numThreads = Math.min(numThreads, Math.round((float) Math.min(specSize, numSpecScannedTogether) / 250));
        if (numThreads == 0)
            numThreads = 1;
        System.out.println("Using " + numThreads + (numThreads == 1 ? " thread." : " threads."));

        SpecDataType specDataType = new SpecDataType(activationMethod, instType, enzyme, protocol);
        int fromIndexGlobal = 0;

        List<MSGFDBResultGenerator.DBMatch> resultList = Collections.synchronizedList(new ArrayList<MSGFDBResultGenerator.DBMatch>());

        while (true) {
            if (fromIndexGlobal >= specSize)
                break;
            int toIndexGlobal = Math.min(specSize, fromIndexGlobal + numSpecScannedTogether);
            System.out.println("Spectrum " + fromIndexGlobal + "-" + (toIndexGlobal - 1) + " (total: " + specSize + ")");

            // Thread pool
            ExecutorService executor = Executors.newFixedThreadPool(numThreads);

            // Partition specKeyList
            int size = toIndexGlobal - fromIndexGlobal;
            int subListSize = size / numThreads;
            int residue = size % numThreads;

            int[] startIndex = new int[numThreads];
            int[] endIndex = new int[numThreads];

            for (int i = 0; i < numThreads; i++) {
                startIndex[i] = i > 0 ? endIndex[i - 1] : fromIndexGlobal;
                endIndex[i] = startIndex[i] + subListSize + (i < residue ? 1 : 0);
            }

            for (int i = 0; i < numThreads; i++) {
                ScoredSpectraMap specScanner = new ScoredSpectraMap(
                        specAcc,
                        Collections.synchronizedList(specKeyList.subList(startIndex[i], endIndex[i])),
                        leftPrecursorMassTolerance,
                        rightPrecursorMassTolerance,
                        numAllowedC13,
                        specDataType,
                        false
                );
                if (doNotDseEdgeScore)
                    specScanner.turnOffEdgeScoring();

                ConcurrentMSGFDB.RunMSGFDB msgfdbExecutor = new ConcurrentMSGFDB.RunMSGFDB(
                        specScanner,
                        sa,
                        enzyme,
                        aaSet,
                        numMatchesPerSpec,
                        minPeptideLength,
                        maxPeptideLength,
                        numAllowedNonEnzymaticTermini,
                        !useTDA,
                        resultList,
                        specFile.getName(),
                        replicateMergedResults
                );
                executor.execute(msgfdbExecutor);
            }

            executor.shutdown();
            while (!executor.isTerminated()) {
            }    // wait until all threads terminate

            fromIndexGlobal += numSpecScannedTogether;
        }

        time = System.currentTimeMillis();
        // Sort search results by spectral probabilities
        Collections.sort(resultList);

        // Write results

        String header =
                "#SpecFile\tSpecIndex\tScan#\t"
                        + "FragMethod\t"
                        + "Precursor\tPMError("
                        + (rightPrecursorMassTolerance.isTolerancePPM() ? "ppm" : "Da")
                        + ")\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecProb\tP-value";

        MSGFDBResultGenerator gen = new MSGFDBResultGenerator(header, resultList);

        if (showFDR && !useTDA && numMatchesPerSpec == 1) {
            PrintStream out = null;
            if (outputFile == null)
                out = System.out;
            else {
                try {
                    out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            System.out.println("Computing EFDRs...");
            gen.computeEFDR();
            System.out.print("Computing EFDRs finished");
            System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - time) / 1000);
            gen.writeResults(out, true, false);
            if (out != System.out)
                out.close();
        } else if (!showFDR || !useTDA) {
            PrintStream out = null;
            if (outputFile == null)
                out = System.out;
            else {
                try {
                    out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            gen.writeResults(out, false, false);
            if (out != System.out)
                out.close();
        } else {
            System.out.println("Computing FDRs...");
            try {
                File tempFile = null;
                if (outputFile != null) {
                    tempFile = new File(outputFile.getAbsolutePath() + ".temp.tsv");
                } else {
                    tempFile = File.createTempFile("MSGFDB", "tempResult");
                    tempFile.deleteOnExit();
                }
                PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(tempFile)));
                gen.writeResults(out, false, false);
                out.flush();
                out.close();
                int specFileCol = 0;
                int specIndexCol = 1;
                int pepCol = 7;
                int dbCol = 8;
                int scoreCol = 11;
                edu.ucsd.msjava.fdr.ComputeFDR.computeFDR(tempFile, null, scoreCol, false, "\t",
                        specFileCol, specIndexCol, pepCol, null, true, showDecoy,
                        true, dbCol, DECOY_PROTEIN_PREFIX,
                        1, 1, outputFile);

            } catch (IOException e) {
                e.printStackTrace();
            }
            System.out.print("Computing FDRs finished");
            System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - time) / 1000);
        }
        return null;
    }
}
