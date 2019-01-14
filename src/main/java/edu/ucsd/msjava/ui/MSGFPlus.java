package edu.ucsd.msjava.ui;

import edu.ucsd.msjava.fdr.ComputeFDR;
import edu.ucsd.msjava.misc.ThreadPoolExecutorWithExceptions;
import edu.ucsd.msjava.msdbsearch.*;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.mzid.MZIdentMLGen;
import edu.ucsd.msjava.mzml.MzMLAdapter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.sequences.Constants;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;


public class MSGFPlus {
    public static final String VERSION = "Release (v2019.01.11)";
    public static final String RELEASE_DATE = "11 January 2019";

    public static final String DECOY_DB_EXTENSION = ".revCat.fasta";
    public static final String DECOY_PROTEIN_PREFIX = "XXX";

    public static void main(String argv[]) {
        long time = System.currentTimeMillis();

        ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
        paramManager.addMSGFPlusParams();

        if (argv.length == 0) {
            paramManager.printUsageInfo();
            return;
        }

        MzMLAdapter.turnOffLogs();

        // Parse parameters
        String errMessage = paramManager.parseParams(argv);
        if (errMessage != null) {
            System.err.println("[Error] " + errMessage);
            System.out.println();
            paramManager.printUsageInfo();
            return;
        }

        // Running MS-GF+
        paramManager.printToolInfo();
        String errorMessage = null;
        try {
            errorMessage = runMSGFPlus(paramManager);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }

        if (errorMessage != null) {
            System.err.println("[Error] " + errorMessage);
            System.out.println();
            System.exit(-1);
        } else
            System.out.format("MS-GF+ complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis() - time) / (float) 1000);
    }

    public static String runMSGFPlus(ParamManager paramManager) {
        SearchParams params = new SearchParams();
        String errorMessage = params.parse(paramManager);
        if (errorMessage != null)
            return errorMessage;
        else {
            List<DBSearchIOFiles> ioList = params.getDBSearchIOList();
            boolean multFiles = false;
            if (ioList.size() >= 2) {
                System.out.println("Processing " + ioList.size() + " spectra");
                for (DBSearchIOFiles ioFiles : ioList)
                    System.out.println("\t" + ioFiles.getSpecFile().getName());
                multFiles = true;
            }


            int ioIndex = -1;
            for (DBSearchIOFiles ioFiles : ioList) {
                ++ioIndex;
                File specFile = ioFiles.getSpecFile();
                SpecFileFormat specFormat = ioFiles.getSpecFileFormat();
                File outputFile = ioFiles.getOutputFile();

                if (multFiles) {
                    if (!outputFile.exists()) {
                        System.out.println("\nProcessing " + specFile.getPath());
                        System.out.println("Writing results to " + outputFile.getPath());
                        String errMsg = runMSGFPlus(ioIndex, specFormat, outputFile, params);
                        if (errMsg != null)
                            return errMsg;
                    } else {
                        System.out.println("\nIgnoring " + specFile.getPath());
                        System.out.println("Output file " + outputFile.getPath() + " exists.");
                    }
                } else {
                    String errMsg = runMSGFPlus(ioIndex, specFormat, outputFile, params);
                    if (errMsg != null)
                        return errMsg;
                }
            }
        }
        return null;
    }

    private static String runMSGFPlus(int ioIndex, SpecFileFormat specFormat, File outputFile, SearchParams params) {
        long time = System.currentTimeMillis();

        // Check the outputFile is valid for writing
        File parent = outputFile.getParentFile();
        if (parent != null && !parent.exists()) {
            return "Cannot create " + outputFile.getPath() + "!";
        }

        // DB file
        File databaseFile = params.getDatabaseFile();

        // PM tolerance
        Tolerance leftParentMassTolerance = params.getLeftParentMassTolerance();
        Tolerance rightParentMassTolerance = params.getRightParentMassTolerance();

        int minIsotopeError = params.getMinIsotopeError();    // inclusive
        int maxIsotopeError = params.getMaxIsotopeError();    // inclusive

        Enzyme enzyme = params.getEnzyme();

        ActivationMethod activationMethod = params.getActivationMethod();
        InstrumentType instType = params.getInstType();
        Protocol protocol = params.getProtocol();

        AminoAcidSet aaSet = params.getAASet();

        int startSpecIndex = params.getStartSpecIndex();
        int endSpecIndex = params.getEndSpecIndex();

        boolean useTDA = params.useTDA();

        int minCharge = params.getMinCharge();
        int maxCharge = params.getMaxCharge();

        int numThreads = params.getNumThreads();
        boolean doNotUseEdgeScore = params.doNotUseEdgeScore();

        int minNumPeaksPerSpectrum = params.getMinNumPeaksPerSpectrum();
        if (minNumPeaksPerSpectrum == -1)    // not specified
        {
            if (instType == InstrumentType.TOF)
                minNumPeaksPerSpectrum = Constants.MIN_NUM_PEAKS_PER_SPECTRUM_TOF;
            else
                minNumPeaksPerSpectrum = Constants.MIN_NUM_PEAKS_PER_SPECTRUM;
        }

        System.out.println("Loading database files...");

        File dbIndexDir = params.getDBIndexDir();
        if (dbIndexDir != null) {

            File newDBFile = new File(Paths.get(dbIndexDir.getPath(), databaseFile.getName()).toString());
            if (!useTDA) {
                if (!newDBFile.exists()) {
                    System.out.println("Creating " + newDBFile.getPath() + ".");
                    ReverseDB.copyDB(databaseFile.getPath(), newDBFile.getPath());
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
                System.out.println("Creating " + concatTargetDecoyDBFile.getPath() + ".");
                if (ReverseDB.reverseDB(databaseFile.getPath(), concatTargetDecoyDBFile.getPath(), true, DECOY_PROTEIN_PREFIX) == false) {
                    return "Cannot create a decoy database file!";
                }
            }
            databaseFile = concatTargetDecoyDBFile;
        }

        DBScanner.setAminoAcidProbabilities(databaseFile.getPath(), aaSet);
        aaSet.registerEnzyme(enzyme);

        CompactFastaSequence fastaSequence = new CompactFastaSequence(databaseFile.getPath());
        if (useTDA) {
            float ratioUniqueProteins = fastaSequence.getRatioUniqueProteins();
            if (ratioUniqueProteins < 0.5f) {
                fastaSequence.printTooManyDuplicateSequencesMessage(databaseFile.getName(), "MS-GF+");
                System.exit(-1);
            }

            float fractionDecoyProteins = fastaSequence.getFractionDecoyProteins();
            if (fractionDecoyProteins < 0.4f || fractionDecoyProteins > 0.6f) {
                System.err.println("Error while reading: " + databaseFile.getName() + " (fraction of decoy proteins: " + fractionDecoyProteins + ")");
                System.err.println("Delete " + databaseFile.getName() + " and run MS-GF+ again.");
                System.exit(-1);
            }
        }

        CompactSuffixArray sa = new CompactSuffixArray(fastaSequence, params.getMaxPeptideLength());
        System.out.print("Loading database finished ");
        System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - time) / 1000);

        System.out.println("Reading spectra...");

        File specFile = params.getDBSearchIOList().get(ioIndex).getSpecFile();

        SpectraAccessor specAcc = new SpectraAccessor(specFile, specFormat);

        if (specAcc.getSpecMap() == null || specAcc.getSpecItr() == null)
            return "Error while parsing spectrum file: " + specFile.getPath();

        ArrayList<SpecKey> specKeyList = SpecKey.getSpecKeyList(specAcc.getSpecItr(),
                startSpecIndex, endSpecIndex, minCharge, maxCharge, activationMethod, minNumPeaksPerSpectrum);
        int specSize = specKeyList.size();
        if (specSize == 0)
            return specFile.getPath() + " does not have any valid spectra";

        System.out.print("Reading spectra finished ");
        System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - time) / 1000);

        if (numThreads <= 0)
            numThreads = 1;

        int spectraPerTaskMinimum = 250;
        int maxThreads = Math.max(1, Math.round((float) specSize / spectraPerTaskMinimum)); // 250 spectra/task(or thread) minimum for efficiency; going smaller slows down processing
        if (maxThreads < numThreads) {
            System.out.println("Note: " + spectraPerTaskMinimum + " spectra per thread minimum; using " + maxThreads + " threads instead of " + numThreads);
            numThreads = maxThreads;
        }

        System.out.println("Using " + numThreads + (numThreads == 1 ? " thread." : " threads."));

        // Print out parameters
        System.out.println("Search Parameters:");
        System.out.println(params.toString());

        SpecDataType specDataType = new SpecDataType(activationMethod, instType, enzyme, protocol);

        List<MSGFPlusMatch> resultList = Collections.synchronizedList(new ArrayList<MSGFPlusMatch>());

        int toIndexGlobal = specSize;
        while (toIndexGlobal < specSize) {
            SpecKey lastSpecKey = specKeyList.get(toIndexGlobal - 1);
            SpecKey nextSpecKey = specKeyList.get(toIndexGlobal);

            if (lastSpecKey.getSpecIndex() == nextSpecKey.getSpecIndex())
                toIndexGlobal++;
            else
                break;
        }

        System.out.println("Spectrum 0-" + (toIndexGlobal - 1) + " (total: " + specSize + ")");

        // Thread pool
        ThreadPoolExecutorWithExceptions executor = ThreadPoolExecutorWithExceptions.newFixedThreadPool(numThreads);
        executor.setTaskName("Search");

        int numTasks = Math.min(numThreads * 3, Math.round((float) specSize / spectraPerTaskMinimum));
        if (numThreads <= 1) {
            numTasks = 1;
        }

        if (params.getNumTasks() != 0) {
            numTasks = params.getNumTasks();
            if (numTasks < 0) {
                numTasks = numThreads * (numTasks * -1);
            }
            if (numTasks < numThreads) {
                System.out.println("Changing specified tasks from " + numTasks + " to " + numThreads + " to provide the minimum of one task per thread.");
                numTasks = numThreads;
            }
        }
        if (numTasks > 1) {
            System.out.println("Splitting work into " + numTasks + " tasks.");
        } else {
            System.out.println("Searching using a single task.");
        }

        // Partition specKeyList
        int size = toIndexGlobal;
        int residue = size % numTasks;

        int[] startIndex = new int[numTasks];
        int[] endIndex = new int[numTasks];

        int subListSize = size / numTasks;
        for (int i = 0; i < numTasks; i++) {
            startIndex[i] = i > 0 ? endIndex[i - 1] : 0;
            endIndex[i] = startIndex[i] + subListSize + (i < residue ? 1 : 0);

            subListSize = size / numTasks;
            while (endIndex[i] < specKeyList.size()) {
                SpecKey lastSpecKey = specKeyList.get(endIndex[i] - 1);
                SpecKey nextSpecKey = specKeyList.get(endIndex[i]);

                if (lastSpecKey.getSpecIndex() == nextSpecKey.getSpecIndex()) {
                    ++endIndex[i];
                    --subListSize;
                } else
                    break;
            }
        }

        try {
            for (int i = 0; i < numTasks; i++) {
                ScoredSpectraMap specScanner = new ScoredSpectraMap(
                        specAcc,
                        Collections.synchronizedList(specKeyList.subList(startIndex[i], endIndex[i])),
                        leftParentMassTolerance,
                        rightParentMassTolerance,
                        minIsotopeError,
                        maxIsotopeError,
                        specDataType,
                        params.outputAdditionalFeatures(),
                        false
                );
                if (doNotUseEdgeScore)
                    specScanner.turnOffEdgeScoring();

                ConcurrentMSGFPlus.RunMSGFPlus msgfdbExecutor = new ConcurrentMSGFPlus.RunMSGFPlus(
                        specScanner,
                        sa,
                        params,
                        resultList,
                        i + 1
                );
                    executor.execute(msgfdbExecutor);
            }
            // Output initial progress report.
            executor.outputProgressReport();

            executor.shutdown();

            try {
                executor.awaitTerminationWithExceptions(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                if (!executor.HasThrownData()) {
                    e.printStackTrace();
                    Logger.getLogger(MSGFPlus.class.getName()).log(Level.SEVERE, e.getMessage(), e);
                }
            }

            // Output completed progress report.
            executor.outputProgressReport();

        } catch (OutOfMemoryError ex) {
            ex.printStackTrace();
            Logger.getLogger(MSGFPlus.class.getName()).log(Level.SEVERE, null, ex);
            executor.shutdownNow();
            int taskMult = numTasks / numThreads;
            return "Task terminated; results incomplete. Please run again with a greater amount of memory, using \"-Xmx4G\", for example.\n" +
                    "\tYou can also use less memory by increasing the number of tasks used for the search, at the cost of more time.\n" +
                    "\tTry doubling the number used for this search with \"-tasks -" + (taskMult * 2) + "\" or \"-tasks " + (numTasks * 2) + "\".";
        } catch (Exception ex) {
            ex.printStackTrace();
            Logger.getLogger(MSGFPlus.class.getName()).log(Level.SEVERE, null, ex);
            executor.shutdownNow();
            return "Task terminated; results incomplete. Please run again.";
        } catch (Throwable ex) {
            ex.printStackTrace();
            Logger.getLogger(MSGFPlus.class.getName()).log(Level.SEVERE, null, ex);
            executor.shutdownNow();
            return "Task terminated; results incomplete. Please run again.";
        }

        time = System.currentTimeMillis();

        if (params.useTDA()) {
            // Compute Q-values
            System.out.println("Computing q-values...");
            ComputeFDR.addQValues(resultList, sa, false);
            System.out.print("Computing q-values finished ");
            System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - time) / 1000);
        }

        // sort by spectral E-values
        time = System.currentTimeMillis();

        System.out.println("Writing results...");
        Collections.sort(resultList);

        MZIdentMLGen mzidGen = new MZIdentMLGen(params, aaSet, sa, specAcc, ioIndex);
        mzidGen.addSpectrumIdentificationResults(resultList);

        mzidGen.writeResults(outputFile);

        System.out.print("Writing results finished ");
        System.out.format("(elapsed time: %.2f sec)\n", (float) (System.currentTimeMillis() - time) / 1000);

        return null;
    }
}
