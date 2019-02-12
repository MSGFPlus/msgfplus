package edu.ucsd.msjava.ui;

import edu.ucsd.msjava.msdbsearch.ConcurrentMSGFDB;
import edu.ucsd.msjava.msdbsearch.ScoredSpectraMap;
import edu.ucsd.msjava.msgf.MSGFDBResultGenerator;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.params.ToleranceParameter;
import edu.ucsd.msjava.sequences.Constants;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

@Deprecated
public class MSGFDBLib {
    public static final String VERSION = "7573";
    public static final String RELEASE_DATE = "04/04/2012";

    public static final String DECOY_DB_EXTENSION = ".revConcat.fasta";

    public static void main(String argv[]) {
        long time = System.currentTimeMillis();

        ParamManager paramManager = new ParamManager("MSGFDBLib", VERSION, RELEASE_DATE, "java -Xmx2000M -jar MSGFDBLib.jar");
        paramManager.addMSGFLibParams();

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

        // Running MS-GFDB
        paramManager.printToolInfo();
        String errorMessage = runMSGFLib(paramManager);
        if (errorMessage != null) {
            System.err.println("[Error] " + errorMessage);
            System.out.println();
        } else
            System.out.format("MS-GFDBLib complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis() - time) / (float) 1000);
    }

    public static String runMSGFLib(ParamManager paramManager) {
        long time = System.currentTimeMillis();

        // Spectrum file
        FileParameter specParam = paramManager.getSpecFileParam();
        File specFile = specParam.getFile();
        SpecFileFormat specFormat = (SpecFileFormat) specParam.getFileFormat();

        // Library file
        File libraryFile = paramManager.getFile("d");

        // PM tolerance
        ToleranceParameter tol = ((ToleranceParameter) paramManager.getParameter("t"));
        Tolerance leftParentMassTolerance = tol.getLeftTolerance();
        Tolerance rightParentMassTolerance = tol.getRightTolerance();

        int numAllowedC13 = paramManager.getIntValue("c13");
        if (rightParentMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f)
            numAllowedC13 = 0;

        File outputFile = paramManager.getOutputFileParam().getFile();

        Enzyme enzyme = paramManager.getEnzyme();
        ActivationMethod activationMethod = paramManager.getActivationMethod();
        InstrumentType instType = paramManager.getInstType();
        if (activationMethod == ActivationMethod.HCD)
            instType = InstrumentType.HIGH_RESOLUTION_LTQ;

        Protocol protocol = paramManager.getProtocol();

        int numMatchesPerSpec = paramManager.getIntValue("n");
        int numThreads = paramManager.getIntValue("thread");

//		System.out.println("Loading the library file...");
//		System.out.print("Loading database finished ");
//		System.out.format("(elapsed time: %.2f sec)\n", (float)(System.currentTimeMillis()-time)/1000);

        System.out.println("Reading spectra...");
        SpectraAccessor specAcc = new SpectraAccessor(specFile, specFormat);

        if (specAcc.getSpecMap() == null || specAcc.getSpecItr() == null)
            return "Error while parsing spectrum file: " + specFile.getPath();

        // determine the number of spectra to be scanned together
        long maxMemory = Runtime.getRuntime().maxMemory() - 1 << 28;

        int avgPeptideMass = 2000;
        int numBytesPerMass = 12;
        int numSpecScannedTogether = (int) ((float) maxMemory / avgPeptideMass / numBytesPerMass);
        ArrayList<SpecKey> specKeyList = SpecKey.getSpecKeyList(specAcc.getSpecItr(), 0, Integer.MAX_VALUE, 0, Integer.MAX_VALUE, activationMethod, Constants.MIN_NUM_PEAKS_PER_SPECTRUM);
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
                        leftParentMassTolerance,
                        rightParentMassTolerance,
                        numAllowedC13,
                        specDataType,
                        false
                );

                ConcurrentMSGFDB.RunMSGFDBLib msgfdbExecutor = new ConcurrentMSGFDB.RunMSGFDBLib(
                        specScanner,
                        numMatchesPerSpec,
                        resultList,
                        specFile.getName(),
                        libraryFile.getPath()
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
                        + (rightParentMassTolerance.isTolerancePPM() ? "ppm" : "Da")
                        + ")\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecProb\tP-value";

        MSGFDBResultGenerator gen = new MSGFDBResultGenerator(header, resultList);
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

        return null;
    }
}
