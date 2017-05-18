package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.misc.ExceptionCapturer;
import edu.ucsd.msjava.misc.ProgressData;
import edu.ucsd.msjava.misc.ProgressReporter;
import edu.ucsd.msjava.misc.ThreadPoolExecutorWithExceptions;
import edu.ucsd.msjava.mzid.MzIDParser;
import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.TimeUnit;

public class AnnotatedSpectra implements ProgressReporter, ExceptionCapturer {
    private File[] resultFiles;
    private File specDir;
    private AminoAcidSet aaSet;
    private float fdrThreshold = 0.01f;
    private ProgressData progress;
    private boolean dropErrorDatasets = false;
    private Throwable exception = null;

    @Override
    public void setProgressData(ProgressData data) {
        progress = data;
    }

    @Override
    public ProgressData getProgressData() {
        return progress;
    }
    
    @Override
    public boolean hasException() {
        return exception != null;
    }
    
    @Override
    public Throwable getException() {
        return exception;
    }
    
    public void setDropErrorDatasets(boolean dropErrors) {
        dropErrorDatasets = dropErrors;
    }

    private SpectraContainer annotatedSpectra;

    public AnnotatedSpectra(File[] resultFiles, File specDir, AminoAcidSet aaSet) {
        this.resultFiles = resultFiles;
        this.specDir = specDir;
        this.aaSet = aaSet;
        this.progress = null;
    }
    
    public AnnotatedSpectra(File[] resultFiles, File specDir, AminoAcidSet aaSet, float fdrThreshold, boolean dropErrors) {
        this.resultFiles = resultFiles;
        this.specDir = specDir;
        this.aaSet = aaSet;
        this.fdrThreshold = fdrThreshold;
        this.dropErrorDatasets = dropErrors;
        this.progress = null;
    }
    
    public static class ConcurrentAnnotatedSpectraParser extends AnnotatedSpectra implements Runnable {
        private List<Spectrum> results;
        private List<String> errors;
        
        public ConcurrentAnnotatedSpectraParser(File[] resultFiles, File specDir, AminoAcidSet aaSet, float fdrThreshold, boolean dropErrors, List<Spectrum> resultList, List<String> errorList) {
            super(resultFiles, specDir, aaSet, fdrThreshold, dropErrors);
            results = resultList;
            errors = errorList;
        }
        
        @Override
        public void run() {
            //String result = parse();
            String result = parse();
            results.addAll(getAnnotatedSpecContainer());
            if (result != null) {
                errors.add(result);
                //System.out.println("ERROR: " + result);
            }
        }
    }

    public AnnotatedSpectra fdrThreshold(float fdrThreshold) {
        this.fdrThreshold = fdrThreshold;
        return this;
    }

    public SpectraContainer getAnnotatedSpecContainer() {
        return annotatedSpectra;
    }
    
    public String parse(int numThreads, boolean dropErrors) {
        if (numThreads <= 1) {
            return parse();
        }

        List<Spectrum> results = Collections.synchronizedList(new ArrayList<Spectrum>());
        List<String> errors = Collections.synchronizedList(new ArrayList<String>());

        // Thread pool
        ThreadPoolExecutorWithExceptions executor = ThreadPoolExecutorWithExceptions.newFixedThreadPool(numThreads);
        executor.setTaskName("Parse");

        try {
            List<List<File>> taskFiles = new ArrayList<List<File>>();
            for (int i = 0; i < numThreads; i++) {
                taskFiles.add(new ArrayList<File>());
            }
            for (int i = 0; i < resultFiles.length; i++)
            {
                // Evenly distribute the files to the threads; doing it in a collated fashion because files
                // of similar size will often have similar names, and we don't want to give one thread all
                // of the largest files.
                taskFiles.get(i % numThreads).add(resultFiles[i]);
            }
            for (int i = 0; i < numThreads; i++) {
                List<File> thisTaskFiles = taskFiles.get(i);
                System.out.println("Task " + (i + 1) + ": " + thisTaskFiles.size() + " files.");
                ConcurrentAnnotatedSpectraParser parser = new ConcurrentAnnotatedSpectraParser(thisTaskFiles.toArray(new File[0]), specDir, aaSet, fdrThreshold, dropErrors, results, errors);
                executor.execute(parser);
            }
            taskFiles.clear();
            
            // Output initial progress report.
            executor.outputProgressReport();

            executor.shutdown();

            try {
                executor.awaitTerminationWithExceptions(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                if (!executor.HasThrownData()) {
                    e.printStackTrace();
                }
            }
            
            // Output completed progress report.
            executor.outputProgressReport();
        } catch (OutOfMemoryError ex) {
            ex.printStackTrace();
            executor.shutdownNow();
            return "Task terminated; results incomplete. Please run again with a greater amount of memory, using \"-Xmx4G\", for example.";
        } catch (Exception ex) {
            ex.printStackTrace();
            executor.shutdownNow();
            return "Task terminated; results incomplete. Please run again.";
        } catch (Throwable ex) {
            ex.printStackTrace();
            executor.shutdownNow();
            return "Task terminated; results incomplete. Please run again.";
        }
        annotatedSpectra = new SpectraContainer();
        annotatedSpectra.addAll(results);
        String errorList = null;
        for (String error : errors) {
            if (errorList == null) {
                errorList = error;
            } else {
                errorList += "\n" + error;
            }
        }
        return errorList;
    }

    public String parse() {
        if (progress == null) {
            progress = new ProgressData();
        }
        annotatedSpectra = new SpectraContainer();

        System.out.println("Using " + resultFiles.length + " result files:");
        for (File resultFile : resultFiles)
            System.out.println("\t" + resultFile.getName());

        
        int count = 0;
        int total = resultFiles.length;
        String aggErrs = null;
        for (File resultFile : resultFiles) {
            String errMsg = parseFile(resultFile);
            count++;
            progress.report(count, total);
            if (errMsg != null){
                String msg = "Error while parsing " + resultFile.getName() + ": " + errMsg;
                if (dropErrorDatasets) {
                    System.out.println(msg);
                    if (aggErrs == null) {
                        aggErrs = msg;
                    } else {
                        aggErrs += "\n" + msg;
                    }
                } else {
                    exception = new Exception(msg);
                    return msg;
                }
            }
        }
        return null;
    }

    public void writeToMgf(PrintStream out) {
        if (annotatedSpectra != null) {
            for (Spectrum spec : annotatedSpectra)
                spec.outputMgf(out);
        }
    }

    public String parseFile(File resultFile) {
        System.out.println("Parsing " + resultFile.getName());
        File tsvResultFile = null;
        if (resultFile.getName().endsWith(".mzid")) {
            String resultFileName = resultFile.getAbsolutePath();
            String tsvResultFileName = resultFileName.substring(0, resultFileName.lastIndexOf('.')) + ".tsv";
            tsvResultFile = new File(tsvResultFileName);
            if (!tsvResultFile.exists()) {
                if (!tsvResultFile.canWrite()) {
                    MzIDParser parser = new MzIDParser(resultFile);
                    parser.writeToTSVFile(tsvResultFile);
                } else {
                    try {
                        System.out.println("Converting " + resultFile.getName());
                        tsvResultFile = File.createTempFile("__AnnotatedSpectra", ".tsv");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                    tsvResultFile.deleteOnExit();

                    MzIDParser parser = new MzIDParser(resultFile);
                    parser.writeToTSVFile(tsvResultFile);
                }
            } else {
                System.out.println(tsvResultFileName + " already exists.");
            }
        } else if (resultFile.getName().endsWith(".tsv")) {
            tsvResultFile = resultFile;
        }


        BufferedLineReader in = null;
        try {
            in = new BufferedLineReader(tsvResultFile.getPath());
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        String s = in.readLine();

        if (!s.startsWith("#")) {
            return "Not a valid tsv result file";
        }

        int specIdCol = -1;
        int specFileCol = -1;
        int pepCol = -1;
        int fdrCol = -1;
        int chargeCol = -1;

        String[] label = s.split("\t");
        for (int i = 0; i < label.length; i++) {
            if (label[i].equalsIgnoreCase("#SpecFile"))
                specFileCol = i;
            else if (label[i].equalsIgnoreCase("SpecID"))
                specIdCol = i;
            else if (label[i].equalsIgnoreCase("Peptide"))
                pepCol = i;
            else if (label[i].equalsIgnoreCase("FDR") || label[i].equalsIgnoreCase("EFDR") || label[i].equalsIgnoreCase("QValue") || label[i].equalsIgnoreCase("SpecQValue"))
                fdrCol = i;
            else if (label[i].equalsIgnoreCase("Charge"))
                chargeCol = i;
        }
        if (specIdCol < 0 || specFileCol < 0 || pepCol < 0)
            return "Not a valid mzid file";
        if (fdrCol < 0)
            return "QValue is missing";

        ArrayList<String> resultList = new ArrayList<String>();
        while ((s = in.readLine()) != null) {
            String[] token = s.split("\t");
            if (token.length <= specIdCol || token.length <= specFileCol || token.length <= pepCol || token.length <= fdrCol)
                continue;

            float fdr = Float.parseFloat(token[fdrCol]);

            if (fdr <= fdrThreshold) {
                resultList.add(s);
            }
        }

        Iterator<String> itr = resultList.iterator();
        List<Spectrum> annotatedResults = new ArrayList<Spectrum>();

        HashMap<String, SpectraAccessor> specAccessorMap = new HashMap<String, SpectraAccessor>();
        while (itr.hasNext()) {
            String str = itr.next();
            String[] token = str.split("\t");

            String pep = token[pepCol];
            if (pep.matches(".\\..+\\.."))
                pep = pep.substring(pep.indexOf('.') + 1, pep.lastIndexOf('.'));

            String specFileName = token[specFileCol];
            specFileName = new File(specFileName).getName();

            int charge = Integer.parseInt(token[chargeCol]);

            SpectraAccessor specAccessor = specAccessorMap.get(specFileName);
            if (specAccessor == null) {
                File specFile = new File(specDir.getPath() + File.separator + specFileName);
                specAccessor = new SpectraAccessor(specFile);
                specAccessorMap.put(specFileName, specAccessor);
            }

            String specId = token[specIdCol];
            Spectrum spec = specAccessor.getSpectrumById(specId);

            if (spec == null)
                return specFileName + ":" + specId + " is not available!";
            else {
                Peptide peptide = new Peptide(pep, aaSet);
                spec.setCharge(charge);

                if (Math.abs(spec.getPeptideMass() - peptide.getMass()) < 5) {
                    spec.setAnnotation(peptide);
                    annotatedResults.add(spec);
                } else {
                    return "parent mass doesn't match " + specFileName + ":" + specId + " " + peptide.toString() + " " + spec.getPeptideMass() + " != " + peptide.getMass();
                }
            }
        }

        try {
            in.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        annotatedSpectra.addAll(annotatedResults);
        return null;
    }
}
