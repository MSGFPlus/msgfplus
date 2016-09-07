package edu.ucsd.msjava.msutil;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import edu.ucsd.msjava.mzid.MzIDParser;
import edu.ucsd.msjava.parser.BufferedLineReader;

public class AnnotatedSpectra {
    private File[] resultFiles;
    private File specDir;
    private AminoAcidSet aaSet;
    private float fdrThreshold = 0.01f;

    private SpectraContainer annotatedSpectra;

    public AnnotatedSpectra(File[] resultFiles, File specDir, AminoAcidSet aaSet) {
        this.resultFiles = resultFiles;
        this.specDir = specDir;
        this.aaSet = aaSet;
    }

    public AnnotatedSpectra fdrThreshold(float fdrThreshold) {
        this.fdrThreshold = fdrThreshold;
        return this;
    }

    public SpectraContainer getAnnotatedSpecContainer() {
        return annotatedSpectra;
    }

    public String parse() {
        annotatedSpectra = new SpectraContainer();

        System.out.println("Using " + resultFiles.length + " result files:");
        for (File resultFile : resultFiles)
            System.out.println("\t" + resultFile.getName());

        for (File resultFile : resultFiles) {
            String errMsg = parseFile(resultFile);
            if (errMsg != null)
                return "Error while parsing " + resultFile.getName() + ": " + errMsg;
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
                try {
                    System.out.println("Converting " + resultFile.getName());
                    tsvResultFile = File.createTempFile("__AnnotatedSpectra", ".tsv");
                } catch (IOException e) {
                    e.printStackTrace();
                }
                tsvResultFile.deleteOnExit();

                MzIDParser parser = new MzIDParser(resultFile);
                parser.writeToTSVFile(tsvResultFile);
            } else {
                System.out.println(tsvResultFileName + " already exists.");
            }
        } else if (resultFile.getName().endsWith(".tsv"))
            tsvResultFile = resultFile;


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
                    annotatedSpectra.add(spec);
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
        return null;
    }
}
