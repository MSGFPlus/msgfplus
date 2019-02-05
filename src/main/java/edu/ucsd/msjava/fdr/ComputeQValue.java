package edu.ucsd.msjava.fdr;

import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.File;
import java.util.ArrayList;

public class ComputeQValue {
    public static final float FDR_REPORT_THRESHOLD = 0.1f;

    public static void main(String argv[]) throws Exception {
        // required
        File targetFile = null;

        // optional
        File outputFile = null;
        boolean isConcatenated = false;
        boolean includeDecoy = false;

        float fdrThreshold = 1;
        float pepFDRThreshold = 1;

        int i = 0;
        while (i < argv.length) {
            // 	-f resultFileName dbCol decoyPrefix or -f targetFileName decoyFileName
            if (argv[i].equalsIgnoreCase("-f")) {
                if (i + 1 >= argv.length)
                    printUsageAndExit("Invalid parameter: " + argv[i]);
                targetFile = new File(argv[i + 1]);
                if (!targetFile.exists())
                    printUsageAndExit(argv[i + 1] + " doesn't exist.");
                else if (!targetFile.isFile())
                    printUsageAndExit(argv[i + 1] + " is not a file.");
                i += 2;
            } else if (argv[i].equalsIgnoreCase("-o")) {
                if (i + 1 >= argv.length)
                    printUsageAndExit("Invalid parameter: " + argv[i]);
                outputFile = new File(argv[i + 1]);
                i += 2;
            } else if (argv[i].equalsIgnoreCase("-decoy")) {
                if (argv[i + 1].equalsIgnoreCase("1"))
                    includeDecoy = true;
                i += 2;
            } else if (argv[i].equalsIgnoreCase("-fdr")) {
                if (i + 1 >= argv.length)
                    printUsageAndExit("Invalid parameter: " + argv[i]);
                try {
                    fdrThreshold = Float.parseFloat(argv[i + 1]);
                } catch (NumberFormatException e) {
                    printUsageAndExit("Invalid pepCol: " + argv[i + 1]);
                }
                i += 2;
            } else if (argv[i].equalsIgnoreCase("-pepfdr")) {
                if (i + 1 >= argv.length)
                    printUsageAndExit("Invalid parameter: " + argv[i]);
                try {
                    pepFDRThreshold = Float.parseFloat(argv[i + 1]);
                } catch (NumberFormatException e) {
                    printUsageAndExit("Invalid pepCol: " + argv[i + 1]);
                }
                i += 2;
            } else {
                printUsageAndExit("Invalid parameter");
            }
        }

        if (targetFile == null)
            printUsageAndExit("Target is missing!");

        computeFDR(targetFile, isConcatenated, includeDecoy, fdrThreshold, pepFDRThreshold, outputFile);
    }

    public static void printUsageAndExit(String message) {
        System.err.println(message);
        System.out.print("Usage: java -cp MSGFPlus.jar fdr.ComputeFDR\n" +
                "\t -f MSGFPlusFileName (*.tsv)\n" +
                "\t [-o outputFileName (default: stdout)]\n" +
                "\t [-fdr fdrThreshold]\n" +
                "\t [-decoy 0/1 (0: don't include decoy (default), 1: include decoy)\n"
                "\t [-pepfdr pepFDRThreshold]\n" +
        );
        System.exit(-1);
    }

    public static void computeFDR(File msgfTsvFile, boolean isConcatenated, boolean includeDecoy,
                                  float fdrThreshold, float pepFDRThreshold, File outputFile) throws Exception {
        // const
        boolean isGreaterBetter = false;
        boolean hasHeader = true;
        File decoyFile = null;
        String delimiter = "\t";
        ArrayList<Pair<Integer, ArrayList<String>>> reqStrList = new ArrayList<Pair<Integer, ArrayList<String>>>();
        String decoyPrefix = "XXX";

        int scoreCol = -1;
        int specFileCol = -1;
        int pepCol = -1;
        int specIndexCol = -1;
        int dbCol = -1;

        BufferedLineReader in = new BufferedLineReader(msgfTsvFile.getPath());
        String header = in.readLine();
        if (header == null) // || (!header.startsWith("#") && !header.startsWith("PSMId")))
        {
            System.out.println("Not a valid MS-GF+ result file!");
            System.exit(0);
        }
        String[] headerToken = header.split("\t");
        for (int i = 0; i < headerToken.length; i++) {
            if (headerToken[i].equalsIgnoreCase("SpecEValue"))
                scoreCol = i;
            if (headerToken[i].equalsIgnoreCase("#SpecFile"))
                specFileCol = i;
            if (headerToken[i].equalsIgnoreCase("Peptide"))
                pepCol = i;
            if (headerToken[i].equalsIgnoreCase("SpecID"))
                specIndexCol = i;
            if (headerToken[i].equalsIgnoreCase("Protein"))
                dbCol = i;
        }

        if (scoreCol < 0) {
            System.out.println("SpecEValue column is missing!");
            System.exit(-1);
        }
        if (specFileCol < 0) {
            System.out.println("SpecFile column is missing!");
            System.exit(-1);
        }
        if (pepCol < 0) {
            System.out.println("Peptide column is missing!");
            System.exit(-1);
        }
        if (specIndexCol < 0) {
            System.out.println("SpecID column is missing!");
            System.exit(-1);
        }
        if (dbCol < 0) {
            System.out.println("Protein column is missing!");
            System.exit(-1);
        }

        ComputeFDR.computeFDR(msgfTsvFile, decoyFile,
                scoreCol, isGreaterBetter,
                delimiter, specFileCol, specIndexCol, pepCol, reqStrList,
                isConcatenated, includeDecoy, hasHeader,
                dbCol, decoyProteinPrefix, fdrThreshold, pepFDRThreshold, outputFile);
    }
}
