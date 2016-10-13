package edu.ucsd.msjava.misc;

import java.io.File;

public class CountSequestIDs {
    public static void main(String argv[]) throws Exception {
        if (argv.length != 1)
            printUsageAndExit("Invalid parameters!");

        File seqDir = new File(argv[0]);
        if (!seqDir.isDirectory())
            printUsageAndExit(argv[0] + " is not a directory!");
    }

    public static void printUsageAndExit(String message) {
        if (message != null)
            System.out.println(message);
        System.exit(-1);
    }

    public static void processPPResults(File seqDir) throws Exception {
        File synFile = null;
        File synPPFile = null;

        for (File f : seqDir.listFiles()) {
            String fileName = f.getName();
            if (fileName.endsWith("_syn.txt") || fileName.endsWith("_syn.tsv"))
                synFile = f;
            else if (fileName.endsWith("_syn_PepProphet.txt") || fileName.endsWith("_syn_PepProphet.tsv"))
                synPPFile = f;
        }

        if (synFile == null)
            printUsageAndExit("_syn.txt file is missing!");
        if (synPPFile == null)
            printUsageAndExit("_syn_PepProphet.txt file is missing!");


    }
}
