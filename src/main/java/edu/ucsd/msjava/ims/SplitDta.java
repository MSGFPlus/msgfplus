package edu.ucsd.msjava.ims;

import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

public class SplitDta {
    public static void main(String argv[]) throws Exception {
        if (argv.length != 1)
            printUsageAndExit("Invalid parameter");

        File dtaFile = new File(argv[0]);
        if (!dtaFile.exists())
            printUsageAndExit("File does not exist.");

        if (dtaFile.isDirectory())
            printUsageAndExit(dtaFile.getName() + " is a directory!");

        if (!dtaFile.getName().endsWith("_dta.txt"))
            printUsageAndExit("Wrong spectrum file format");

        split(dtaFile);
    }

    public static void printUsageAndExit(String message) {
        if (message != null)
            System.out.println(message);
        System.out.println("Usage: java DtaToMSGFInputDB DTAPath(*_dta.txt)");
        System.exit(-1);
    }

    public static void split(File dtaFile) throws Exception {
        int specIndex = 0;
        String s;

        BufferedLineReader in = new BufferedLineReader(dtaFile.getPath());
        PrintStream out = null;

        String filePath = dtaFile.getPath();
        String prefix = filePath.substring(0, filePath.lastIndexOf("_dta.txt"));
        int fileNum = 0;

        while ((s = in.readLine()) != null) {
            if (s.startsWith("===")) {
                if (specIndex % 1000000 == 0) {
                    if (out != null)
                        out.close();
                    fileNum++;
                    out = new PrintStream(new BufferedOutputStream(new FileOutputStream(prefix + fileNum + "_dta.txt")));
                }
                specIndex++;
            }
            out.println(s);
        }

        if (out != null)
            out.close();
        in.close();
    }
}
