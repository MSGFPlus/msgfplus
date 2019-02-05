package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.ui.MSGFPlus;

import java.io.*;

public class ReverseDB {

    public static void main(String argv[]) {
        if (argv.length != 2)
            printUsageAndExit();

        String ext1 = argv[0].substring(argv[0].lastIndexOf('.') + 1);
        String ext2 = argv[1].substring(argv[1].lastIndexOf('.') + 1);
        if (!ext1.equalsIgnoreCase("fasta") || !ext2.equalsIgnoreCase("fasta")) {
            System.out.println(ext1 + "," + ext2);
            printUsageAndExit();
        }
        String decoyProteinPrefix;
        if (argv.length > 2)
            decoyProteinPrefix = argv[2].trim();
        else
            decoyProteinPrefix = MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX;

        reverseDB(argv[0], argv[1], false, decoyProteinPrefix);

    }

    public static void printUsageAndExit() {
        System.out.println("usage: java ReverseDB input.fasta output.fasta [DecoyProteinPrefix]");
        System.exit(0);
    }

    public static boolean reverseDB(String inFileName, String outFileName, boolean concat, String revPrefix) {
        BufferedReader in = null;
        PrintStream out = null;
        try {
            out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFileName)));
        } catch (FileNotFoundException e1) {
            e1.printStackTrace();
        }

        if (revPrefix == null || revPrefix.trim().isEmpty())
            revPrefix = MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX;

        // Make sure that revPrefix does not end in an underscore, since we add it below
        while (revPrefix.endsWith("_")) {
            revPrefix = revPrefix.substring(0, revPrefix.length() - 1);
        }

        if (revPrefix.trim().isEmpty())
            revPrefix = MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX;

        String s;
        if (concat) {
            try {
                in = new BufferedReader(new FileReader(inFileName));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
            try {
                while ((s = in.readLine()) != null) {
                    out.println(s);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        try {
            in = new BufferedReader(new FileReader(inFileName));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        StringBuffer protein = null;
        String annotation = null;
        try {
            while ((s = in.readLine()) != null) {
                if (s.startsWith(">"))    // start of a protein
                {
                    if (annotation != null) {
                        StringBuffer rev = new StringBuffer();
                        for (int i = protein.length() - 1; i >= 0; i--)
                            rev.append(protein.charAt(i));
                        out.println(">" + revPrefix + "_" + annotation);
                        out.println(rev.toString().trim());
                    }
                    annotation = s.substring(1);
                    protein = new StringBuffer();
                } else
                    protein.append(s);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        if (protein != null && annotation != null) {
            out.println(">" + revPrefix + "_" + annotation);
            out.println(protein.reverse().toString().trim());
        }
        try {
            in.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        out.close();

        return true;
    }

    public static boolean copyDB(String inFileName, String outFileName) {
        BufferedReader in = null;
        PrintStream out = null;
        try {
            out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFileName)));
        } catch (FileNotFoundException e1) {
            e1.printStackTrace();
            return false;
        }

        String s;
        try {
            in = new BufferedReader(new FileReader(inFileName));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return false;
        }

        try {
            while ((s = in.readLine()) != null) {
                out.println(s);
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }

        out.flush();
        out.close();

        return true;
    }
}
