package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;

public class IPA {
    public static void main(String argv[]) throws Exception {
        testIPA();
        System.out.println("Done");
    }

    public static void testIPA() throws Exception {
        String iqResultFile = "C:\\cygwin\\home\\kims336\\Research\\Data\\IQ\\QCShew\\Gordon\\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11_results.tsv";
        BufferedLineReader in2 = new BufferedLineReader(iqResultFile);
        String[] header2 = in2.readLine().split("\t");

        HashMap<String, Float> pepFitMap = new HashMap<String, Float>();
        String s2;
        while ((s2 = in2.readLine()) != null) {
            String[] token = s2.split("\t");
            if (token.length < 20)
                continue;
            int targetID = Integer.parseInt(token[1]);
            String peptide = token[2];
            float fit = Float.parseFloat(token[18]);
            float iScore = Float.parseFloat(token[19]);

            Float prevFit = pepFitMap.get(peptide);
            if (prevFit == null || prevFit > fit)
                pepFitMap.put(peptide, fit);
        }
        in2.close();

        String outputFile = "C:\\cygwin\\home\\kims336\\Research\\Data\\IQ\\QCShew\\Gordon\\IPA.tsv";
        PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

        String msgfResultFile = "C:\\cygwin\\home\\kims336\\Research\\Data\\IQ\\QCShew\\Target_TD.txt";
//		HashMap<String,Float> msgfMap = new HashMap<String,Float>();
        BufferedLineReader in = new BufferedLineReader(msgfResultFile);
        String header = in.readLine();
        out.println(header + "\tFit");
        String s;
        while ((s = in.readLine()) != null) {
            String[] token = s.split("\t");
            if (token.length < 16)
                continue;
            String peptide = token[8];
            float fit = pepFitMap.get(peptide);
            out.println(s + "\t" + fit);
        }
        in.close();
        out.close();
    }

}
