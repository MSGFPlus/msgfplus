package edu.ucsd.msjava.ims;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.parser.BufferedLineReader;

public class Summarize {
    public static void main(String argv[]) throws Exception {
        summarize();
    }

    public static void summarize() throws Exception {
        File dir = new File("/Users/kims336/Research/Data/IMS/5milTargetsRev/");
        File resultFile = new File(dir.getPath() + File.separator + "bestPerPeptides_Rev.tsv");
        String s;
        BufferedLineReader in = new BufferedLineReader(resultFile.getPath());

        String headerLine = in.readLine();
        String[] header = headerLine.split("\t");
        int specProbCol = -1;
        for (int i = 0; i < header.length; i++) {
            if (header[i].equals("SpecProb"))
                specProbCol = i;
        }
        if (specProbCol < 0) {
            System.out.println("No SpecProb column");
            System.exit(-1);
        }

        PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(dir.getPath() + File.separator + "goodPeptides.tsv")));
        out.println(headerLine);
        float[] threshold = {1e-10f, 1e-11f, 1e-12f, 1e-13f};
        int[] numID = new int[threshold.length];
        while ((s = in.readLine()) != null) {
            if (s.startsWith("#") || s.length() == 0)
                continue;
            String[] token = s.split("\t");
            if (token.length != header.length)
                continue;
            float specProb = Float.parseFloat(token[specProbCol]);
            if (specProb < threshold[0]) {
                out.print(token[0]);
                for (int i = 1; i < token.length; i++) {
                    if (i == 2)
                        out.print("\t" + token[i].replaceAll("\\+15\\.995", "@").replaceAll("C", "C!"));
                    else
                        out.print("\t" + token[i]);
                }
                out.println();
            }

            for (int i = 0; i < threshold.length; i++) {
                if (specProb < threshold[i]) {
                    numID[i]++;
                }
            }
        }

        System.out.println("Threshold\tNumID");
        for (int i = 0; i < threshold.length; i++) {
            System.out.println(threshold[i] + "\t" + numID[i]);
        }
        in.close();
        out.close();
    }
}
