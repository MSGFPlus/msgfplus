package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;

import edu.ucsd.msjava.parser.BufferedLineReader;

public class GetProteinLength {
    public static void main(String argv[]) throws Exception {
        File tsvFile = new File("/Users/kims336/Research/Data/SNU/SpecCounts.tsv");
        File outputFile = new File("/Users/kims336/Research/Data/SNU/QSpecInput.tsv");
        File dbFile = new File("/Users/kims336/Research/Data/CommonContaminants/H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta");
        convert(tsvFile, outputFile, dbFile);
    }

    public static void convert(File tsvFile, File outputFile, File dbFile) throws Exception {
        BufferedLineReader in = new BufferedLineReader(tsvFile.getPath());
        PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
        HashMap<String, Integer> protLengthMap = MSGFDBToQSpec.getAnnotationProtLengthMap(dbFile.getPath());

        String s;
        in.readLine();    // header
        out.println("protid\tprotLen\t0\t1");
        while ((s = in.readLine()) != null) {
            String[] token = s.split("\t");
            if (token.length != 3)
                continue;
            String protId = token[0];
            Integer length = protLengthMap.get(protId);
            if (length == null) {
                System.out.println(protId + " doesn't exist in the database!");
                System.exit(-1);
            }
            out.println(protId + "\t" + length + "\t" + token[1] + "\t" + token[2]);
        }

        out.close();
        in.close();
        System.out.println("Done");
    }
}
