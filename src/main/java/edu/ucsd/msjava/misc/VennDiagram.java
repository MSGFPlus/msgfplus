package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.parser.TSVResultParser;

import java.io.File;
import java.util.Set;

public class VennDiagram {
    public static void main(String argv[]) throws Exception {
        if (argv.length != 2 && argv.length != 3)
            printUsageAndExit();
        File resultFile1 = new File(argv[0]);
        File resultFile2 = new File(argv[1]);
        float threshold = 0.01f;
        if (argv.length == 3)
            threshold = Float.parseFloat(argv[2]);
        vennDiagram(resultFile1, resultFile2, threshold);
    }

    public static void printUsageAndExit() {
        System.out.println("usage: java VennDiagram result1 result2 [threshold]");
        System.exit(-1);
    }

    public static void vennDiagram(File result1, File result2, float fdrThreshold) throws Exception {
        TSVResultParser parser1 = new TSVResultParser(result1);
        String err;
        if ((err = parser1.parse(fdrThreshold)) != null) {
            System.out.println(err);
            System.exit(-1);
        }

        TSVResultParser parser2 = new TSVResultParser(result2);
        if ((err = parser2.parse(fdrThreshold)) != null) {
            System.out.println(err);
            System.exit(-1);
        }

        int[] vennScan = getVennDiagram(parser1.getScanSet(), parser2.getScanSet());
        System.out.println("PSM: [ " + vennScan[0] + " [" + vennScan[2] + "] " + vennScan[1] + "]");
        int[] vennPep = getVennDiagram(parser1.getPepSet(), parser2.getPepSet());
        System.out.println("Peptide: [" + vennPep[0] + " [" + vennPep[2] + "] " + vennPep[1] + "]");
    }

    public static int[] getVennDiagram(Set<String> set1, Set<String> set2) {
        int shared = 0;
        for (String item1 : set1) {
            if (set2.contains(item1))
                shared++;
        }

        return new int[]{set1.size() - shared, set2.size() - shared, shared};
    }
}
