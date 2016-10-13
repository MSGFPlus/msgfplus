package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.msscorer.NewRankScorer;

import java.io.File;

public class ParamToTxt {
    public static void main(String argv[]) throws Exception {
        if (argv.length != 1 || !new File(argv[0]).exists()) {
            pringUsageAndExit(null);
        }

        File paramFile = new File(argv[0]);
        String ext = paramFile.getName().substring(paramFile.getName().lastIndexOf('.'));
        if (!ext.equalsIgnoreCase(".param"))
            pringUsageAndExit("Invalid file format: " + paramFile.getAbsolutePath());

        convert(paramFile);
    }

    public static void pringUsageAndExit(String message) {
        if (message != null)
            System.out.println(message);
        System.out.println("java ParamToTxt *.param");
        System.exit(-1);
    }

    public static void convert(File paramFile) {
        String filePath = paramFile.getAbsolutePath();
        File outputFile = new File(filePath.substring(0, filePath.lastIndexOf('.')) + ".txt");
        NewRankScorer scorer = new NewRankScorer(paramFile.getPath());
        scorer.writeParametersPlainText(outputFile);
    }
}
