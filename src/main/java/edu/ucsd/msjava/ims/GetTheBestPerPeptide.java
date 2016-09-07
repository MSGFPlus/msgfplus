package edu.ucsd.msjava.ims;

import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class GetTheBestPerPeptide {
    public static void main(String argv[]) throws Exception {
        if (argv.length != 1)
            printUsageAndExit("Illegal parameter.");
        File tsvFile = new File(argv[0]);
        if (!tsvFile.exists())
            printUsageAndExit("File does not exist.");
        getTheBest(tsvFile);
    }

    public static void printUsageAndExit(String message) {
        if (message != null)
            System.out.println(message);
        System.out.println("Usage: java GetTheBestPerScan TSVFile");
        System.exit(-1);
    }

    public static void getTheBest(File tsvFile) throws Exception {
        BufferedLineReader in = new BufferedLineReader(tsvFile.getPath());
        String header = in.readLine();
        String[] headerToken = header.split("\t");

        int annotationIndexCol = -1;
        int specProbCol = -1;
        for (int i = 0; i < headerToken.length; i++) {
            if (headerToken[i].equalsIgnoreCase("Annotation"))
                annotationIndexCol = i;
            if (headerToken[i].equalsIgnoreCase("SpecProb"))
                specProbCol = i;
        }

        if (annotationIndexCol == -1) {
            System.err.println("DtaIndex column does not exist.");
            System.exit(-1);
        }
        if (specProbCol == -1) {
            System.err.println("SpecProb column does not exist.");
            System.exit(-1);
        }

        Map<String, String> table = new LinkedHashMap<String, String>();

        String s;
        while ((s = in.readLine()) != null) {
            if (s.startsWith("#"))
                continue;
            String[] token = s.split("\t");
            String annotation = token[annotationIndexCol];
            String prev = table.get(annotation);
            if (prev == null)
                table.put(annotation, s);
            else {
                String[] tokenPrev = prev.split("\t");
                float prevSpecProb = Float.parseFloat(tokenPrev[specProbCol]);
                float specProb = Float.parseFloat(token[specProbCol]);
                if (specProb < prevSpecProb)
                    table.put(annotation, s);
            }
        }

        in.close();

        List<String> annotationList = new ArrayList<String>(table.keySet());
        System.out.println(header);
        for (String annotation : annotationList)
            System.out.println(table.get(annotation));
    }
}
