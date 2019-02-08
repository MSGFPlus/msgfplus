package edu.ucsd.msjava.params;

import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;


public class UserParam {
    public static ArrayList<String> parseFromFile(String fileName, int tokenLength) {
        ArrayList<String> paramLines = new ArrayList<String>();
        BufferedLineReader reader = null;
        try {
            reader = new BufferedLineReader(fileName);
        } catch (IOException e) {
            e.printStackTrace();
        }

        String s;
        while ((s = reader.readLine()) != null) {
            String trimmedLine = s.trim();
            if (trimmedLine.startsWith("#") || trimmedLine.length() == 0) {
                continue;
            }

            String[] token = trimmedLine.split(",");
            if (token.length < tokenLength) {
                continue;
            }

            paramLines.add(trimmedLine);
        }
        return paramLines;
    }

}
