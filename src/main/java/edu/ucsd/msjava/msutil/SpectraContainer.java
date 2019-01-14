package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.parser.SpectrumParser;

import java.io.*;
import java.util.ArrayList;


public class SpectraContainer extends ArrayList<Spectrum> {
    /**
     *
     */
    private static final long serialVersionUID = 1L;

    public SpectraContainer() {
    }

    public SpectraContainer(String fileName, SpectrumParser parser) {
        SpectraIterator iterator = null;
        try {
            iterator = new SpectraIterator(fileName, parser);
        } catch (IOException e) {
            e.printStackTrace();
        }
        while (iterator.hasNext())
            this.add(iterator.next());
    }

    public void outputMgfFile(String fileName) {
        PrintStream out = null;
        try {
            out = new PrintStream(new BufferedOutputStream(new FileOutputStream(fileName)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        for (Spectrum spec : this) {
            spec.outputMgf(out);
            out.println();
        }
        out.close();
    }
}
