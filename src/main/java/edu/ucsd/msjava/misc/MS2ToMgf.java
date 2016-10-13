package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class MS2ToMgf {
    public static void main(String argv[]) throws Exception {
        if (argv.length != 1)
            printUsageAndExit("Invalid parameters");

        File ms2File = new File(argv[0]);
        if (!ms2File.exists() || !ms2File.getName().endsWith(".ms2"))
            printUsageAndExit("Invalid ms2 input: " + argv[0]);

        String ms2FileName = ms2File.getAbsolutePath();
        File mgfFile = new File(ms2FileName.substring(0, ms2FileName.lastIndexOf('.')) + ".mgf");
        convert(ms2File, mgfFile);
    }

    public static void printUsageAndExit(String message) throws Exception {
        System.out.println(message);
        System.out.println("Usage: java MS2ToMgf *.ms2");
        System.exit(-1);
    }

    public static void convert(File ms2File, File mgfFile) throws Exception {
        System.out.println("Converting " + ms2File.getAbsolutePath() + " into " + mgfFile.getAbsolutePath());

        PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(mgfFile)));

        BufferedLineReader in = new BufferedLineReader(ms2File.getPath());
        String s;

        int startScanNum = -1;
        int endScanNum = -1;
        int ms2Index = 0;

        List<Spectrum> specList = null;
        while ((s = in.readLine()) != null) {
            String[] token = s.split("\\s+");
            if (s.startsWith("H") || s.startsWith("I"))
                continue;
            else if (s.startsWith("S")) {
                ms2Index++;
                if (specList != null) {
                    for (Spectrum spec : specList)
                        spec.outputMgf(out);
                    specList.clear();
                }

                startScanNum = Integer.parseInt(token[1]);
                endScanNum = Integer.parseInt(token[2]);
                specList = new ArrayList<Spectrum>();
            } else if (s.startsWith("Z")) {
                int charge = Integer.parseInt(token[1]);
                float precursorMH = Float.parseFloat(token[2]);
                float precursorMz = ((precursorMH - (float) Composition.ChargeCarrierMass()) + charge * (float) Composition.ChargeCarrierMass()) / charge;

                Spectrum spec = new Spectrum();
                spec.setStartScanNum(startScanNum);
                spec.setEndScanNum(endScanNum);
                spec.setPrecursor(new Peak(precursorMz, 0, charge));
                spec.setTitle("ms2Index:" + ms2Index + " Z:" + token[1] + " MH:" + token[2]);
                specList.add(spec);
            } else if (token.length == 2)    // a peak
            {
                assert (specList != null);
                float mass = Float.parseFloat(token[0]);
                float intensity = Float.parseFloat(token[1]);
                for (Spectrum spec : specList)
                    spec.add(new Peak(mass, intensity, 1));
            }
        }

        if (specList != null) {
            for (Spectrum spec : specList)
                spec.outputMgf(out);
        }

        in.close();
        out.flush();
        out.close();

        System.out.println("Done.");
    }
}
