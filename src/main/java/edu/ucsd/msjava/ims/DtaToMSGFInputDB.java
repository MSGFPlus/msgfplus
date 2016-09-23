package edu.ucsd.msjava.ims;

import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class DtaToMSGFInputDB {
    public static void main(String argv[]) throws Exception {
        if (argv.length != 1)
            printUsageAndExit("Illegal parameter");

        File dtaPath = new File(argv[0]);
        if (!dtaPath.exists())
            printUsageAndExit("File does not exist.");

        List<File> dtaFileList = new ArrayList<File>();

        if (dtaPath.isDirectory()) {
            for (File f : dtaPath.listFiles()) {
                if (f.getName().endsWith("_dta.txt"))
                    dtaFileList.add(f);
            }
        } else {
            if (dtaPath.getName().endsWith("_dta.txt"))
                dtaFileList.add(dtaPath);
        }

        if (dtaFileList.size() == 0)
            printUsageAndExit("No _dta.txt file!");

        makeMSGFInput(dtaFileList);
    }

    public static void printUsageAndExit(String message) {
        if (message != null)
            System.out.println(message);
        System.out.println("Usage: java DtaToMSGFInputDB DTAPath");
        System.exit(-1);
    }

    public static void makeMSGFInput(List<File> dtaFileList) throws Exception {
        for (File dtaFile : dtaFileList) {
            String dtaFilePath = dtaFile.getAbsolutePath();
            String msgfInputFilePath = dtaFilePath.substring(0, dtaFilePath.lastIndexOf("_dta.txt")) + "_msgfInput.txt";
            File msgfInputFile = new File(msgfInputFilePath);
            System.out.println(dtaFile.getName() + "->" + msgfInputFile.getName());
            makeMSGFInput(dtaFile, msgfInputFile);
        }
    }

    public static void makeMSGFInput(File dtaFile, File msgfInputFile) throws Exception {
        PrintStream msgfOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(msgfInputFile)));
        String header = "#SpectrumFile\tScan#\tAnnotation\tPrecursorMz\tCharge\tPrecursorIntensity\tFrameNum\tFromScan\tToScan";
        msgfOut.println(header);

//		String dtaFilePath = dtaFile.getAbsolutePath();
//		File mgfFile = new File(dtaFilePath.substring(0, dtaFilePath.lastIndexOf("_dta.txt")) + ".mgf");
//		String mgfFileName = mgfFile.getName();

        BufferedLineReader in = new BufferedLineReader(dtaFile.getPath());
        String s;

        int lineNum = 0;
        int origSpecIndex = 0;
        while ((s = in.readLine()) != null) {
            lineNum++;
            if (s.startsWith("===")) {
                ++origSpecIndex;
                String metaInfo = s.substring(s.indexOf('"') + 1, s.lastIndexOf(".dta"));
                String[] token = metaInfo.split("\\.");
                if (token.length != 9) {
                    System.out.println("Syntax Error in Line " + lineNum + ": " + s);
                    System.exit(-1);
                }
                String annotation = token[0] + "." + token[1].replaceAll("!", "").replaceAll("@", "+15.995").replaceAll("\\*", "+15.995") + "." + token[2];
                int charge = Integer.parseInt(token[3]);
                int frameNum = Integer.parseInt(token[4]);
                int fromScan = Integer.parseInt(token[5]);
                int toScan = Integer.parseInt(token[6]);
                float precursorIntensity = Float.parseFloat(token[7] + token[8]);

                String precursorStr = in.readLine();
                if (precursorStr == null)
                    break;
                String[] precursorToken = precursorStr.split("\\s+");
                if (precursorToken.length < 2) {
                    System.out.println("Syntax Error in line " + lineNum + ": " + precursorStr);
                    System.exit(-1);
                }
                float precursorMH = Float.parseFloat(precursorToken[0]);
                float precursorMz = (precursorMH - (float) Composition.ChargeCarrierMass()) / charge + (float) Composition.ChargeCarrierMass();

                msgfOut.println(dtaFile.getName() + "\t" + origSpecIndex + "\t" + annotation + "\t" + precursorMz + "\t" + charge + "\t" + precursorIntensity + "\t" + frameNum + "\t" + fromScan + "\t" + toScan);
            }
        }

        in.close();
        msgfOut.close();
    }
}
