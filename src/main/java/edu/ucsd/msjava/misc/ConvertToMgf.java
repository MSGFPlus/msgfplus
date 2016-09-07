package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Iterator;


public class ConvertToMgf {
    public static void main(String argv[]) throws Exception {
        boolean writeActivationMethod = false;
        ActivationMethod activationMethod = null;
        File source = null;
        File target = null;
        int specIndex = -1;
        int scanNum = -1;
        int charge = -1;
        String id = null;
        boolean eraseCharge = false;

        for (int i = 0; i < argv.length; i += 2) {
            if (!argv[i].startsWith("-") || i + 1 >= argv.length)
                printUsageAndExit("Illegal parameters");
            else if (argv[i].equalsIgnoreCase("-s")) {
                source = new File(argv[i + 1]);
                if (!source.exists())
                    printUsageAndExit(argv[i + 1] + " doesn't exist!");
            } else if (argv[i].equalsIgnoreCase("-t")) {
                target = new File(argv[i + 1]);
                if (!target.getName().endsWith(".mgf"))
                    printUsageAndExit(argv[i + 1] + " should end with .mgf!");
            } else if (argv[i].equalsIgnoreCase("-w")) {
                if (argv[i + 1].equals("0"))
                    writeActivationMethod = false;
                else if (argv[i + 1].equals("1"))
                    writeActivationMethod = true;
            } else if (argv[i].equalsIgnoreCase("-m"))    // Fragmentation method
            {
                // (0: written in the spectrum, 1: CID , 2: ETD, 3: HCD)
                if (argv[i + 1].equalsIgnoreCase("0")) {
                    activationMethod = null;
                } else if (argv[i + 1].equalsIgnoreCase("1")) {
                    activationMethod = ActivationMethod.CID;
                } else if (argv[i + 1].equalsIgnoreCase("2")) {
                    activationMethod = ActivationMethod.ETD;
                } else if (argv[i + 1].equalsIgnoreCase("3")) {
                    activationMethod = ActivationMethod.HCD;
                } else
                    printUsageAndExit("Illegal activation method: " + argv[i + 1]);
            } else if (argv[i].equalsIgnoreCase("-index")) {
                specIndex = Integer.parseInt(argv[i + 1]);

            } else if (argv[i].equalsIgnoreCase("-scan")) {
                scanNum = Integer.parseInt(argv[i + 1]);

            } else if (argv[i].equalsIgnoreCase("-id")) {
                id = argv[i + 1];
            } else if (argv[i].equalsIgnoreCase("-charge")) {
                charge = Integer.parseInt(argv[i + 1]);
            } else if (argv[i].equalsIgnoreCase("-eraseCharge")) {
                if (argv[i + 1].equals("0"))
                    eraseCharge = false;
                else if (argv[i + 1].equals("1"))
                    eraseCharge = true;
            } else {
                printUsageAndExit("Invalid parameter: " + argv[i]);
            }
        }

        if (source == null || target == null)
            printUsageAndExit("Invalid parameters!");
        convert(source, target, writeActivationMethod, activationMethod, id, scanNum, specIndex, charge, eraseCharge);
    }

    public static void printUsageAndExit(String message) {
        if (message != null)
            System.out.println(message);
        System.out.println("Usage: java ConvertToMgf\n" +
                "\t-s SourceFile or Directory\n" +
                "\t-t TargetFileName (*.mgf)\n" +
                "\t[-w 0/1] (0: don't write ActivationMethod (default), 1: write ActivationMethod)\n" +
                "\t[-m FragmentationMethodID] (0: convert all (Default), 1: CID, 2: ETD, 3: HCD)\n" +
                "\t[-index SpecIndex] (only write the spectrum of the specified index will be converted)\n" +
                "\t[-scan ScanNum] (only write the spectrum of the specified scan number will be converted)\n" +
                "\t[-id id] (only write the spectrum of the specified id will be converted)\n" +
                "\t[-charge charge] (only write the spectrum with the specified charge)\n" +
                "\t[-eraseCharge 0/1] (0: Don't erase charge, 1: Erase precursor charge)\n"
        );
        System.exit(-1);
    }

    public static void convert(File source, File target, boolean writeActivationMethod,
                               ActivationMethod activationMethod, String id, int scanNum, int specIndex, int charge,
                               boolean eraseCharge) throws Exception {
        PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(target)));

        File[] fileList;
        if (!source.isDirectory()) {
            fileList = new File[1];
            fileList[0] = source;
        } else
            fileList = source.listFiles();

        int numFileConverted = 0;
        for (File sourceFile : fileList) {
            Iterator<Spectrum> specItr = new SpectraAccessor(sourceFile).getSpecItr();
            if (specItr != null) {
                System.out.print(sourceFile.getName() + ": ");
                int numSpecs = convertFile(specItr, target, writeActivationMethod, activationMethod, id, scanNum, specIndex, charge, eraseCharge, out);
                System.out.println(numSpecs + " spectra converted.");
                numFileConverted++;
            }
        }
        out.close();
        System.out.println("Converted " + numFileConverted + " files.");
    }

    public static int convertFile(Iterator<Spectrum> specItr, File target, boolean writeActivationMethod,
                                  ActivationMethod activationMethod, String id, int scanNum, int specIndex, int charge, boolean eraseCharge,
                                  PrintStream out) throws Exception {
        int numSpecs = 0;
        while (specItr.hasNext()) {
            Spectrum spec = specItr.next();
            if (id != null && spec.getID().endsWith(id))
                continue;
            if (specIndex > 0 && spec.getSpecIndex() != specIndex)
                continue;
            if (scanNum > 0 && spec.getScanNum() != scanNum)
                continue;
            if (activationMethod != null && spec.getActivationMethod() != activationMethod)
                continue;
            if (charge >= 0 && spec.getCharge() != charge)
                continue;
            if (eraseCharge) {
                if (spec.getCharge() <= 4) {
                    float precursorMz;
                    if (spec.getIsolationWindowTargetMz() != null)
                        precursorMz = spec.getIsolationWindowTargetMz();
                    else
                        precursorMz = spec.getPrecursorPeak().getMz();
                    spec.setPrecursor(new Peak(precursorMz, 0, 0));
                }
            }
            spec.outputMgf(out, writeActivationMethod);
//			System.out.println(spec.getID() + " " + spec.getScanNum());
            numSpecs++;
        }
        return numSpecs;
    }
}
