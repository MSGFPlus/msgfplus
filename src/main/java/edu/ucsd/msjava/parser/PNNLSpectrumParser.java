package edu.ucsd.msjava.parser;

import edu.ucsd.msjava.msutil.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

public class PNNLSpectrumParser implements SpectrumParser {

    public static final String SCAN_TYPE_FILE_EXTENSION = "_ScanType.txt";

    public Spectrum readSpectrum(LineReader lineReader) {
        Spectrum spec = null;

        String buf;
        float prevMass = 0;
        boolean isSorted = true;

        while ((buf = lineReader.readLine()) != null) {
            if (buf.length() == 0) {
                if (spec != null) {
                    if (!isSorted)
                        Collections.sort(spec);
                    return spec;
                } else
                    continue;
            } else if (buf.startsWith("==")) {
                if (spec != null) {
                    System.out.println("There must be at least one empty line between spectra: " + buf);
                    System.exit(-1);
                }
                int lastDotIndex = buf.lastIndexOf('.');
                int secondLastDotIndex = buf.lastIndexOf('.', lastDotIndex - 1);
                int thirdLastDotIndex = buf.lastIndexOf('.', secondLastDotIndex - 1);
                int fourthLastDotIndex = buf.lastIndexOf('.', thirdLastDotIndex - 1);

                int scanNum = Integer.parseInt(buf.substring(fourthLastDotIndex + 1, thirdLastDotIndex));

                String annotation = buf;
                // first line of a spectrum
                buf = lineReader.readLine();
                if (buf == null || buf.trim().length() == 0) {
                    System.out.println("Error while parsing _Dta.txt file: " + annotation);
                    System.out.println("No spectrum!");
                    System.exit(-1);
                }

                spec = new Spectrum();
                String[] token = buf.split("\\s+");
                float mPlusH = Float.parseFloat(token[0]);
                int charge = Integer.parseInt(token[1].substring(token[1].indexOf('=') + 1));
                float precursorMz = (mPlusH - (float) Composition.ChargeCarrierMass()) / charge + (float) Composition.ChargeCarrierMass();
                spec.setPrecursor(new Peak(precursorMz, 0, charge));
                spec.setScanNum(scanNum);
            } else if (Character.isDigit(buf.charAt(0)))    // peak
            {
                if (spec == null) {
                    System.out.println("Error while parsing _Dta.txt file.");
                    System.out.println("Header line is missing: " + buf);
                    System.exit(-1);
                }
                String[] token2 = buf.split("\\s+");
                if (token2.length != 2)
                    continue;
                float mass = Float.parseFloat(token2[0]);
                if (isSorted && mass < prevMass)
                    isSorted = false;

                float intensity = Float.parseFloat(token2[1]);
                spec.add(new Peak(mass, intensity, 1));
                prevMass = mass;
            }
        }
        return spec;
    }

    @Override
    public Map<Integer, SpectrumMetaInfo> getSpecMetaInfoMap(BufferedRandomAccessLineReader lineReader) {
        Hashtable<Integer, SpectrumMetaInfo> specIndexMap = new Hashtable<Integer, SpectrumMetaInfo>();
        String buf;
        long offset = 0;
        int specIndex = 0;
        while ((buf = lineReader.readLine()) != null) {
            if (buf.startsWith("==")) {
//				specIndexMap.put(++specIndex, offset);
                ++specIndex;
                int lastDotIndex = buf.lastIndexOf('.');
                int secondLastDotIndex = buf.lastIndexOf('.', lastDotIndex - 1);
                int thirdLastDotIndex = buf.lastIndexOf('.', secondLastDotIndex - 1);
                int fourthLastDotIndex = buf.lastIndexOf('.', thirdLastDotIndex - 1);

                String annotation = buf;
                // first line of a spectrum
                buf = lineReader.readLine();
                if (buf == null || buf.trim().length() == 0) {
                    System.out.println("Error while parsing _Dta.txt file: " + annotation);
                    System.out.println("No spectrum!");
                    System.exit(-1);
                }

                String[] token = buf.split("\\s+");
                float mPlusH = Float.parseFloat(token[0]);
                int charge = Integer.parseInt(token[1].substring(token[1].indexOf('=') + 1));
                float precursorMz = (mPlusH - (float) Composition.ChargeCarrierMass()) / charge + (float) Composition.ChargeCarrierMass();

                SpectrumMetaInfo metaInfo = new SpectrumMetaInfo();
                metaInfo.setID("index=" + (specIndex - 1));
                metaInfo.setPrecursorMz(precursorMz);
                metaInfo.setPosition(offset);
                specIndexMap.put(specIndex, metaInfo);
            }
            offset = lineReader.getPosition();
        }
        return specIndexMap;
    }

//	static class ScanType
//	{
//		public ScanType(ActivationMethod activationMethod,
//				boolean isHighPrecision) {
//			this.activationMethod = activationMethod;
//			this.isHighPrecision = isHighPrecision;
//		}
//		
//		ActivationMethod getActivationMethod() {
//			return activationMethod;
//		}
//		boolean isHighPrecision() {
//			return isHighPrecision;
//		}
//
//		private ActivationMethod activationMethod;
//		private boolean isHighPrecision;
//	}

    static HashMap<Integer, ScanType> getScanTypeMap(String fileName) {
        File specFile = new File(fileName);
        String scanTypeFileName =
                specFile.getAbsoluteFile().getParentFile().getPath()
                        + File.separator
                        + specFile.getName().substring(0, specFile.getName().lastIndexOf('_'))
                        + PNNLSpectrumParser.SCAN_TYPE_FILE_EXTENSION;
        File scanTypeFile = new File(scanTypeFileName);

        if (!scanTypeFile.exists())
            return null;

        HashMap<Integer, ScanType> scanNumScanTypeMap = new HashMap<Integer, ScanType>();

        BufferedLineReader in = null;
        try {
            in = new BufferedLineReader(scanTypeFile.getPath());
        } catch (IOException e) {
            e.printStackTrace();
        }

        String s;

        s = in.readLine();    // header
        boolean hasScanTimes = false;
        String[] hTokens = s.split("\t");
        if (hTokens.length > 3 && hTokens[3].toLowerCase().contains("time")) {
            hasScanTimes = true;
        }

        while ((s = in.readLine()) != null) {
            String[] token = s.split("\t");
            if (token.length < 3)
                continue;

            int scanNum = Integer.parseInt(token[0]);
            String scanType = token[1].toLowerCase();

            ActivationMethod method = null;
            if (scanType.contains("etcid"))
                method = ActivationMethod.ETD;
            else if (scanType.contains("ethcd"))
                method = ActivationMethod.ETD;
            else if (scanType.contains("cid"))
                method = ActivationMethod.CID;
            else if (scanType.contains("etd"))
                method = ActivationMethod.ETD;
            else if (scanType.contains("hcd"))
                method = ActivationMethod.HCD;
            else if (scanType.contains("pqd"))
                method = ActivationMethod.PQD;

            boolean isHighPrecision = false;
            if (scanType.contains("hms"))
                isHighPrecision = true;

            int msLevel = Integer.parseInt(token[2]);

            float scanTime = -1;
            if (hasScanTimes && token.length > 3) {
                scanTime = Float.parseFloat(token[3]);
            }

            if (method != null) {
                scanNumScanTypeMap.put(scanNum, new ScanType(method, isHighPrecision, msLevel, scanTime));
            }
        }

        if (in != null) {
            try {
                in.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return scanNumScanTypeMap;
    }

    public static void main(String argv[]) throws Exception {
        long time = System.currentTimeMillis();
        String fileName = System.getProperty("user.home") + "/Research/ToolDistribution/PNNLTest/QC_Shew_08_04_pt5_b_22Jan09_Owl_09-01-04_dta.txt";
        SpectraIterator itr = new SpectraIterator(fileName, new PNNLSpectrumParser());
        int numSpecs = 0;
        HashSet<Integer> scanNumSet = new HashSet<Integer>();
        while (itr.hasNext()) {
            Spectrum spec = itr.next();
            numSpecs++;
            if (scanNumSet.contains(spec.getScanNum())) {
                System.out.println(spec.getScanNum());
            } else
                scanNumSet.add(spec.getScanNum());
//			System.out.println(spec+ "\t" + spec.getScanNum()+"\t"+(spec.getParentMass()+(float)Composition.ChargeCarrierMass)+"\t"+spec.getCharge());
        }
        System.out.println("NumSpecs: " + numSpecs);
        System.out.println("Time: " + (System.currentTimeMillis() - time));

        time = System.currentTimeMillis();
        SpectraMap map = new SpectraMap(fileName, new PNNLSpectrumParser());
        numSpecs = 0;
        for (int specIndex : map.getSpecIndexList()) {
            Spectrum spec = map.getSpectrumBySpecIndex(specIndex);
            numSpecs++;
//			System.out.println(spec+ "\t" + spec.getScanNum()+"\t"+(spec.getParentMass()+(float)Composition.ChargeCarrierMass)+"\t"+spec.getCharge());
        }
        System.out.println("NumSpecs: " + numSpecs);
        System.out.println("Time: " + (System.currentTimeMillis() - time));
    }
}
