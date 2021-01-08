package edu.ucsd.msjava.parser;

import edu.ucsd.msjava.msutil.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Map;

import static edu.ucsd.msjava.misc.TextParsingUtils.isInteger;

/**
 * This class enables to parse spectrum file with mgf format.
 *
 * @author sangtaekim
 */
public class MgfSpectrumParser implements SpectrumParser {

    private long linesRead;

    private long negativePolarityWarningCount;

    private long scanMissingWarningCount;

    /**
     * Number of scans where we could not determine the scan number
     * This method is required by interface SpectrumParser
     * @return
     */
    public long getScanMissingWarningCount()
    {
        return scanMissingWarningCount;
    }

    /**
     * Amino acid set to be used to parse "SEQ="
     */
    private AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();

    /**
     * Specify amino acid set to be used to parse "SEQ=" field.
     *
     * @param aaSet amino acid set.
     * @return this object.
     */
    public MgfSpectrumParser aaSet(AminoAcidSet aaSet) {
        this.aaSet = aaSet;
        linesRead = 0;
        negativePolarityWarningCount = 0;
        scanMissingWarningCount = 0;
        return this;
    }

    /**
     * Implementation of readSpectrum method. Implicitly lineReader points to the start of a spectrum.
     * Reads mgf file line by line until the spectrum ends, generate a Spectrum object and returns it.
     * If it cannot read a spectrum, it returns null.
     *
     * @param lineReader a LineReader object points to the start of a spectrum
     * @return a spectrum object. null if no spectrum can be generated.
     */
    public Spectrum readSpectrum(LineReader lineReader) {
        Spectrum spec = null;
        String title = null;

        float precursorMz = 0;
        float precursorIntensity = 0;
        int precursorCharge = 0;
        ActivationMethod activation = null;
        float elutionTimeSeconds = 0;
//		Float toleranceVal = null;
//		Tolerance.Unit toleranceUnit = null;

        String buf;
        boolean parse = false;   // parse only after the BEGIN IONS
        boolean sorted = true;
        float prevMass = 0;

        while (true) {
            String dataLine = (buf = lineReader.readLine());
            if (dataLine == null)
                break;

            if (linesRead == 0) {
                buf = BufferedRandomAccessLineReader.stripBOM(buf);
            }
            linesRead++;

            if (buf.length() == 0)
                continue;

            if (buf.startsWith("BEGIN IONS")) {
                parse = true;
                spec = new Spectrum();
            } else if (parse) {
                if (Character.isDigit(buf.charAt(0))) {
                    assert (spec != null);
                    String[] token = buf.split("\\s+");
                    if (token.length < 2)
                        continue;
                    float mass = Float.parseFloat(token[0]);
                    if (sorted && mass < prevMass)
                        sorted = false;
                    else
                        prevMass = mass;
                    float intensity = Float.parseFloat(token[1]);
                    spec.add(new Peak(mass, intensity, 1));
                } else if (buf.startsWith("TITLE")) {
                    title = buf.substring(buf.indexOf('=') + 1);
                    spec.setTitle(title);
//  				spec.setID(title);
                } else if (buf.startsWith("CHARGE")) {
                    // Charge state, e.g. CHARGE=2+
                    // Extract the text after the equals sign
                    String chargeStr = buf.substring(buf.indexOf("=") + 1).trim();

                    // Only use the charge state if there is a single value listed
                    // We will leave precursorCharge as 0 if the mgf file has lines like this:
                    //  CHARGE=2+ and 3+
                    //  CHARGE=2+,3+
                    // First split on whitespace
                    String[] chargeStrToken = chargeStr.split("\\s+");
                    if (chargeStrToken.length == 1) {
                        // Only one charge state is listed
                        // Now split on commas
                        String[] multipleChargeToken = chargeStr.split(",");
                        if (chargeStr.length() > 0 && multipleChargeToken.length == 1) {
                            // Only one value is present
                            if (chargeStr.startsWith("+")) {
                                // The charge is listed as +2 (which is non-standard)
                                // Remove the plus sign
                                chargeStr = chargeStr.substring(1);
                            } else if (chargeStr.charAt(chargeStr.length() - 1) == '+') {
                                // The charge is listed as 2+ (which is standard)
                                // Remove the plus sign
                                chargeStr = chargeStr.substring(0, chargeStr.length() - 1);
                            } else if (chargeStr.startsWith("-")) {
                                // The charge is listed as -2
                                // This is a negative charge, which means negative scan polarity
                                // MS-GF+ does not yet support this, but we'll store the charge anyway (as a positive number)
                                warnNegativePolarity(buf);
                                chargeStr = chargeStr.substring(1);
                                spec.setScanPolarity(Spectrum.Polarity.NEGATIVE);
                            } else if (chargeStr.charAt(chargeStr.length() - 1) == '-') {
                                // The charge is listed as 2-
                                // This is a negative charge, which means negative scan polarity
                                // MS-GF+ does not yet support this, but we'll store the charge anyway (as a positive number)
                                warnNegativePolarity(buf);
                                chargeStr = chargeStr.substring(0, chargeStr.length() - 1);
                                spec.setScanPolarity(Spectrum.Polarity.NEGATIVE);
                            }

                            // We should now have an integer to parse
                            precursorCharge = Integer.valueOf(chargeStr);
                        }
                    }
                } else if (buf.startsWith("SEQ")) {
                    String annotationStr = buf.substring(buf.lastIndexOf('=') + 1);
                    if (spec.getAnnotation() == null)
                        spec.setAnnotation(new Peptide(annotationStr, aaSet));
                    spec.addSEQ(annotationStr);
                } else if (buf.startsWith("PEPMASS")) {
                    String[] token = buf.substring(buf.indexOf("=") + 1).split("\\s+");
                    precursorMz = Float.valueOf(token[0]);
                } else if (buf.startsWith("SCANS")) {
                    if (buf.matches(".+=\\d+-\\d+"))    // e.g. SCANS=953-959
                    {
                        int startScanNum = Integer.parseInt(buf.substring(buf.indexOf('=') + 1, buf.lastIndexOf('-')));
                        int endScanNum = Integer.parseInt(buf.substring(buf.lastIndexOf('-') + 1));
                        spec.setStartScanNum(startScanNum);
                        spec.setEndScanNum(endScanNum);
                    } else {
                        // Look for a single integer after the equals sign
                        try {
                            int scanNum = Integer.valueOf(buf.substring(buf.indexOf("=") + 1));
                            spec.setScanNum(scanNum);
                        } catch (NumberFormatException e) {
                            // Not an integer; the scan number will be the zero based sequence number of the spectrum
                        }
                    }
                } else if (buf.startsWith("ACTIVATION")) {
                    String activationName = buf.substring(buf.indexOf("=") + 1);
                    activation = ActivationMethod.get(activationName);
                    spec.setActivationMethod(activation);
                } else if (buf.startsWith("RTINSECONDS")) {
                    String[] token = buf.substring(buf.indexOf("=") + 1).split("\\s+");
                    elutionTimeSeconds = Float.valueOf(token[0]);
                }
//  			else if(buf.startsWith("TOL="))
//  			{
//  				String tolStr = buf.substring(buf.indexOf("=")+1);
//  				float toleranceValue = Float.parseFloat(tolStr);
//  				if(toleranceValue > 0)
//  				{
//  					toleranceVal = toleranceValue;
//  				}
//  			}
//  			else if(buf.startsWith("TOLU="))
//  			{
//  				String tolUnitStr = buf.substring(buf.indexOf("=")+1);
//  				if(tolUnitStr.equalsIgnoreCase("ppm"))
//  					toleranceUnit = Tolerance.Unit.PPM;
//  				else if(tolUnitStr.equalsIgnoreCase("Da"))
//  					toleranceUnit = Tolerance.Unit.Da;
//  			}
                else if (buf.startsWith("END IONS")) {
                    assert (spec != null);
                    if (spec.getScanNum() < 0 && title != null) {
                        if (title.matches("Scan:\\d+\\s.+")) {
                            // Title line is of the form Scan:ScanNumber AdditionalText
                            // for example, "Scan:8492 Charge:2"
                            // Extract the integer after "Scan:"
                            // Split on spaces
                            String[] token = title.split("\\s++");
                            int scanNum = Integer.parseInt(token[0].substring("Scan:".length()));
                            spec.setScanNum(scanNum);

                        } else if (title.matches(".+\\.\\d+\\.\\d+\\.\\d+$") ||
                                title.matches(".+\\.\\d+\\.\\d+\\.$")) {
                            // Title line is of the form DatasetName.ScanStart.ScanEnd.Charge or DatasetName.ScanStart.ScanEnd.
                            // for example, DatasetName.8492.8492.2
                            extractScanRangeFromTitle(spec, title);

                        } else if (title.contains(".") && title.contains(" ")) {
                            // Remove text after the first space and try to match DatasetName.ScanStart.ScanEnd.Charge
                            // Split on periods
                            String titleStart = title.substring(0, title.indexOf(' '));
                            extractScanRangeFromTitle(spec, titleStart);
                        } else {
                            warnScanNotFoundInTitle(title);
                        }

                        //Match result = dtaStyleMatcher.matcher(spec.Title)
                    }
                    spec.setPrecursor(new Peak(precursorMz, precursorIntensity, precursorCharge));
                    if (elutionTimeSeconds > 0) {
                        spec.setRt(elutionTimeSeconds);
                        spec.setRtIsSeconds(true);
                    }
//  				if(toleranceVal != null && toleranceUnit != null)
//  				{
//  					Tolerance precursorTolerance = new Tolerance(toleranceVal, toleranceUnit);
//  					spec.setPrecursorTolerance(precursorTolerance);
//  				}
                    if (!sorted)
                        Collections.sort(spec);

                    return spec;
                }
            }
        }
        return null;
    }

    /**
     * Extract start and end scan from the title if it is of the form:
     * DatasetName.ScanStart.ScanEnd.Charge
     *
     * @param spec  Spectrum
     * @param title Title line
     */
    private void extractScanRangeFromTitle(Spectrum spec, String title) {
        // Split on periods
        String[] token = title.split("\\.");
        String candidateStartScan;
        String candidateEndScan;

        if (token.length > 3) {
            // Assume DatasetName.ScanStart.ScanEnd.Charge
            // For example: DatasetName.10418.10418.4
            candidateStartScan = token[token.length - 3];
            candidateEndScan = token[token.length - 2];
        } else if (token.length == 3 && title.endsWith(".")) {
            // Charge not specified, but title does end with a period
            // In this case, .split() only returns 3 items

            // Assume DatasetName.ScanStart.ScanEnd.
            // For example: DatasetName.40193.40193.
            candidateStartScan = token[token.length - 2];
            candidateEndScan = token[token.length - 1];
        } else {
            warnScanNotFoundInTitle(title);
            return;
        }

        boolean success = false;
        if (isInteger(candidateStartScan)) {
            int startScanNum = Integer.parseInt(candidateStartScan);
            spec.setStartScanNum(startScanNum);
            success = true;
        }

        if (isInteger(candidateEndScan)) {
            int endScanNum = Integer.parseInt(candidateEndScan);
            spec.setEndScanNum(endScanNum);
        }

        if (!success) {
            warnScanNotFoundInTitle(title);
        }
    }

    /**
     * Implementation of getSpecIndexMap object. Reads the entire spectrum file and
     * generates a map from a spectrum index to the file position of the spectrum.
     *
     * @param lineReader a LineReader object that points to the start of a file.
     * @return A map from spectrum indexes to the spectrum meta information.
     */
    public Map<Integer, SpectrumMetaInfo> getSpecMetaInfoMap(BufferedRandomAccessLineReader lineReader) {
        Hashtable<Integer, SpectrumMetaInfo> specIndexMap = new Hashtable<Integer, SpectrumMetaInfo>();
        String buf;
        long offset = 0;
        int specIndex = 0;
        SpectrumMetaInfo metaInfo = null;
        while (true) {
            String dataLine = (buf = lineReader.readLine());
            if (dataLine == null)
                break;

            if (offset == 0 && lineReader.getBOMLength() > 0) {
                offset += lineReader.getBOMLength();
            }

            if (buf.startsWith("BEGIN IONS")) {
                specIndex++;
                metaInfo = new SpectrumMetaInfo();
                metaInfo.setPosition(offset);
                metaInfo.setID("index=" + String.valueOf(specIndex - 1));
                specIndexMap.put(specIndex, metaInfo);
            } else if (buf.startsWith("TITLE")) {
                String title = buf.substring(buf.indexOf('=') + 1);
                metaInfo.setAdditionalInfo("title", title);
            } else if (buf.startsWith("PEPMASS")) {
                String[] token = buf.substring(buf.indexOf("=") + 1).split("\\s+");
                float precursorMz = Float.valueOf(token[0]);
                metaInfo.setPrecursorMz(precursorMz);
            }

            offset = lineReader.getPosition();
        }
        return specIndexMap;
    }

    private void warnNegativePolarity(String currentLine) {
        negativePolarityWarningCount++;
        if (negativePolarityWarningCount > MAX_NEGATIVE_POLARITY_WARNINGS)
            return;

        if (negativePolarityWarningCount == 1) {
            System.out.println(
                    "Warning: negative precursor charge found, indicating a negative polarity spectrum; " +
                    "you likely need to use a negative charge carrier");
        }
        System.out.println("Negative charge found on line " + Long.toString(linesRead) + ": " + currentLine);

        if (negativePolarityWarningCount == MAX_NEGATIVE_POLARITY_WARNINGS) {
            System.out.println("Additional warnings regarding negative polarity will not be shown");
        }
    }

    void warnScanNotFoundInTitle(String title) {
        scanMissingWarningCount++;
        if (scanMissingWarningCount <= MAX_SCAN_MISSING_WARNINGS) {
            System.out.println("Unable to extract the scan number from the title: " + title);
            if (scanMissingWarningCount == 1) {
                System.out.println("Expected format is DatasetName.ScanStart.ScanEnd.Charge");
            }
        }
    }

    // test code
    public static void main(String argv[]) throws Exception {
        long time = System.currentTimeMillis();
        String mgfFile = "/Users/sangtaekim/Research/Data/PNNL/IPYS_TD_Scere010_Orbitrap_001a.mgf";
//		String mgfFile = "/Users/sangtaekim/Research/Data/AgilentQTOF/notAnnotatedAgilentQTOF.mgf";

	    /*
        // SpectraIterator test
		MgfSpectrumParser parser = new MgfSpectrumParser();
	    SpectraIterator itr = new SpectraIterator(mgfFile, parser);
	    int size = 0;
	    while(itr.hasNext())
	    {
	    	Spectrum spec = itr.next();
	    	size++;
	    	System.out.println(spec.getScanNum()+" "+spec.getPrecursorPeak());
	    }
	    System.out.println("Size: " + size);
	    */
        //  SpectraMap test

	    /*	SpectraMap test
	    SpectraMap map = new SpectraMap(mgfFile, new MgfSpectrumParser());
	    Spectrum spec = map.getSpectrumByScanNum(1585);
	    System.out.println(spec.getScanNum() + " " + spec.getPrecursorPeak());
	    */

//	    SpectraContainer container = new SpectraContainer(mgfFile, new MgfSpectrumParser());
//	    for(Spectrum spec : container)
//	    	System.out.println(spec.getScanNum() + " " + spec.getPrecursorPeak());
        ArrayList<Spectrum> specContainer = new ArrayList<Spectrum>();
        SpectraIterator iterator = new SpectraIterator(mgfFile, new MgfSpectrumParser());
        while (iterator.hasNext())
            specContainer.add(iterator.next());
        System.out.println("Time: " + (System.currentTimeMillis() - time));
    }
}
