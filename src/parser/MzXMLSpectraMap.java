package parser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import msutil.ActivationMethod;
import msutil.Peak;
import msutil.Spectrum;
import msutil.SpectrumAccessorByScanNum;

import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.Scan;
import org.systemsbiology.jrap.stax.ScanHeader;

/**
 * A data structure that allows random access of the mzXML file.
 * @author jung
 *
 */
public class MzXMLSpectraMap implements SpectrumAccessorByScanNum {
	private	MSXMLParser parser;

	// the pattern to extract the retention time string from text field
	private static Pattern rtPattern = Pattern.compile("\\D*(\\d+\\.\\d*|\\d*\\.\\d+|\\d+)(m|M|s|S)?");

	// if(maxMSLevel >= minMSLevel > 0) only spectra within [minMSLevel, maxMSLevel] will be returned
	private int minMSLevel = 0;		// inclusive
	private int maxMSLevel = 0;		// inclusive

	/***** CONSTRUCTORS *****/
	/**
	 * Constructor taking the file path to the mzXML file.
	 * @param fileName the path to the file.
	 */
	public MzXMLSpectraMap(String fileName) {
		parser = new MSXMLParser(fileName);
	}

	/**
	 * Setter to set msLevel.
	 * @param minMSLevel minimum msLevel to be considered (inclusive).
	 * @param maxMSLevel maximum msLevel to be considered (inclusive).
	 * @return this object.
	 */
	public MzXMLSpectraMap msLevel(int minMSLevel, int maxMSLevel) { this.minMSLevel = minMSLevel; this.maxMSLevel = maxMSLevel; return this; }

	/**
	 * Get the spectrum by scan number. The scanNumber is the absolute position 
	 * of this scan in the file starting at 1. This is different than the num in
	 * the mzXML file.
	 * @param scanNumber the scan number of the spectrum to look for. 
	 * @return the Spectrum object or null if not found.
	 */
	public Spectrum getSpectrumByScanNum(int scanNumber) {
		Scan scanObj = parser.rap(scanNumber);

		if (scanObj==null)  return null;
		int msLevel = scanObj.getHeader().getMsLevel();
		if (minMSLevel > 0 && minMSLevel <= maxMSLevel && (msLevel < minMSLevel || msLevel > maxMSLevel))
			return null;
		// get peak list array (mass, intensities) pairs
		double[][] peakList = scanObj.getMassIntensityList();

		ScanHeader header = scanObj.getHeader();
		int precursorCharge = header.getPrecursorCharge();
		precursorCharge = (precursorCharge < 0) ? 0 : precursorCharge; 

		Spectrum spec = new Spectrum(header.getPrecursorMz(), precursorCharge, header.getPrecursorIntensity());
		spec.setScanNum(header.getNum());

		// parse retention time. Note that retention time is a required field
		String rtStr = header.getRetentionTime();
		if (rtStr!=null) {
			Matcher matcher = rtPattern.matcher(rtStr);
			if (matcher.find() && matcher.groupCount() > 0) {
				float rtFloat = Float.parseFloat(rtStr.substring(matcher.start(1), matcher.end(1)));
				if (matcher.groupCount() > 1) {
					String timeScale = rtStr.substring(matcher.start(2), matcher.end(2));
					if (timeScale.equals("M") || timeScale.equals("m"))    rtFloat *= 60.0;
				}
				spec.setRt(rtFloat);
			}
		}

		// set ms level
		spec.setMsLevel(header.getMsLevel());

		// add activation method
		String activationName = header.getActivationMethod();
		if(activationName != null)
		{
			ActivationMethod method = ActivationMethod.get(activationName);
			if(method == null)
			{
				ActivationMethod.register(activationName);
				method = ActivationMethod.get(activationName);
			}
			spec.setActivationMethod(method);
		}

		// add peaks
		boolean sorted = true;
		float prevMass = 0;
		for (int j=0; j<peakList[0].length; j++) {
			float mass = (float)peakList[0][j];
			float intensity =  (float)peakList[1][j];
			spec.add(new Peak(mass, intensity, 1));
			if(sorted && mass < prevMass)
				sorted = false;
			else
				prevMass = mass;
		}
		if(!sorted)    Collections.sort(spec);

		return spec;
	}



	/**
	 * Get the number of scans in this file.
	 * @return the number of total scans.
	 */
	public int getScanCount() {
		return parser.getScanCount();
	}

	private ArrayList<Integer> scanNumList = null;
	
	@Override
	public ArrayList<Integer> getScanNumList() {
		if(scanNumList == null)
		{
			scanNumList = new ArrayList<Integer>();
			for(int scanNumber = 0; scanNumber<=parser.getScanCount(); scanNumber++)
			{
				Scan scanObj = parser.rap(scanNumber);
				if (scanObj != null && scanObj.getHeader() != null && scanObj.getHeader().getMsLevel() == 2)
					scanNumList.add(scanNumber);
			}
		}
		return scanNumList;
	}
}
