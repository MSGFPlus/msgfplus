package edu.ucsd.msjava.mzml;

import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.SpectrumParser;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

import java.io.File;
import java.util.Iterator;


public class MzMLSpectraIterator implements Iterator<edu.ucsd.msjava.msutil.Spectrum>, Iterable<edu.ucsd.msjava.msutil.Spectrum> {
    private final MzMLUnmarshaller unmarshaller;
    private final int minMSLevel;        // inclusive
    private final int maxMSLevel;        // exclusive

    private MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> itr;
    private boolean hasNext;
    private edu.ucsd.msjava.msutil.Spectrum currentSpectrum = null;
    private long negativeChargeWarningCount;

    public MzMLSpectraIterator(MzMLAdapter mzmlAdapter) {
        unmarshaller = mzmlAdapter.getUnmarshaller();
        minMSLevel = mzmlAdapter.getMinMSLevel();
        maxMSLevel = mzmlAdapter.getMaxMSLevel();
        negativeChargeWarningCount = 0;

        itr = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", uk.ac.ebi.jmzml.model.mzml.Spectrum.class);
        currentSpectrum = parseNextSpectrum();
        hasNext = currentSpectrum != null;
    }

    public boolean hasNext() {
        return hasNext;
    }

    /**
     * Get next spectrum.
     *
     * @return the next spectrum.
     */
    public edu.ucsd.msjava.msutil.Spectrum next() {
        Spectrum curSpecCopy = currentSpectrum;
        currentSpectrum = parseNextSpectrum();
        if (currentSpectrum == null)
            hasNext = false;

        if (curSpecCopy.getScanPolarity() == Spectrum.Polarity.NEGATIVE) {
            warnNegativeCharge(curSpecCopy);
        }
        return curSpecCopy;
    }

    public edu.ucsd.msjava.msutil.Spectrum parseNextSpectrum() {
        edu.ucsd.msjava.msutil.Spectrum spec = null;
        uk.ac.ebi.jmzml.model.mzml.Spectrum jmzSpec = null;

        while (itr.hasNext()) {
            jmzSpec = itr.next();
            spec = SpectrumConverter.getSpectrumFromJMzMLSpec(jmzSpec);
            if (spec.getMSLevel() < minMSLevel || spec.getMSLevel() > maxMSLevel)
                continue;
            else
                return spec;
        }

        return null;
    }

    public void remove() {
        throw new UnsupportedOperationException("SpectraIterator.remove() not implemented");
    }

    public Iterator<edu.ucsd.msjava.msutil.Spectrum> iterator() {
        return this;
    }


    public static void main(String argv[]) throws Exception {
//		List<Logger> loggers = Collections.<Logger>list(LogManager.getCurrentLoggers());
//		loggers.add(LogManager.getRootLogger());
//		for ( Logger logger : loggers ) {
//		    logger.setLevel(Level.OFF);
//		}		
        test();
    }

    public static void test() throws Exception {
        File xmlFile = new File("/cygwin/home/kims336/Research/Data/JMzReader/example.mzML");
        xmlFile = new File("/cygwin/home/kims336/Research/Data/JMzReader/small.pwiz.1.1.mzML");

        MzMLAdapter adapter = new MzMLAdapter(xmlFile);
        MzMLSpectraIterator itr = new MzMLSpectraIterator(adapter);
        while (itr.hasNext()) {
            edu.ucsd.msjava.msutil.Spectrum spec = itr.next();
            System.out.println("-----------");
            System.out.println(spec.getID());
            System.out.println(spec.getSpecIndex());
            System.out.println(spec.getMSLevel());
            if (spec.getMSLevel() == 2) {
                System.out.println(spec.getPrecursorPeak().getMz());
                System.out.println(spec.getPrecursorPeak().getCharge());
                System.out.println(spec.getActivationMethod().getName());
            }
        }
    }

    private void warnNegativeCharge(Spectrum currentSpectrum) {
        negativeChargeWarningCount++;
        if (negativeChargeWarningCount > SpectrumParser.MAX_NEGATIVE_CHARGE_WARNINGS)
            return;

        if (negativeChargeWarningCount == 1) {
            System.out.println("Warning: MS-GF+ does not support negative mode precursor ions (i.e. negative scan polarity)");
        }
        System.out.println("Ignoring the charge state defined for scan " + Long.toString(currentSpectrum.getScanNum()));

        if (negativeChargeWarningCount == SpectrumParser.MAX_NEGATIVE_CHARGE_WARNINGS) {
            System.out.println("Additional warnings regarding negative mode precursor ions will not be shown");
        }
    }

}