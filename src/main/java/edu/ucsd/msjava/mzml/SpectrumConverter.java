package edu.ucsd.msjava.mzml;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Peak;
import uk.ac.ebi.jmzml.model.mzml.*;

import java.util.Collections;

public class SpectrumConverter {
    public static edu.ucsd.msjava.msutil.Spectrum getSpectrumFromJMzMLSpec(uk.ac.ebi.jmzml.model.mzml.Spectrum jmzMLSpec) {
        edu.ucsd.msjava.msutil.Spectrum spec = new edu.ucsd.msjava.msutil.Spectrum();

        // ID
        String id = jmzMLSpec.getId();
        spec.setID(id);

        // scan number
        String[] idToken = id.split("\\s+");
        if (idToken.length > 0 && idToken[idToken.length - 1].matches("scan=\\d+")) {
            int scanNum = Integer.parseInt(idToken[idToken.length - 1].substring(5));
            spec.setScanNum(scanNum);
        }

        // MS Level
        CVParam msLevelParam = null;
        Boolean isCentroided = false;
        for (CVParam cvParam : jmzMLSpec.getCvParam()) {
            if (cvParam.getAccession().equals("MS:1000511"))    // MS level
            {
                msLevelParam = cvParam;
            } else if (cvParam.getAccession().equals("MS:1000127"))    // centroid spectrum
            {
                isCentroided = true;
            } else if (cvParam.getAccession().equals("MS:1000128"))    // profile spectrum
            {
                isCentroided = false;
            }
        }

        spec.setIsCentroided(isCentroided);

        float precursorMz = -1;
        float scanStartTime = -1;
        boolean scanStartTimeIsSeconds = true;

        // Scan list to get monoisotopic m/z
        ScanList scanList = jmzMLSpec.getScanList();
        if (scanList != null && scanList.getCount().intValue() > 0
                && scanList.getScan().size() > 0 && scanList.getScan().get(0).getUserParam().size() > 0) {
            for (CVParam cvParam : scanList.getScan().get(0).getCvParam()) {
                if (cvParam.getAccession().equals("MS:1000016")) {
                    scanStartTime = Float.parseFloat(cvParam.getValue());
                    if (cvParam.getUnitAccession().equals("UO:0000031")) {
                        // in minutes
                        scanStartTimeIsSeconds = false;
                    } else if (cvParam.getUnitAccession().equals("UO:0000010")) {
                        // in seconds
                        scanStartTimeIsSeconds = true;
                    }
                }
            }
            for (UserParam param : scanList.getScan().get(0).getUserParam()) {
                if (param.getName().equals("[Thermo Trailer Extra]Monoisotopic M/Z:")) {
                    precursorMz = Float.parseFloat(param.getValue());
                    break;
                }
            }
        }
        spec.setRt(scanStartTime);
        spec.setRtIsSeconds(scanStartTimeIsSeconds);

        int msLevel = msLevelParam != null ? Integer.parseInt(msLevelParam.getValue()) : 0;
        spec.setMsLevel(msLevel);

        // Precursor
        PrecursorList precursorList = jmzMLSpec.getPrecursorList();
        if (precursorList != null && precursorList.getCount().intValue() > 0 && precursorList.getPrecursor().get(0).getSelectedIonList() != null) {
            Precursor precursor = precursorList.getPrecursor().get(0);    // consider only the first precursor

            ParamGroup isolationWindowParams = precursor.getIsolationWindow();
            if (isolationWindowParams != null && isolationWindowParams.getCvParam() != null) {
                Float isolationWindowTargetMz = null;
                for (CVParam param : isolationWindowParams.getCvParam()) {
                    if (param.getAccession().equals("MS:1000827"))    // selected ion m/z
                    {
                        isolationWindowTargetMz = Float.parseFloat(param.getValue());    // assume that unit is m/z (MS:1000040)
                    }
                }
                spec.setIsolationWindowTargetMz(isolationWindowTargetMz);
            }

            // precursor mz, charge
            int precursorCharge = 0;
            float precursorIntensity = 0;

            ParamGroup paramGroup = precursor.getSelectedIonList().getSelectedIon().get(0);
            for (CVParam param : paramGroup.getCvParam()) {
                if (precursorMz < 0.01 && param.getAccession().equals("MS:1000744"))    // selected ion m/z
                {
                    precursorMz = Float.parseFloat(param.getValue());    // assume that unit is m/z (MS:1000040)
                } else if (param.getAccession().equals("MS:1000041"))    // charge state
                {
                    precursorCharge = Integer.parseInt(param.getValue());
                } else if (param.getAccession().equals("MS:1000042"))    // peak intensity
                {
                    precursorIntensity = Float.parseFloat(param.getValue());
                }
            }

            spec.setPrecursor(new Peak(precursorMz, precursorIntensity, precursorCharge));

            // activation method
            ParamGroup actMethodParams = precursor.getActivation();
            boolean isETD = false;
            for (CVParam param : actMethodParams.getCvParam()) {
                ActivationMethod am = ActivationMethod.getByCV(param.getAccession());
                if (am != null) {
                    if (am == ActivationMethod.ETD) {
                        isETD = true;
                        break;
                    }
                    if (spec.getActivationMethod() == null)
                        spec.setActivationMethod(am);
                }
            }
            if (isETD)
                spec.setActivationMethod(ActivationMethod.ETD);
        }

        // Peak list
        BinaryDataArray mzArray = null, intenArray = null;

        if (jmzMLSpec.getBinaryDataArrayList() != null && jmzMLSpec.getBinaryDataArrayList().getBinaryDataArray() != null) {
            for (BinaryDataArray array : jmzMLSpec.getBinaryDataArrayList().getBinaryDataArray()) {
                if (array.getEncodedLength() == 0)
                    continue;
                // check the cvParams
                for (CVParam param : array.getCvParam()) {
                    if (param.getAccession().equals("MS:1000514")) {
                        mzArray = array;
                        break;
                    }
                    if (param.getAccession().equals("MS:1000515")) {
                        intenArray = array;
                        break;
                    }
                }
                if (mzArray != null && intenArray != null)
                    break;
            }
        }

        if (mzArray != null && intenArray != null) {
            Number mzNumbers[] = mzArray.getBinaryDataAsNumberArray();
            Number intenNumbers[] = intenArray.getBinaryDataAsNumberArray();

            if (mzNumbers.length != intenNumbers.length) {
                System.err.println("Different sizes for m/z and intensity value arrays for spectrum" + jmzMLSpec.getId());
                System.exit(-1);
            }

            for (int i = 0; i < mzNumbers.length; i++)
                spec.add(new Peak(mzNumbers[i].floatValue(), intenNumbers[i].floatValue(), 1));
        }

        // SpecIndex
        spec.setSpecIndex(jmzMLSpec.getIndex() + 1);    // 1-based spectrum index

        // sort peaks by increasing order of m/z
        Collections.sort(spec);

        spec.determineIsCentroided();

        // ScanNum is currently missing
        return spec;
    }

    public static Float getPrecursorMzFromJMzMLSpec(uk.ac.ebi.jmzml.model.mzml.Spectrum jmzMLSpec) {
        PrecursorList precursorList = jmzMLSpec.getPrecursorList();
        if (precursorList != null && precursorList.getCount().intValue() > 0 && precursorList.getPrecursor().get(0).getSelectedIonList() != null) {
            Precursor precursor = precursorList.getPrecursor().get(0);    // consider only the first precursor

            // precursor mz
            float precursorMz = 0;

            ParamGroup paramGroup = precursor.getSelectedIonList().getSelectedIon().get(0);
            for (CVParam param : paramGroup.getCvParam()) {
                if (param.getAccession().equals("MS:1000744"))    // selected ion m/z
                {
                    precursorMz = Float.parseFloat(param.getValue());    // assume that unit is m/z (MS:1000040)
                    return precursorMz;
                }
            }
        }
        return null;
    }
}
