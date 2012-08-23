package edu.ucsd.msjava.mzml;

import java.util.Collections;

import uk.ac.ebi.jmzml.model.mzml.BinaryDataArray;
import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.ParamGroup;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.PrecursorList;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Peak;

public class SpectrumConverter {
	public static edu.ucsd.msjava.msutil.Spectrum getSpectrumFromJMzMLSpec(uk.ac.ebi.jmzml.model.mzml.Spectrum jmzMLSpec)
	{
		edu.ucsd.msjava.msutil.Spectrum spec = new edu.ucsd.msjava.msutil.Spectrum();

        // ID
		String id = jmzMLSpec.getId();
		spec.setID(id);

        // MS Level
		CVParam msLevelParam = null;
		for(CVParam cvParam : jmzMLSpec.getCvParam())
		{
			if(cvParam.getAccession().equals("MS:1000511"))	// MS level
			{
				msLevelParam = cvParam;
				break;
			}
		}
		int msLevel = msLevelParam != null ? Integer.parseInt(msLevelParam.getValue()) : 0;
		spec.setMsLevel(msLevel);
		
		// Precursor
		PrecursorList precursorList = jmzMLSpec.getPrecursorList();
		if(precursorList != null && precursorList.getCount().intValue() > 0 && precursorList.getPrecursor().get(0).getSelectedIonList() != null)
		{
			Precursor precursor = precursorList.getPrecursor().get(0);	// consider only the first precursor
			
			// precursor mz, charge
			float precursorMz = 0;
			int precursorCharge = 0;
			float precursorIntensity = 0;
			
			ParamGroup paramGroup = precursor.getSelectedIonList().getSelectedIon().get(0);
			for(CVParam param : paramGroup.getCvParam())
			{
				if(param.getAccession().equals("MS:1000744"))	// selected ion m/z
				{
					precursorMz = Float.parseFloat(param.getValue());	// assume that unit is m/z (MS:1000040)
				}
				else if(param.getAccession().equals("MS:1000041"))	// charge state
				{
					precursorCharge = Integer.parseInt(param.getValue());
				}
				else if(param.getAccession().equals("MS:1000042"))	// peak intensity
				{
					precursorIntensity = Float.parseFloat(param.getValue());
				}	//MS:1000511
			}

			spec.setPrecursor(new Peak(precursorMz, precursorIntensity, precursorCharge));
			
			// activation method
			ParamGroup actMethodParams = precursor.getActivation();
			for(CVParam param : actMethodParams.getCvParam())
			{
				ActivationMethod am = ActivationMethod.getByCV(param.getAccession());
				if(am != null)
				{
					spec.setActivationMethod(am);
					break;
				}
			}
		}

		// Peak list
        BinaryDataArray mzArray = null, intenArray = null;

        for (BinaryDataArray array : jmzMLSpec.getBinaryDataArrayList().getBinaryDataArray()) {
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

        if (mzArray != null && intenArray != null)
        {
            Number mzNumbers[] = mzArray.getBinaryDataAsNumberArray();
            Number intenNumbers[] = intenArray.getBinaryDataAsNumberArray();
            
            if(mzNumbers.length != intenNumbers.length)
            {
            	System.err.println("Different sizes for m/z and intensity value arrays for spectrum" + jmzMLSpec.getId());
            	System.exit(-1);
            }
            
            for(int i=0; i<mzNumbers.length; i++)
            	spec.add(new Peak(mzNumbers[i].floatValue(), intenNumbers[i].floatValue(), 1));
        }
		
		// SpecIndex
		spec.setSpecIndex(jmzMLSpec.getIndex()+1);	// 1-based spectrum index
		
		// sort peaks by increasing order of m/z
		Collections.sort(spec);
		
		// ScanNum is currently missing
		return spec;
	}
	
	public static Float getPrecursorMzFromJMzMLSpec(uk.ac.ebi.jmzml.model.mzml.Spectrum jmzMLSpec)
	{
		PrecursorList precursorList = jmzMLSpec.getPrecursorList();
		if(precursorList != null && precursorList.getCount().intValue() > 0 && precursorList.getPrecursor().get(0).getSelectedIonList() != null)
		{
			Precursor precursor = precursorList.getPrecursor().get(0);	// consider only the first precursor
			
			// precursor mz
			float precursorMz = 0;
			
			ParamGroup paramGroup = precursor.getSelectedIonList().getSelectedIon().get(0);
			for(CVParam param : paramGroup.getCvParam())
			{
				if(param.getAccession().equals("MS:1000744"))	// selected ion m/z
				{
					precursorMz = Float.parseFloat(param.getValue());	// assume that unit is m/z (MS:1000040)
					return precursorMz;
				}
			}
		}
		return null;
	}	
}
