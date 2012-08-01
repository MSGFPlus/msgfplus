package jmzparser;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Collections;
import java.util.Iterator;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import msutil.ActivationMethod;
import msutil.Peak;

import uk.ac.ebi.jmzml.model.mzml.*;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;


public class MzMLSpectraIterator implements Iterator<msutil.Spectrum>, Iterable<msutil.Spectrum> {
	private MzMLUnmarshaller unmarshaller;
	private MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> itr;
	
	public MzMLSpectraIterator(File specFile) throws FileNotFoundException
	{
		unmarshaller = new MzMLUnmarshaller(specFile);
		itr = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", uk.ac.ebi.jmzml.model.mzml.Spectrum.class);
	}
	
	public boolean hasNext() 
	{
		return itr.hasNext();
	}

	public msutil.Spectrum next() 
	{
		uk.ac.ebi.jmzml.model.mzml.Spectrum jmzSpec = itr.next();
		return getSpectrumFromJMzMLSpec(jmzSpec);
	}

	public void remove() 
	{
		throw new UnsupportedOperationException("SpectraIterator.remove() not implemented");
	}
	
	public Iterator<msutil.Spectrum> iterator() {
		return this;
	}
	
	public static msutil.Spectrum getSpectrumFromJMzMLSpec(uk.ac.ebi.jmzml.model.mzml.Spectrum jmzMLSpec)
	{
		msutil.Spectrum spec = new msutil.Spectrum();

        // ID
		String id = jmzMLSpec.getId();
		spec.setID(id);

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
					break;
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

		// SpecIndex
		spec.setSpecIndex(jmzMLSpec.getIndex()+1);	// 1-based spectrum index
		
		// sort peaks by increasing order of m/z
		Collections.sort(spec);
		
		// ScanNum is currently missing
		return spec;
	}
	
	public static void main(String argv[]) throws Exception
	{
		Logger.getRootLogger().setLevel(Level.OFF);
		test();
	}
	
	public static void test() throws Exception
	{
		File xmlFile = new File("/cygwin/home/kims336/Research/Data/JMzReader/example.mzML");
		xmlFile = new File("/cygwin/home/kims336/Research/Data/JMzReader/small.pwiz.1.1.mzML");
		
		MzMLSpectraIterator itr = new MzMLSpectraIterator(xmlFile);
		while(itr.hasNext())
		{
			msutil.Spectrum spec = itr.next();
			System.out.println("-----------");
			System.out.println(spec.getID());
			System.out.println(spec.getSpecIndex());
			System.out.println(spec.getMSLevel());
			if(spec.getMSLevel() == 2)
			{
				System.out.println(spec.getPrecursorPeak().getMz());
				System.out.println(spec.getPrecursorPeak().getCharge());
				System.out.println(spec.getActivationMethod().getName());
			}
		}
	}
	
}