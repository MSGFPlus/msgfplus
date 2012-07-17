package jmzparser;

import java.util.Collections;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

import msutil.Peak;

public class JmzSpectrumConverter {
	public static msutil.Spectrum convert(uk.ac.ebi.pride.tools.jmzreader.model.Spectrum jmzSpec)
	{
		msutil.Spectrum spec = new msutil.Spectrum();
		
		Double precursorMz = jmzSpec.getPrecursorMZ();
		Double precursorIntensity = jmzSpec.getPrecursorIntensity();
		Integer precursorCharge = jmzSpec.getPrecursorCharge();

		// setup precursor
		Peak precursor;
		if(precursorMz == null || precursorIntensity == null || precursorCharge == null)
			precursor = null;
		else
			precursor = new Peak(precursorMz.floatValue(), precursorIntensity.floatValue(), precursorCharge);
		spec.setPrecursor(precursor);
		
		// setup peak list
		Map<Double, Double> peakList = jmzSpec.getPeakList();
		Iterator<Entry<Double, Double>> itr = peakList.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<Double, Double> entry = itr.next();
			Double mz = entry.getValue();
			Double intensity = entry.getKey();
			if(mz != null && intensity != null)
				spec.add(new Peak(mz.floatValue(), intensity.floatValue(), 1));
		}

		String id = jmzSpec.getId();
		spec.setID(id);
		
		int msLevel = jmzSpec.getMsLevel();
		spec.setMsLevel(msLevel);

		// should set scanNum, specIndex, activationMethod
		
		// sort peaks by increasing order of m/z
		Collections.sort(spec);
		
		return spec;
	}
	
//	public static msutil.Spectrum convert(uk.ac.ebi.jmzml.model.mzml.Spectrum jmzSpec)
//	{
//		msutil.Spectrum spec = new msutil.Spectrum();
//		
//		Double precursorMz = jmzSpec.getPrecursorMZ();
//		Double precursorIntensity = jmzSpec.getPrecursorIntensity();
//		Integer precursorCharge = jmzSpec.getPrecursorCharge();
//
//		// setup precursor
//		Peak precursor;
//		if(precursorMz == null || precursorIntensity == null || precursorCharge == null)
//			precursor = null;
//		else
//			precursor = new Peak(precursorMz.floatValue(), precursorIntensity.floatValue(), precursorCharge);
//		spec.setPrecursor(precursor);
//		
//		// setup peak list
//		Map<Double, Double> peakList = jmzSpec.getPeakList();
//		Iterator<Entry<Double, Double>> itr = peakList.entrySet().iterator();
//		while(itr.hasNext())
//		{
//			Entry<Double, Double> entry = itr.next();
//			Double mz = entry.getValue();
//			Double intensity = entry.getKey();
//			if(mz != null && intensity != null)
//				spec.add(new Peak(mz.floatValue(), intensity.floatValue(), 1));
//		}
//
//		String id = jmzSpec.getId();
//		spec.setID(id);
//		
//		int msLevel = jmzSpec.getMsLevel();
//		spec.setMsLevel(msLevel);
//
//		// should set scanNum, specIndex, activationMethod
//		
//		// sort peaks by increasing order of m/z
//		Collections.sort(spec);
//		
//		return spec;
//		return null;
//	}
	
}
