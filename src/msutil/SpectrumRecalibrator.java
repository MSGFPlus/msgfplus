package msutil;

import msgf.LinearCalibration;
import msgf.Tolerance;

/**
 * This class creates an object that re-calibrates Spectrum objects. 
 * Currently, it tries to find alpha and beta that converts a mass m into alpha*m + beta.
 * It uses the least squares method.
 * @author sangtaekim
 *
 */
public class SpectrumRecalibrator implements Reshape {
	Tolerance tolerance;
	/**
	 * Constructor.
	 * @param tolerance	Tolerance.
	 */
	public SpectrumRecalibrator(Tolerance tolerance)
	{
		this.tolerance = tolerance;
	}

	public SpectrumRecalibrator()
	{
		this.tolerance = null;
	}
	/**
	 * @param s	Spectrum
	 * @return	Re-calibrated spectrum 
	 */
	public Spectrum apply(Spectrum s) {
		return null;
	}

	/**
	 * Use rescaled masses for recalibration.
	 */
	public Spectrum recalibrateUsingRescaling(Spectrum s) {
		LinearCalibration calibrator = new LinearCalibration();
		WindowFilter filter = new WindowFilter(5, 50);
		Spectrum filteredSpec = filter.apply(s);
		for(Peak p : filteredSpec)
		{
			float rescaledMass = p.getMz()*Constants.INTEGER_MASS_SCALER;
			calibrator.addData(rescaledMass, Math.round(rescaledMass));
		}
		
		Spectrum newSpec = s.getCloneWithoutPeakList();
		for(Peak p : s)
		{
			newSpec.add(p.getShiftedPeak(calibrator.predict(p.getMz())));
		}
		
		return newSpec;
	}
	
	/**
	 * Re-calibrate the given spectrum using the given annotation as a template.
	 * @param s	Spectrum
	 * @param annotation	Identification of s
	 * @return	Re-calibrated spectrum
	 */
	public Spectrum apply(Spectrum s, Peptide annotation)
	{
		LinearCalibration calibrator = new LinearCalibration();
		
		// extracting b/y peaks
		float prmMass = 0;
		float prmSum = annotation.getMass();
		for(AminoAcid aa : annotation)
		{
			prmMass += aa.getMass();
			
			float bMass = IonType.B.getMz(prmMass);
				
			Peak p = s.getPeakByMass(bMass, tolerance);
			if(p != null)
				calibrator.addData(p.getMz(), bMass);
			float yMass = IonType.Y.getMz(prmSum - prmMass);
			p = s.getPeakByMass(yMass, tolerance);
			if(p != null)
				calibrator.addData(p.getMz(), yMass);
		}
		Spectrum newSpec = s.getCloneWithoutPeakList();
		for(Peak p : s)
		{
			newSpec.add(p.getShiftedPeak(calibrator.predict(p.getMz())));
		}
		
		return newSpec;
	}
}
