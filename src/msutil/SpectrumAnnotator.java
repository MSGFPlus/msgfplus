package msutil;

import java.util.ArrayList;

import msgf.Tolerance;
import msutil.IonType.PrefixIon;

/**
 * For a given peptide spectrum pair, retrieve annotation information.
 * @author sangtaekim
 *
 */

public class SpectrumAnnotator {

	Spectrum spec;
	Peptide peptide;
	/**
	 * Constructor.
	 * @param spec spectrum
	 * @param peptide peptide
	 */
	public SpectrumAnnotator(Spectrum spec, Peptide peptide)
	{
		this.spec = spec;
		this.peptide = peptide;
	}
	
	/**
	 * Gets list of peaks annotated as given ion type.
	 * @param ion ion type
	 * @return peaks with given ion type, empty list if no peak exists
	 */
	public ArrayList<Peak> getPeakListOfIon(IonType ionType, Tolerance tolerance)
	{
		ArrayList<Peak> peakList = new ArrayList<Peak>();
		float prm = 0;
		float prmSum = peptide.getMass();
		for(AminoAcid aa : peptide)
		{
			prm += aa.getMass();
			float mass;
			if(ionType instanceof PrefixIon)
				mass = ionType.getMz(prm);
			else
				mass = ionType.getMz(prmSum - prm);
			Peak p = spec.getPeakByMass(mass, tolerance);
			if(p != null)
				peakList.add(p);
		}
		return peakList;
	}
	
	/**
	 * Gets list of peak errors of given ion type.
	 * @param ion ion type
	 * @return peaks with given ion type, empty list if no peak exists
	 */
	public ArrayList<Float> getPeakErrorPPMOfIon(IonType ionType, Tolerance tolerance)
	{
		ArrayList<Float> errorList = new ArrayList<Float>();
		float prm = 0;
		float prmSum = peptide.getMass();
		for(AminoAcid aa : peptide)
		{
			prm += aa.getMass();
			float mass;
			if(ionType instanceof PrefixIon)
				mass = ionType.getMz(prm);
			else
				mass = ionType.getMz(prmSum - prm);
			Peak p = spec.getPeakByMass(mass, tolerance);
			
			if(p != null)
				errorList.add((p.getMz() - mass)*1e6f/mass);
		}
		return errorList;
	}
	
	/**
	 * Gets list of peak errors of given ion type.
	 * @param ion ion type
	 * @return peaks with given ion type, empty list if no peak exists
	 */
	public ArrayList<Float> getPeakErrorOfIon(IonType ionType, Tolerance tolerance)
	{
		ArrayList<Float> errorList = new ArrayList<Float>();
		float prm = 0;
		float prmSum = peptide.getMass();
		for(AminoAcid aa : peptide)
		{
			prm += aa.getMass();
			float mass;
			if(ionType instanceof PrefixIon)
				mass = ionType.getMz(prm);
			else
				mass = ionType.getMz(prmSum - prm);
			Peak p = spec.getPeakByMass(mass, tolerance);
			
			if(p != null)
				errorList.add(p.getMz() - mass);
		}
		return errorList;
	}

	
	
	/**
	 * Gets the indices of cleavages where peaks are observed
	 * @param ionTypes array of ion types
	 * @return indices of cleavages
	 */
	public ArrayList<Integer> getCoveredCleavages(IonType[] ionTypes, Tolerance tolerance)
	{
		float prm = 0;
		float prmSum = peptide.getMass();
		ArrayList<Integer> indices = new ArrayList<Integer>();
		for(int i=0; i<peptide.size()-1; i++)
		{
			prm += peptide.get(i).getMass();
			for(IonType ionType : ionTypes)
			{
				float mass;
				if(ionType instanceof PrefixIon)
					mass = ionType.getMz(prm);
				else
					mass = ionType.getMz(prmSum - prm);
				Peak p = spec.getPeakByMass(mass, tolerance);
				if(p != null)
				{
					indices.add(i+1);
					break;
				}
			}
		}
		return indices;
	}		
	/**
	 * Gets the coverage of the peptide
	 * @param ionTypes array of ion types
	 * @return coverage of given ion types
	 */
	public float getCoverage(IonType[] ionTypes, Tolerance tolerance)
	{
		return getCoveredCleavages(ionTypes, tolerance).size()/(float)(peptide.size()-1);
	}	
	
}
