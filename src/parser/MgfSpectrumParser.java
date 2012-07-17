package parser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Peak;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;

/**
 * This class enables to parse spectrum file with mgf format.
 * @author sangtaekim
 *
 */
public class MgfSpectrumParser implements SpectrumParser {

	/**
	 * Amino acid set to be used to parse "SEQ="
	 */
	private AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
	
	/**
	 * Specify amino acid set to be used to parse "SEQ=" field.
	 * @param aaSet amino acid set.
	 * @return this object.
	 */
	public MgfSpectrumParser aaSet(AminoAcidSet aaSet)
	{
		this.aaSet = aaSet;
		return this;
	}
	/**
	 * Implementation of readSpectrum method. Implicitly lineReader points to the start of a spectrum. 
	 * Reads mgf file line by line until the spectrum ends, generate a Spectrum object and returns it. 
	 * If it cannot read a spectrum, it returns null.
	 * @param lineReader a LineReader object points to the start of a spectrum
	 * @return a spectrum object. null if no spectrum can be generated.
	 */
	public Spectrum readSpectrum(LineReader lineReader)
	{
		Spectrum spec = null;
		String title = null;

		float precursorMz = 0;
		float precursorIntensity = 0;
		int precursorCharge = 0;
		ActivationMethod activation = null;

		String buf;
		boolean parse = false;   // parse only after the BEGIN IONS
		boolean sorted = true;
		float prevMass = 0;
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.length() == 0)
				continue;
			
			if(buf.startsWith("BEGIN IONS")) {  
                parse = true;
                spec = new Spectrum();
            }
			else if(parse){
                if(Character.isDigit(buf.charAt(0)))
                {
                    assert(spec != null);
                    String[] token = buf.split("\\s+");
                    if(token.length < 2)
                        continue;
                    float mass = Float.parseFloat(token[0]);
                    if(sorted && mass < prevMass)
                    	sorted = false;
                    else
                    	prevMass = mass;
                    float intensity = Float.parseFloat(token[1]);
                    spec.add(new Peak(mass, intensity, 1));
                }
  			else if(buf.startsWith("TITLE"))
  			{
  				title = buf.substring(buf.indexOf('=')+1);
  				spec.setTitle(title);
  			}
  			else if(buf.startsWith("CHARGE"))
  			{
  				String chargeStr = buf.substring(buf.indexOf("=")+1).trim();
  				if(chargeStr.startsWith("+"))
  					chargeStr = chargeStr.substring(1);
  				if(chargeStr.charAt(chargeStr.length()-1) == '+')
  					chargeStr = chargeStr.substring(0, chargeStr.length()-1);
  				precursorCharge = Integer.valueOf(chargeStr);
  				if(precursorCharge == 0)
  					precursorCharge = 2;
  			}
  			else if(buf.startsWith("SEQ"))
  			{
  				String annotationStr = buf.substring(buf.lastIndexOf('=')+1);
  				if(spec.getAnnotation() == null)
  	  				spec.setAnnotation(new Peptide(annotationStr, aaSet));
  				spec.addSEQ(annotationStr);
  			}
  			else if(buf.startsWith("PEPMASS"))
  			{
  				String[] token = buf.substring(buf.indexOf("=")+1).split("\\s+");
  				precursorMz = Float.valueOf(token[0]);
  			}
  			else if(buf.startsWith("SCANS"))
  			{
  				if(buf.matches(".+=\\d+-\\d+"))	// e.g. SCANS=953-959
  				{
  					int startScanNum = Integer.parseInt(buf.substring(buf.indexOf('=')+1, buf.lastIndexOf('-')));
  					int endScanNum = Integer.parseInt(buf.substring(buf.lastIndexOf('-')+1));
  					spec.setStartScanNum(startScanNum);
  					spec.setEndScanNum(endScanNum);
  				}
  				else
  				{
  	  				// for mgf files, scan number is set as the zero based sequence number of the spectrum
  	  				int scanNum = Integer.valueOf(buf.substring(buf.indexOf("=")+1));
  	  				spec.setScanNum(scanNum);
  				}
  			}
  			else if(buf.startsWith("ACTIVATION"))
  			{
  				String activationName = buf.substring(buf.indexOf("=")+1);
  				activation = ActivationMethod.get(activationName);
  				spec.setActivationMethod(activation);
  			}
  			else if(buf.startsWith("END IONS"))
  			{
  				assert(spec != null);
  				spec.setPrecursor(new Peak(precursorMz, precursorIntensity, precursorCharge));
  				if(!sorted)
  					Collections.sort(spec);
  				return spec;
  			}
		  }
		}
		return null;		
	}	

	/**
	 * Implementation of getSpecIndexMap object. Reads the entire spectrum file and
	 * generates a map from a spectrum index of a spectrum and the position of the spectrum.
	 * @param lineReader a LineReader object that points to the start of a file.
	 * @return A hashtable maps spectrum indexes to spectral positions in the file.
	 */
	public Hashtable<Integer, Long> getSpecIndexMap(BufferedRandomAccessLineReader lineReader)
	{
		Hashtable<Integer, Long> specIndexMap = new Hashtable<Integer, Long>();
		String buf;
		long offset = 0;
		long specOffset = 0;
		int specIndex = 0;
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.startsWith("BEGIN IONS"))
			{
				specIndex++;
				specOffset = offset;
				specIndexMap.put(specIndex, specOffset);
			}
			offset = lineReader.getPosition();
		}
		return specIndexMap;
	}
	
	// test code
	public static void main(String argv[]) throws Exception
	{
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
	    while(iterator.hasNext())
	    	specContainer.add(iterator.next());
		System.out.println("Time: " + (System.currentTimeMillis() - time));
	}
}
