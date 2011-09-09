package parser;

import java.util.HashSet;
import java.util.Hashtable;

import msutil.Composition;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.SpectraMap;
import msutil.Spectrum;

public class PNNLSpectrumParser implements SpectrumParser {

	public Spectrum readSpectrum(LineReader lineReader) 
	{
		Spectrum spec = null;

		String buf;
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.length() == 0)
			{
				if(spec != null)
					return spec;
				else
					continue;
			}
			else if(buf.startsWith("=="))
			{
				if(spec != null)
				{
					System.out.println("There must be at least one empty line between spectra: " + buf);
					System.exit(-1);
				}
				int lastDotIndex = buf.lastIndexOf('.');
				int secondLastDotIndex = buf.lastIndexOf('.', lastDotIndex-1);
				int thirdLastDotIndex = buf.lastIndexOf('.', secondLastDotIndex-1);
				int fourthLastDotIndex = buf.lastIndexOf('.', thirdLastDotIndex-1);
				
				int scanNum = Integer.parseInt(buf.substring(fourthLastDotIndex+1, thirdLastDotIndex));
				
				String annotation = buf;
				// first line of a spectrum
				buf = lineReader.readLine();
				if(buf == null || buf.trim().length() == 0)
				{
					System.out.println("Error while parsing _Dta.txt file: " + annotation);
					System.out.println("No spectrum!");
					System.exit(-1);
				}
				
				spec = new Spectrum();
				String[] token = buf.split("\\s+");
				float mPlusH = Float.parseFloat(token[0]);
				int charge = Integer.parseInt(token[1].substring(token[1].indexOf('=')+1));
				float precursorMz = (mPlusH-(float)Composition.PROTON)/charge+(float)Composition.PROTON;
				spec.setPrecursor(new Peak(precursorMz, 0, charge));
				spec.setScanNum(scanNum);
			}
			else if(Character.isDigit(buf.charAt(0)))	// peak
			{
				if(spec == null)
				{
					System.out.println("Error while parsing _Dta.txt file.");
					System.out.println("Header line is missing: " + buf);
					System.exit(-1);
				}
				String[] token2 = buf.split("\\s+");
				if(token2.length != 2)
					continue;
				float mass = Float.parseFloat(token2[0]);
				float intensity = Float.parseFloat(token2[1]);
				spec.add(new Peak(mass, intensity, 1));
			}
		}
		return spec;
	}

	@Override
	public Hashtable<Integer, Long> getSpecIndexMap(BufferedRandomAccessLineReader lineReader)
	{
		Hashtable<Integer, Long> specIndexMap = new Hashtable<Integer, Long>();
		String buf;
		long offset = 0;
		int specIndex = 0;
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.startsWith("=="))
			{
				specIndexMap.put(++specIndex, offset);
			}
			offset = lineReader.getPosition();
		}
		return specIndexMap;
	}
	
	public static void main(String argv[]) throws Exception
	{
		long time = System.currentTimeMillis();
		String fileName = System.getProperty("user.home")+"/Research/ToolDistribution/PNNLTest/QC_Shew_08_04_pt5_b_22Jan09_Owl_09-01-04_dta.txt";
		SpectraIterator itr = new SpectraIterator(fileName, new PNNLSpectrumParser());
		int numSpecs = 0;
		HashSet<Integer> scanNumSet = new HashSet<Integer>();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			numSpecs++;
			if(scanNumSet.contains(spec.getScanNum()))
			{
				System.out.println(spec.getScanNum());
			}
			else
				scanNumSet.add(spec.getScanNum());
//			System.out.println(spec+ "\t" + spec.getScanNum()+"\t"+(spec.getParentMass()+(float)Composition.PROTON)+"\t"+spec.getCharge());
		}
		System.out.println("NumSpecs: " + numSpecs);
		System.out.println("Time: " + (System.currentTimeMillis()-time));
		
		time = System.currentTimeMillis();
		SpectraMap map = new SpectraMap(fileName, new PNNLSpectrumParser());
		numSpecs = 0;
		for(int specIndex : map.getSpecIndexList())
		{
			Spectrum spec = map.getSpectrumBySpecIndex(specIndex);
			numSpecs++;
//			System.out.println(spec+ "\t" + spec.getScanNum()+"\t"+(spec.getParentMass()+(float)Composition.PROTON)+"\t"+spec.getCharge());
		}
		System.out.println("NumSpecs: " + numSpecs);
		System.out.println("Time: " + (System.currentTimeMillis()-time));
	}
}
