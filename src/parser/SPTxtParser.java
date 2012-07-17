package parser;

import java.util.Hashtable;

import msutil.AminoAcidSet;
import msutil.Pair;
import msutil.Peak;
import msutil.Peptide;
import msutil.SpectraMap;
import msutil.SpectraMapByTitle;
import msutil.Spectrum;

public class SPTxtParser implements SpectrumParserWithTitle {

	public Spectrum readSpectrum(LineReader lineReader) 
	{
		Spectrum spec = null;

		String buf;
		
		buf = lineReader.readLine();	// Name: n[43]GAAA....MAR/1
		String[] nameToken = buf.split("\\s+");
		String name = nameToken[1];
		Pair<String,Integer> namePair = parseSPTXTName(name);
		
		String pepSeq = namePair.getFirst();
		int precursorCharge = namePair.getSecond();
		
        spec = new Spectrum();
        Peptide pep = new Peptide(pepSeq, AminoAcidSet.getStandardAminoAcidSet());
        spec.setAnnotation(pep);
        spec.setTitle(namePair.getFirst()+":"+namePair.getSecond());

		float precursorMz = 0;
		boolean parse = false;  
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.startsWith("NumPeaks:")) {  
                parse = true;
            }
			else if(buf.startsWith("PrecursorMZ"))
			{
				String[] token = buf.split("\\s+");
				precursorMz = Float.parseFloat(token[1]);
			}	
  			else if(buf.trim().length() == 0)
  			{
  				assert(spec != null);
  				spec.setPrecursor(new Peak(precursorMz, 0, precursorCharge));
  				return spec;
  			}
			else if(parse && Character.isDigit(buf.charAt(0)))
			{
				String[] token = buf.split("\\s+");
				if(token.length < 2)
					continue;
				float mass = Float.parseFloat(token[0]);
				float intensity = Float.parseFloat(token[1]);
				spec.add(new Peak(mass, intensity, 1));
			}
		}
		return null;	
	}

	public Hashtable<Integer, Long> getSpecIndexMap(BufferedRandomAccessLineReader lineReader) 
	{
		Hashtable<Integer, Long> specIndexMap = new Hashtable<Integer, Long>();
		String buf;
		long offset = 0;
		long specOffset = 0;
		int specIndex = 0;
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.startsWith("Name:"))
			{
				specIndex++;
				specOffset = offset;
				specIndexMap.put(specIndex, specOffset);
			}
			offset = lineReader.getPosition();
		}
		return specIndexMap;	
	}
	
	public Hashtable<String, Integer> getTitleToSpecIndexMap(BufferedRandomAccessLineReader lineReader) 
	{
		Hashtable<String, Integer> titleToSpecIndexMap = new Hashtable<String, Integer>();
		String buf;
		int specIndex = 0;
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.startsWith("Name:"))
			{
				specIndex++;
				Pair<String,Integer> pair = parseSPTXTName(buf.split("\\s+")[1]);
				titleToSpecIndexMap.put(pair.getFirst()+":"+pair.getSecond(), specIndex);
			}
		}
		return titleToSpecIndexMap;	
	}
	
	public static Pair<String,Integer> parseSPTXTName(String name)
	{
		String annotationStr = name.substring(0, name.lastIndexOf('/'));
		StringBuffer pepBuf = new StringBuffer();
		int startIndex=0;
		if(annotationStr.startsWith("n[43]"))
		{
			pepBuf.append("+42");
			startIndex = 5;
		}
		char prevAA = '\0';
		for(int i=startIndex; i<annotationStr.length(); i++)
		{
			char c = annotationStr.charAt(i);
			if(Character.isUpperCase(c))
				pepBuf.append(c);
			else if(c == '[')
			{
				StringBuffer massBuf = new StringBuffer();
				while(annotationStr.charAt(++i) != ']')
					massBuf.append(annotationStr.charAt(i));
				int mass = Integer.parseInt(massBuf.toString());
				int residueMass = AminoAcidSet.getStandardAminoAcidSet().getAminoAcid(prevAA).getNominalMass();
				int delMass = mass-residueMass;
				if(delMass > 0)
					pepBuf.append("+");
				pepBuf.append(delMass);
			}
			prevAA = c;
		}		
		
		int charge = Integer.parseInt(name.substring(name.lastIndexOf('/')+1));
		
		return new Pair<String,Integer>(pepBuf.toString(), charge);
	}	
	
	public static void main(String argv[]) throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/NISTLib/human_targetdecoy_spectrast.sptxt";
		SpectraMapByTitle map = new SpectraMapByTitle(fileName, new SPTxtParser());
		System.out.println("Parsing complete.");
		Spectrum spec = map.getSpectrumByTitle("+42AAAAAAGAGPEM+16VRGQVFDVGPR:3");
		System.out.println(spec.getSpecIndex()+"\t"+spec.size());
		
	}
}
