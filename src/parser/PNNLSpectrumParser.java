package parser;

import java.util.Hashtable;

import msutil.Peak;
import msutil.Spectrum;

public class PNNLSpectrumParser implements SpectrumParser {

	public Spectrum readSpectrum(LineReader lineReader) 
	{
		String buf;
		Spectrum spec = new Spectrum();

		// first line: M+H scan=### cs=# scan=### cs=###
		// first line: M+H, charge
		buf = lineReader.readLine();	
		String[] token = buf.split("\\s+");
		assert(token.length > 2);
		float precursorMZ = Float.parseFloat(token[0]);
		int charge = Integer.parseInt(token[1]);
//		int scanNum = Integer.parseInt(token[2].substring(token[2].indexOf('=')+1));
		spec.setPrecursor(new Peak(precursorMZ, 0, charge));
//		spec.setScanNum(scanNum);
		
		while((buf = lineReader.readLine()) != null) 
		{
			if(buf.length() <= 1)
				continue;
			if(Character.isDigit(buf.charAt(0)))
			{
				String[] token2 = buf.split("\\s+");
				if(token2.length != 2)
					continue;
				float mass = Float.parseFloat(token2[0]);
				float intensity = Float.parseFloat(token2[1]);
				spec.add(new Peak(mass, intensity, 1));
			}
			else if(buf.startsWith("="))
			{
				return spec;
			}
		}
		return spec;	
	}

	@Override
	public Hashtable<Integer, Long> getScanNumMap(BufferedRandomAccessLineReader lineReader)
	{
		Hashtable<Integer, Long> scanNumMap = new Hashtable<Integer, Long>();
		String buf;
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.startsWith("=="))
			{
				/*
				int specOffset = lineReader.getPosition();
				buf = lineReader.readLine();
				String[] token = buf.split("\\s+");
				assert(token.length > 4);
				int scanNum = Integer.parseInt(token[2].substring(token[2].indexOf('=')+1));
				*/
				int startIndex = buf.indexOf('.');
				int endIndex = buf.indexOf('.', startIndex+1);
				int scanNum = Integer.parseInt(buf.substring(startIndex+1, endIndex));
				long specOffset = lineReader.getPosition();
				scanNumMap.put(scanNum, specOffset);
			}
		}
		return scanNumMap;
	}
}
