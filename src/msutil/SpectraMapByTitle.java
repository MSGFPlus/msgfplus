package msutil;

import java.util.Hashtable;

import parser.BufferedRandomAccessLineReader;
import parser.SpectrumParserWithTitle;

public class SpectraMapByTitle extends SpectraMap implements SpectrumAccessorByTitle {

	private Hashtable<String,Integer> titleToSpecIndex = null; 	// key: specIndex, value: filePos
	
	public SpectraMapByTitle(String fileName, SpectrumParserWithTitle parser) 
	{
		super(fileName, parser);
		lineReader.seek(0);
		titleToSpecIndex = parser.getTitleToSpecIndexMap(super.lineReader);
	}
	
	public Spectrum getSpectrumByTitle(String title)
	{
		Integer specIndex = titleToSpecIndex.get(title);
		if(specIndex == null)
			return null;
		else
			return super.getSpectrumBySpecIndex(specIndex);
	}

}
