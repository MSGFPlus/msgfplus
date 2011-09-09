package parser;

import java.io.FileNotFoundException;

import msutil.SpectraIterator;

public class PNNLSpectraIterator extends SpectraIterator {

	private String scanTypeFileName;
	
	public PNNLSpectraIterator(String fileName, SpectrumParser parser) throws FileNotFoundException 
	{
		super(fileName, parser);
	}

	
}
