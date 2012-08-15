package edu.ucsd.msjava.msutil;

import java.io.FileNotFoundException;
import java.util.Iterator;

import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.LineReader;
import edu.ucsd.msjava.parser.SpectrumParser;


public class SpectraIterator implements Iterator<Spectrum>, Iterable<Spectrum> {
	private String[] filenames=null;
	private String nextSpecFilePath;
	private int nextFileIndex;
	private SpectrumParser parser;
	private boolean hasNext;
	protected Spectrum currentSpectrum;
	LineReader lineReader;
	private int specIndex;		
	
	public SpectraIterator(String fileName, SpectrumParser parser) throws FileNotFoundException
	{
		nextSpecFilePath=fileName;
		lineReader = new BufferedLineReader(fileName);
		this.parser = parser;
		specIndex = 0;
		parseFirstSpectrum();
	}
	
	/**
	 * Added by Louis
	 * Enables iterator to read multiple files seamlessly.
	 * @param filenames List of filenames to process
	 * @param parser
	 * @throws FileNotFoundException thrown only if no files are found
	 */
	public SpectraIterator(String[] filenames, SpectrumParser parser) throws FileNotFoundException
	{
		this.filenames=filenames;
        nextFileIndex=0;
        this.parser=parser;
        specIndex = 0;
        if (!nextFile()) throw new FileNotFoundException("No files found.");
	}
	
	public boolean hasNext() 
	{
		return hasNext;
	}

	public Spectrum next() 
	{
		Spectrum curSpecCopy = currentSpectrum;
		currentSpectrum = parser.readSpectrum(lineReader);
		if(currentSpectrum == null) { // Means file has ended
			if (filenames==null || !nextFile()) hasNext = false;
		}
		else
		{
			currentSpectrum.setSpecIndex(++specIndex);
			currentSpectrum.setID("index="+String.valueOf(specIndex-1));
		}
		return curSpecCopy;
	}

	public void remove() 
	{
		throw new UnsupportedOperationException("SpectraIterator.remove() not implemented");
	}
	
	public Iterator<Spectrum> iterator() {
		return this;
	}
	
	/**
	 * 
	 * @return Filename of source file of next spectrum to be returned by next(). Returns null if last spectrum in last file was returned.
	 */
	private String getNextSpectrumFilePath() {
		return nextSpecFilePath;
	}
	
	private boolean nextFile() {
		lineReader=null;
		nextSpecFilePath=null;
		while(nextFileIndex<filenames.length) {
			try {
				nextSpecFilePath=filenames[nextFileIndex++];
				lineReader = new BufferedLineReader(getNextSpectrumFilePath());
				break;
			} catch (FileNotFoundException e) { 
				// Suppress file not found error - when files in directory has disappeared while reading other files
			}
		}
		if (lineReader==null)
			return false;
		else {
	        parseFirstSpectrum();
	        return true;
		}
	}
	
	private void parseFirstSpectrum()
	{
		currentSpectrum = parser.readSpectrum(lineReader);
		if (currentSpectrum==null) throw new Error("Error while parsing spectrum");
		if(currentSpectrum != null)
		{
			hasNext = true;
			currentSpectrum.setSpecIndex(++specIndex);
			currentSpectrum.setID("index="+String.valueOf(specIndex-1));
		}
		else
			hasNext = false;
	}
}