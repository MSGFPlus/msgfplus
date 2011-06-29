package msutil;

import java.io.FileNotFoundException;
import java.util.Iterator;

import parser.BufferedLineReader;
import parser.LineReader;
import parser.SpectrumParser;

public class SpectraIterator implements Iterator<Spectrum>, Iterable<Spectrum> {
	private String[] filenames=null;
	private String nextSpecFilePath;
	private int nextFileIndex;
	private SpectrumParser parser;
	private boolean hasNext;
	protected Spectrum currentSpectrum;
	LineReader lineReader;
	private boolean hasScanNum;
	private int sequentialScanNum;		// this will be used when spectra don't have scan numbers
	
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
        sequentialScanNum = 0;
        if (!nextFile()) throw new FileNotFoundException("No files found.");
	}
	
	/**
	 * 
	 * @return Filename of source file of next spectrum to be returned by next(). Returns null if last spectrum in last file was returned.
	 */
	public String getNextSpectrumFilePath() {
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
	
	public SpectraIterator(String fileName, SpectrumParser parser) throws FileNotFoundException
	{
		nextSpecFilePath=fileName;
		lineReader = new BufferedLineReader(fileName);
		this.parser = parser;
		sequentialScanNum = 0;
		parseFirstSpectrum();
	}
	
	public SpectraIterator(LineReader lineReader, SpectrumParser parser) throws FileNotFoundException
	{
		this.lineReader = lineReader;
		this.parser = parser;
		sequentialScanNum = 0;
		parseFirstSpectrum();
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
		if (hasNext && !hasScanNum) // currentSpectrum possible set to first spec from new file in nextFile()
			currentSpectrum.setScanNum(sequentialScanNum++);
		return curSpecCopy;
	}

	public void remove() 
	{
		throw new UnsupportedOperationException("SpectraIterator.remove() not implemented");
	}
	
	private void parseFirstSpectrum()
	{
		currentSpectrum = parser.readSpectrum(lineReader);
		if (currentSpectrum==null) throw new Error("Error while parsing spectrum");
		if(currentSpectrum.getScanNum() < 0)
		{
			hasScanNum = false;
			currentSpectrum.setScanNum(sequentialScanNum++);
		}
		else
			hasScanNum = true;
		if(currentSpectrum != null)
			hasNext = true;
		else
			hasNext = false;
	}

	@Override
	public Iterator<Spectrum> iterator() {
		return this;
	}
}