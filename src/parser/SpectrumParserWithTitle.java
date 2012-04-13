package parser;

import java.util.Hashtable;


public interface SpectrumParserWithTitle extends SpectrumParser {
	Hashtable<String, Integer> getTitleToSpecIndexMap(BufferedRandomAccessLineReader lineReader);	// title -> specIndex
}
