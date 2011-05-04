package msutil;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;

import parser.SpectrumParser;

public class SpectraContainer extends ArrayList<Spectrum> {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public SpectraContainer() {}
	public SpectraContainer(String fileName, SpectrumParser parser) 
	{
		SpectraIterator iterator = null;
		try {
			iterator = new SpectraIterator(fileName, parser);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		while(iterator.hasNext())
			this.add(iterator.next());
	}
	
	public void outputMgfFile(String fileName)
	{
		PrintStream out = null; 
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(fileName)));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		for(Spectrum spec : this)
		{
			spec.outputMgf(out);
			out.println();
		}
		out.close();
	}
}
