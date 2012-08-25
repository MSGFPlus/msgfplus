package edu.ucsd.msjava.fdr;

import java.io.File;
import java.io.PrintStream;
import java.util.List;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

public class MzIDPSMSet extends PSMSet {

	private final File mzIDFile;
	private final String scoreName;
	private final boolean isGreaterBetter;
	
	private MzIdentMLUnmarshaller unmarshaller;
	
	public MzIDPSMSet(File mzIDFile, String scoreName, boolean isGreaterBetter)
	{
		this.mzIDFile = mzIDFile;
		this.scoreName = scoreName;
		this.isGreaterBetter = isGreaterBetter;
	}
	
	
	@Override
	public boolean isGreaterBetter() 
	{
		return isGreaterBetter;
	}

	@Override
	public void writeResults(TargetDecoyAnalysis tda, PrintStream out,
			float fdrThreshold, float pepFDRThreshold, float scoreThreshold) {
	}

	@Override
	public void read() 
	{
		unmarshaller = new MzIdentMLUnmarshaller(mzIDFile);
	}

}
