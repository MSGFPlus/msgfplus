package edu.ucsd.msjava.fdr;

import java.io.File;
import java.io.PrintStream;

public class MzIDPSMSet extends PSMSet {

	private final File mzIDFile;
	private final String scoreName;
	private final boolean isGreaterBetter;
	
	public MzIDPSMSet(File mzIDFile, String scoreName, boolean isGreaterBetter)
	{
		this.mzIDFile = mzIDFile;
		this.scoreName = scoreName;
		this.isGreaterBetter = isGreaterBetter;
	}
	
	public void read()
	{
		
	}
	
	@Override
	public boolean isGreaterBetter() {
		return isGreaterBetter;
	}

	@Override
	public void writeResults(TargetDecoyAnalysis tda, PrintStream out,
			float fdrThreshold, float pepFDRThreshold, float scoreThreshold) {
	}

}
