package edu.ucsd.msjava.msutil;

import java.io.File;

public class DBSearchIOFiles {
	private File specFile;
	private SpecFileFormat specFileFormat;
	private File outputFile;
	
	public DBSearchIOFiles(File specFile, SpecFileFormat specFileFormat, File outputFile) 
	{
		this.specFile = specFile;
		this.specFileFormat = specFileFormat;
		this.outputFile = outputFile;
	}
	
	public File getSpecFile() 
	{
		return specFile;
	}
	
	public SpecFileFormat getSpecFileFormat() 
	{
		return specFileFormat;
	}
	
	public File getOutputFile() 
	{
		return outputFile;
	}
}
