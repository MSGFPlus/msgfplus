package msutil;

public class SpecFileFormat extends FileFormat {
	private SpecFileFormat(String[] suffixes)
	{
		super(suffixes);
	}
	
	private SpecFileFormat(String suffix)
	{
		super(suffix);
	}
	
	public static final SpecFileFormat MGF;
	public static final SpecFileFormat MZXML;
	public static final SpecFileFormat MZML;
	public static final SpecFileFormat MS2;
	public static final SpecFileFormat PKL;
	public static final SpecFileFormat DTA_TXT;
	
	static {
		MGF = new SpecFileFormat(".mgf");
		MZXML = new SpecFileFormat(".mzXML");
		MZML = new SpecFileFormat(".mzML");
		MS2 = new SpecFileFormat(".ms2");
		PKL = new SpecFileFormat(".pkl");
		DTA_TXT = new SpecFileFormat("_dta.txt");
	}
}
