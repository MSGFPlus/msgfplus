package msutil;

public class FileFormat {
	private final String[] suffixes;
	private boolean isCaseSensitive = false;
	
	public FileFormat(String[] suffixes)
	{
		this.suffixes = suffixes;
	}
	
	public FileFormat(String suffix)
	{
		this.suffixes = new String[1];
		suffixes[0] = suffix;
	}
	
	protected void setCaseSensitive()
	{
		this.isCaseSensitive = true;
	}
	
	public boolean isCaseSensitive() {
		return isCaseSensitive;
	}
	
	public String[] getSuffixes() {
		return suffixes;
	}
}
