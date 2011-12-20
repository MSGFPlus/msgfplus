package msutil;

public class FileFormat {
	private final String[] suffixes;
	private boolean isCaseSensitive = false;
	
	protected FileFormat(String[] suffixes)
	{
		this.suffixes = suffixes;
	}
	
	protected FileFormat(String suffix)
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
