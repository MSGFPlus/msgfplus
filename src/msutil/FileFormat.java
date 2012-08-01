package msutil;

public class FileFormat {
	public static final FileFormat DIRECTORY = new FileFormat("__DIRECTORY__");
	
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
	
	public FileFormat setCaseSensitive()
	{
		this.isCaseSensitive = true;
		return this;
	}
	
	public boolean isCaseSensitive() {
		return isCaseSensitive;
	}
	
	public String[] getSuffixes() {
		return suffixes;
	}
	
	public String toString()
	{
		if(suffixes == null || suffixes.length == 0)
			return "null";
		StringBuffer buf = new StringBuffer();
		buf.append("["+suffixes[0]);
		for(int i=1; i<suffixes.length; i++)
			buf.append(","+suffixes[i]);
		buf.append("]");
		return buf.toString();
	}
}
