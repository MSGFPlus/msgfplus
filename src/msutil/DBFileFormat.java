package msutil;

public class DBFileFormat extends FileFormat {
	private DBFileFormat(String[] suffixes)
	{
		super(suffixes);
	}
	
	private DBFileFormat(String suffix)
	{
		super(suffix);
	}

	public static final DBFileFormat FASTA = new DBFileFormat(new String[] {".fa", ".fasta"});
}
