package params;

import java.io.File;
import java.util.ArrayList;

import msutil.FileFormat;

public class FileParameter extends Parameter {

	private boolean mustExist = false;
	private boolean mustNotExist = false;
	
	private boolean mustBeADirectory = false;
	private boolean mustBeAFile = false;
	
	private ArrayList<FileFormat> fileFormats = new ArrayList<FileFormat>();	// available file format; if empty, all files are allowed.
	
	private File file;
	private FileFormat fileFormat;
	
	public FileParameter(String key, String name, String description) {
		super(key, name, description);
	}

	public FileParameter setAsOptional()
	{
		super.setOptional();
		return this;
	}
	
	public FileParameter fileMustExist()
	{
		this.mustExist = true;
		return this;
	}
	
	public FileParameter fileMustNotExist()
	{
		this.mustNotExist = true;
		return this;
	}

	public FileParameter mustBeADirectory()
	{
		this.mustBeADirectory = true;
		return this;
	}

	public FileParameter mustBeAFile()
	{
		this.mustBeAFile = true;
		return this;
	}
	
	public FileParameter addFileFormat(FileFormat fileFormat)
	{
		fileFormats.add(fileFormat);
		return this;
	}
	
	public boolean isSupported(FileFormat fileFormat)
	{
		if(fileFormats == null)
			return false;
		else
			return fileFormats.contains(fileFormat); 
	}
	
	@Override
	public String parse(String value) 
	{
		File path = new File(value);
		
		if(path.isDirectory())
		{
			if(this.mustBeAFile)
				return "must not be a directory";
		}
		else	// path is a file
		{
			if(this.mustBeADirectory)
				return "must be a directory";
		}
		
		if(!fileFormats.isEmpty())
		{
			if(path.isDirectory() && fileFormats.contains(FileFormat.DIRECTORY))
			{
				this.fileFormat = FileFormat.DIRECTORY;
			}
			else
			{
				this.fileFormat = null;
				String fileName = path.getName();
				
				for(FileFormat format : fileFormats)
				{
					if(!format.isCaseSensitive())
						fileName = fileName.toLowerCase();
					for(String suffix : format.getSuffixes())
					{
						if(!format.isCaseSensitive())
							suffix = suffix.toLowerCase();
						if(fileName.endsWith(suffix))
						{
							this.fileFormat = format;
							break;
						}
					}
				}
			}
			if(this.fileFormat == null)
				return "extension does not match";
		}
		
		if(this.mustExist && !path.exists())
			return "file does not exist";
			
		if(this.mustNotExist && path.exists())
			return "file already exists";
		
		this.file = path;
		
		return null;
	}
	
	public File getFile()
	{
		return file;
	}
	
	public FileFormat getFileFormat()
	{
		return fileFormat;
	}

	@Override
	public String getValueAsString() {
		if(file == null)
			return null;
		return file.getPath();
	}
}
