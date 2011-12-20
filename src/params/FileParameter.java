package params;

import java.io.File;
import java.util.HashSet;

public class FileParameter extends Parameter {

	private boolean mustExist = false;
	private boolean mustNotExist = false;
	
	private boolean mustBeADirectory = false;
	private boolean mustBeAFile = false;
	
	private HashSet<String> extensions = new HashSet<String>();	// if not empty, only those files with extensions contained here will be considered.
	
	private File file;
	
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
	
	public FileParameter addExtension(String ext)
	{
		extensions.add(ext.toLowerCase());
		return this;
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
		
		
		if(!extensions.isEmpty())
		{
			String fileName = path.getName();
			String ext = fileName.substring(fileName.lastIndexOf('.')+1);
			if(!extensions.contains(ext.toLowerCase()))
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

	@Override
	public String getValueAsString() {
		if(file == null)
			return null;
		return file.getPath();
	}
}
