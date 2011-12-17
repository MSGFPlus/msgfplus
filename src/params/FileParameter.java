package params;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;

public class FileParameter extends Parameter {

	private boolean mustExist = false;
	private boolean mustNotExist = false;
	private boolean canBeADirectory = false;
	private HashSet<String> extensions = new HashSet<String>();	// if not empty, only those files with extensions contained here will be considered.
	
	private File[] files;
	
	public FileParameter(String key, String name, String description, boolean isOptional) {
		super(key, name, description, isOptional);
	}

	public FileParameter mustExist()
	{
		this.mustExist = true;
		return this;
	}
	
	public FileParameter mustNotExist()
	{
		this.mustNotExist = true;
		return this;
	}

	public FileParameter canBeADirectory()
	{
		this.canBeADirectory = true;
		return this;
	}
	
	public FileParameter addExtension(String ext)
	{
		extensions.add(ext.toLowerCase());
		return this;
	}
	
	@Override
	public boolean parse(String value) 
	{
		File path = new File(value);
		
		File[] filesToBeConsidered;
		if(path.isDirectory())
		{
			if(!this.canBeADirectory)
			{
				return false;
			}
			filesToBeConsidered = path.listFiles();
		}
		else
		{
			filesToBeConsidered = new File[1];
			filesToBeConsidered[0] = path;
		}
		
		ArrayList<File> fileList = new ArrayList<File>();
		for(File f : filesToBeConsidered)
		{
			if(extensions.isEmpty())
			{
				if(this.mustExist && !f.exists())
					continue;
				if(this.mustNotExist && f.exists())
					continue;
				fileList.add(f);
			}
			else
			{
				String fileName = f.getName();
				String ext = fileName.substring(fileName.lastIndexOf('.')+1);
				if(extensions.contains(ext.toLowerCase()))
				{
					if(this.mustExist && !f.exists())
						continue;
					if(this.mustNotExist && f.exists())
						continue;
					fileList.add(f);
				}
			}
		}
		if(fileList.size() == 0)
		{
			return false;
		}
		
		files = fileList.toArray(new File[0]);
		return true;
	}
	
	public File[] getFiles()
	{
		return files;
	}
	
	public File getFile()
	{
		if(this.canBeADirectory || files.length != 1)
			return null;
		return files[0];
	}

	@Override
	public String getValueAsString() {
		if(files == null)
			return null;
		StringBuffer output = new StringBuffer();
		if(files.length == 0)
			return output.toString();
		output.append(files[0].getPath());
		for(int i=1; i<files.length; i++)
			output.append(","+files[i].getPath());
		return output.toString();
	}
}
