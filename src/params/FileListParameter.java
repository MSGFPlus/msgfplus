package params;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;

public class FileListParameter extends Parameter {

	private HashSet<String> extensions = new HashSet<String>();	// if not empty, only those files with extensions contained here will be considered.
	
	private File[] files;
	
	public FileListParameter(String key, String name, String description) {
		super(key, name, description);
	}

	public FileListParameter setAsOptional()
	{
		super.setOptional();
		return this;
	}
	
	public FileListParameter addExtension(String ext)
	{
		extensions.add(ext.toLowerCase());
		return this;
	}
	
	@Override
	public String parse(String value) 
	{
		File path = new File(value);
		
		if(!path.isDirectory())
			return "must be a directory";
		
		ArrayList<File> fileList = new ArrayList<File>();
		for(File f : path.listFiles())
		{
			if(extensions.isEmpty())
			{
				fileList.add(f);
			}
			else
			{
				String fileName = f.getName();
				String ext = fileName.substring(fileName.lastIndexOf('.')+1);
				if(extensions.contains(ext.toLowerCase()))
				{
					fileList.add(f);
				}
			}
		}
		if(fileList.size() == 0)
		{
			return "file does not exist";
		}
		
		files = fileList.toArray(new File[0]);
		return null;
	}
	
	public File[] getFiles()
	{
		return files;
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
