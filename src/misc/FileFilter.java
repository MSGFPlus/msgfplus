package misc;

import java.io.File;

public class FileFilter {
	public static class FileExtFilter implements java.io.FileFilter {
		String ext;
		
		public FileExtFilter(String ext)
		{
			this.ext = ext;
		}
		
		public boolean accept(File f) {
			if(f.getName().endsWith(ext))
				return true;
			return false;
		}
	}
	
	public static class FilePrefixFilter implements java.io.FileFilter {
		String prefix;
		
		public FilePrefixFilter(String ext)
		{
			this.prefix = ext;
		}
		
		public boolean accept(File f) {
			if(f.getName().startsWith(prefix))
				return true;
			return false;
		}

	}
	
}
