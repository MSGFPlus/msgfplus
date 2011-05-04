package msgap;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import parser.BufferedLineReader;

public class CountAAinFasta {
	public static long countAA(String filename) throws IOException{
		File folder = new File(filename);
		File[] files;
		long totalAA = 0;
		
		if(folder.isDirectory()){
			files = folder.listFiles();
		}else{
			files = new File[1];
			files[0] = folder; 
		}
		
		for(File file : files){
			if(!file.getName().endsWith(".fasta") && !file.getName().endsWith(".fa")) continue;
			BufferedLineReader in = new BufferedLineReader(file.getPath());
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith(">")) continue;
				else totalAA += s.length();
			}
			in.close();
		}
		
		return (totalAA);
	}
}
