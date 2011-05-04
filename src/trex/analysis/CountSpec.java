package trex.analysis;

import java.io.IOException;

import parser.BufferedLineReader;

public class CountSpec {
	public static void main(String[] argv) throws IOException{
		BufferedLineReader in = new BufferedLineReader("/home/kwj/workspace/outputs/Mas/deNovo.out");
		
		String s;
		int scanNum = -1, numSpec = 0;
		
		while((s=in.readLine()) != null){
			if(s.startsWith("#") || s.startsWith("Spec Title")) continue;
			
			int sn = Integer.parseInt(s.split("\t")[1]); 
			if(scanNum != sn){
				scanNum = sn;
				numSpec++;
			}
		}
		
		System.out.println(numSpec);
		in.close();
	}
}
