package ims;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import org.junit.Test;

import edu.ucsd.msjava.parser.BufferedLineReader;
public class IMSResultProcessor {
	@Test
	public void filterIMSResults() throws Exception
	{
		File result = new File("C:\\cygwin\\home\\kims336\\Data\\IMS_Sarc\\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Precursors_Removed_Collision_Energy_Collapsed_0801.tsv");
		File output = new File("C:\\cygwin\\home\\kims336\\Data\\IMS_Sarc\\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Precursors_Removed_Collision_Energy_Collapsed_0801_filtered.tsv");

		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(output)));
		
		BufferedLineReader in = new BufferedLineReader(result.getPath());
		
		// header
		out.println(in.readLine());
		
		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 14)
				continue;
			double score = Double.parseDouble(token[12]);
			if(score > -20)
				out.println(s);
		}
		
		in.close();
		out.close();
		
		System.out.println("Done");
	}
}
