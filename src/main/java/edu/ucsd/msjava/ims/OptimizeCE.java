package edu.ucsd.msjava.ims;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import edu.ucsd.msjava.parser.BufferedLineReader;

public class OptimizeCE {
	public static void main(String argv[]) throws Exception
	{
		makeSummaryTable();
	}
	
	public static void makeSummaryTable() throws Exception
	{
		File dir = new File("/Users/kims336/Research/Data/IMS/Sarc_DTAs");
		
		int[] ceArr = {18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 74};

		class ID {
			public ID(float specProb, int frameNum, int fromScan, int toScan) {
				this.specProb = specProb;
				this.frameNum = frameNum;
				this.fromScan = fromScan;
				this.toScan = toScan;
			}
			float specProb = 1;
			int frameNum = -1;
			int fromScan = -1;
			int toScan = -1;
		}
		
		Map<String,ID[]> table = new LinkedHashMap<String,ID[]>();
		
		// PrecursorMz\tPrevSpecProb
		Map<String,String> metaInfoTable = new HashMap<String,String>();
		
		for(int i=0; i<ceArr.length; i++)
		{
			int ce = ceArr[i];
			String fileName = "imsDtaCE"+ce+"_msgfOutput.txt";
			File result = new File(dir.getPath()+File.separator+fileName);
			if(!result.exists())
			{
				System.out.println(result.getName() + " doesn't exist!");
				System.exit(-1);
			}
			
			String s;
			BufferedLineReader in = new BufferedLineReader(result.getPath());
			in.readLine();	// header
			while((s=in.readLine()) != null)
			{
				String[] token = s.split("\t");
				if(token.length != 10)
					continue;
				String annotation = token[2];
				float precursorMz = Float.parseFloat(token[3]);
				String charge = token[4];
				int frameNum = Integer.parseInt(token[5]);
				int fromScan = Integer.parseInt(token[6]);
				int toScan = Integer.parseInt(token[7]);
				String prevSpecProb = token[8];
				float specProb = Float.parseFloat(token[9]);
				String key = annotation+":"+charge;
				
				ID[] idList = table.get(key);
				if(idList == null)
				{
					idList = new ID[ceArr.length];
					table.put(key, idList);
					metaInfoTable.put(key, precursorMz+"\t"+prevSpecProb);
				}
				idList[i] = new ID(specProb, frameNum, fromScan, toScan);
			}
			
			in.close();
		}
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(dir.getPath()+File.separator+"imsDtaCE_Summary.tsv")));
		
		out.print("Annotation\tCharge\tPrecursorMz\tPrevSpecProb");
		for(int i=0; i<ceArr.length; i++)
			out.print("\tSpecProbCE"+ceArr[i]+"\tFrameCE"+ceArr[i]+"\tFromScanCE"+ceArr[i]+"\tToScanCE"+ceArr[i]);
		out.println("\tBestSpecProb\tBestEnergy");
		for(String key : table.keySet())
		{
			String token[] = key.split(":");
			String annotation = token[0];
			int charge = Integer.parseInt(token[1]);
			String metaInfo = metaInfoTable.get(key);
			out.print(annotation+"\t"+charge+"\t"+metaInfo);
			
			ID[] idArr = table.get(key);
			int bestIDIndex = -1;
			for(int i=0; i<idArr.length; i++)
			{
				ID id = idArr[i];
				if(id != null)
				{
					out.print("\t"+id.specProb+"\t"+id.frameNum+"\t"+id.fromScan+"\t"+id.toScan);
					if(bestIDIndex == -1 || id.specProb < idArr[bestIDIndex].specProb)
						bestIDIndex = i;
				}
				else
					out.print("\t"+1+"\t"+-1+"\t"+-1+"\t"+-1);
			}
			
			ID bestID = idArr[bestIDIndex];
			out.print("\t"+bestID.specProb+"\t"+ceArr[bestIDIndex]);
			out.println();
		}
		
		out.close();
		
		System.out.println("Done");
	}
}
