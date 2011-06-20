package misc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import parser.BufferedLineReader;
import parser.PSM;
import parser.PSMList;
import parser.PepXMLParser;

public class MSBlenderTest {
	public static void main(String argv[]) throws Exception
	{
		sequestTest();
//		isMSBlenderBetterThanMSGFDB();
	}
	
	public static void sequestTest() throws Exception
	{
		String xmlFileName = "/Users/sangtaekim/Research/ToolDistribution/Taejoon/20090731_SMPAO1_1_2.sequest.pepxml.pepxml";
		PSMList<PSM> psmList = PepXMLParser.parse(xmlFileName);
		ArrayList<Float> target = new ArrayList<Float>();
		ArrayList<Float> decoy = new ArrayList<Float>();

		System.out.println("NumPSMs: " + psmList.size());
		for(PSM psm : psmList)
		{
			String protein = psm.getProtein();
			float xcorr = psm.getScore("xcorr");
			if(!protein.startsWith("xf_"))
				target.add(xcorr);
			else
				decoy.add(xcorr);
		}
		System.out.println(getNumID(target, decoy));
	}
	
	public static void isMSBlenderBetterThanMSGFDB() throws Exception
	{
		String fileName = "/Users/sangtaekim/Research/ToolDistribution/Taejoon/MSBlender/SMPAO1_1.MSGFDB.msblender_in.msblender_out";
		fileName = "/Users/sangtaekim/Research/ToolDistribution/Taejoon/MSBlender/SMPAO1_1.sequest.msblender_in.msblender_out";
		fileName = "/Users/sangtaekim/Research/ToolDistribution/Taejoon/MSBlender/SMPAO1_1.all.msblender_in.msblender_out";
		fileName = "/Users/sangtaekim/Research/ToolDistribution/Taejoon/MSBlender/SMPAO1_1.no_MSGFDB.msblender_in.msblender_out";
		int scoreCol = 5;
		
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		in.readLine();	// header
		
		HashMap<String, String> best = new HashMap<String, String>();
		
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 4)
				continue;
			String key = token[0].substring(0, token[0].indexOf('.', token[0].indexOf('.')+1));
			float score = Float.parseFloat(token[scoreCol]);
			if(best.get(key) == null)
				best.put(key, s);
			else
			{
				float prevBestScore = Float.parseFloat(best.get(key).split("\t")[scoreCol]);
				if(score > prevBestScore)
					best.put(key, s);
			}
		}
		
		ArrayList<Float> target = new ArrayList<Float>();
		ArrayList<Float> decoy = new ArrayList<Float>();
		
		for(String key : best.keySet())
		{
			s = best.get(key);
			String[] token = s.split("\t");
			float score = Float.parseFloat(token[scoreCol]);
			if(token[1].equalsIgnoreCase("F"))
				target.add(score);
			else
				decoy.add(score);
		}
		
		int numID = getNumID(target, decoy);
		System.out.println(numID);
	}
	
	public static int getNumID(ArrayList<Float> target, ArrayList<Float> decoy)
	{
		TreeMap<Float,Float> fdrMap = fdr.TargetDecoyPSMSet.getFDRMap(target, decoy, true, true, 1);
		
		int numID = 0;
		for(float score : target)
		{
			float psmFDR = fdrMap.lowerEntry(score).getValue();
			if(psmFDR > 0.01f)
				continue;
			numID++;
		}
		return numID;
	}
}
