package msgap.analysis;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import msgap.ResultProcessor;
import msutil.AminoAcidSet;
import parser.BufferedLineReader;
import parser.InsPecTPSM;
import parser.InsPecTParser;
import parser.MSGappedDictionaryPSM;
import parser.MSGappedDictionaryParser;
import parser.PSMList;

public class CompareWithInsPecT {
	
	static HashSet<String> insPep = new HashSet<String>();
	static HashSet<String> gapPep = new HashSet<String>();
	
	
	private static int getNum(String fn, float threshold) throws IOException{
		BufferedLineReader in = new BufferedLineReader(fn);

		String s;
		int scanBefore = -1 , num = 0, numMP = 0;
		HashSet<String> peptides = new HashSet<String>();
		
		while((s = in.readLine()) != null){
			if(s.startsWith("#")) continue;
			
			String[] token = s.split("\t");
			
			int scanNum = Integer.parseInt(token[1]);
			
			if(scanNum == scanBefore) continue; // single peptide per spectrum
			scanBefore = scanNum;
			if(Float.parseFloat(token[14]) < threshold) continue;
			//System.out.println(scanNum);
			peptides.add(token[2].substring(token[2].indexOf('.')+1, token[2].lastIndexOf('.')));
			
		//	if(token[2].contains("+"))
			//	System.out.println(token[2] + "\t" + token[0] + "\t" + token[1]);
				
				
			num++;
		}
		System.out.println(num + "peps : " + "\t" + peptides.size());
		
		for(String pep : peptides)
			if(pep.contains("+")){ 
				numMP++;
				//insPep.add(pep.replace("M+16", "m"));
			//	System.out.println(pep);
			}
		System.out.println("MP " +  numMP);
		
		
		return num;
	}
	
	public static void main(String[] argv) throws IOException{
		//String targetOutFileName = "/home/kwj/workspace/outputs/PTM/SF/InsPecT/woPTM/out_2.5Da_wo_ptm.txt";
		//String decoyOutFileName = "/home/kwj/workspace/outputs/PTM/SF/InsPecT/woPTM/out_decoy_2.5Da_wo_ptm.txt";
		String targetOutFileName = "/home/kwj/workspace/outputs/cavebear/inspect/out.txt";
		String decoyOutFileName = "/home/kwj/workspace/outputs/cavebear/inspect/out_rev.txt";
		
		
		float threshold = 14f;
		
		getNum(targetOutFileName, threshold);
		getNum(decoyOutFileName, threshold);
		
		//7981peps : 	4825
		//MP 545
		
		
	//	String aaFileName = "/home/kwj/workspace/inputs/AA_shewanella.txt";
		String aaFileName = "/home/kwj/workspace/inputs/AASetWithHydroProline.txt";
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(aaFileName);
	//	targetOutFileName = "/home/kwj/workspace/outputs/PTM/SF/MSGD/out.txt";
		targetOutFileName = "/home/kwj/workspace/outputs/cavebear/new.txt_before_filter";
		decoyOutFileName = "/home/kwj/workspace/outputs/cavebear/new.txt.decoy_out_before_filter";
		ResultProcessor.filterSearchResultWithFDR(targetOutFileName, decoyOutFileName,  0.2f, aaSet);

			
		/*
		PSMList<MSGappedDictionaryPSM> targetPSM = MSGappedDictionaryParser.parse(targetOutFileName, aaSet).getDistinctiveSpectralSet();
		//.getDistinctivePeptideSet()
		
		System.out.println("num spec" + targetPSM.size());
		int n=0;
		for(MSGappedDictionaryPSM psm : targetPSM)
			if(psm.isPeptideModified()) n++;
		System.out.println("num spec modified" + n);
		targetPSM = targetPSM.getDistinctivePeptideSet();
		
		System.out.println(targetPSM.size());
		n=0;
		for(MSGappedDictionaryPSM psm : targetPSM)
			if(psm.isPeptideModified()){
				n++;
				gapPep.add(psm.getPeptideStr());
			}
		
		for(String g : insPep){
			if(!gapPep.contains(g))
				System.out.println(g);
		}
		
		System.out.println(n);
		*/
	}
}
