package msgap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import msutil.AminoAcidSet;
import parser.BufferedLineReader;
import parser.MSGappedDictionaryPSM;
import parser.MSGappedDictionaryParser;
import parser.PSM;
import parser.PSMList;

public class ResultProcessor {
	static public void rewriteResultWithThreshold(String fileName, String outFileName, float specProbThreshold){
		BufferedLineReader in;
		
		//new File(fileName).renameTo(new File(fileNameBeforeRewrite));
		
		try {
			in = new BufferedLineReader(fileName);
			BufferedWriter out = new BufferedWriter(new FileWriter(outFileName));
			
			String s;
			
			while((s=in.readLine()) != null){
				if(s.startsWith("#") || s.length() == 0){
					out.write(s+"\n");
					continue;
				}
				String[] token = s.split("\t");
				if(Float.parseFloat(token[9]) > specProbThreshold) continue;
				out.write(s+"\n");
			}
			
			in.close();
			out.close();
			
			
		} catch (IOException e) {
			 System.err.println("IOException caught when rewritng the result file " + fileName);
			 e.printStackTrace();
			 System.exit(-1);
		}
	}
	
	
	// per spec file name
	static public void filterSearchResultWithFDR(String targetOutFileName, String decoyOutFileName,  float FDR, AminoAcidSet aaSet){
		PSMList<MSGappedDictionaryPSM> targetPSM = MSGappedDictionaryParser.parse(targetOutFileName, aaSet).getDistinctiveSpectralSet().getDistinctivePeptideSet();
		PSMList<MSGappedDictionaryPSM> decoyPSM = MSGappedDictionaryParser.parse(decoyOutFileName, aaSet).getDistinctiveSpectralSet().getDistinctivePeptideSet();
		PSMList<MSGappedDictionaryPSM> significantPSMList = PSMList.selectUsingFDR(targetPSM, decoyPSM, FDR);

		PSMList<MSGappedDictionaryPSM> significantPSMPeptideList = significantPSMList.getDistinctivePeptideSet();
		
		MSGappedDictionaryPSM worst = (MSGappedDictionaryPSM)PSMList.getWorstPSM(significantPSMPeptideList);
		
		float threshold = 1;
		
		if(worst != null)
			threshold = worst.getProbScore();
			
		new File(targetOutFileName).renameTo(new File(targetOutFileName + "_before_filter"));
		new File(decoyOutFileName).renameTo(new File(decoyOutFileName + "_before_filter"));
		
		rewriteResultWithThreshold(targetOutFileName + "_before_filter", targetOutFileName, threshold);
		rewriteResultWithThreshold(decoyOutFileName + "_before_filter", decoyOutFileName, threshold);
		
		System.out.println("\nThreshold spec prob: " + threshold);

	}
	
	static public float getFDR(String targetOutFileName, String decoyOutFileName, AminoAcidSet aaSet){
		PSMList<MSGappedDictionaryPSM> targetPSM = MSGappedDictionaryParser.parse(targetOutFileName, aaSet).getDistinctiveSpectralSet().getDistinctivePeptideSet();
		PSMList<MSGappedDictionaryPSM> decoyPSM = MSGappedDictionaryParser.parse(decoyOutFileName, aaSet).getDistinctiveSpectralSet().getDistinctivePeptideSet();
		
		return (float) decoyPSM.size() / (targetPSM.size() + decoyPSM.size());
	}
}
