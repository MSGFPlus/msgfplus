package trex.analysis;

import parser.MSGappedDictionaryPSM;
import parser.MSGappedDictionaryParser;
import parser.PSM;
import parser.PSMList;
import msgap.ResultProcessor;
import msgap.results.DeprecatedSpectrumMatches;
import msgf.Tolerance;
import msutil.AminoAcidSet;

public class newTrex {
	public static void main(String[] argv){
		String outFileName = "/home/kwj/workspace/outputs/newTrex/analysis/match.txt";
		//String outFileName = "/home/kwj/workspace/outputs/newTrex/analysis/match.decoy_out";
		String aaFileName = "/home/kwj/workspace/inputs/newTrex/AASetWithHydroProline.txt";
		
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(aaFileName);
		ResultProcessor.rewriteResultWithThreshold(outFileName, outFileName + "_filtered", 1e-9f);
		
		PSMList<MSGappedDictionaryPSM> filteredTargetPSM = MSGappedDictionaryParser.parse(outFileName + "_filtered", aaSet).getDistinctiveSpectralSet();
		
		System.out.println("# spectra: " + filteredTargetPSM.size());
		
		filteredTargetPSM = filteredTargetPSM.getDistinctivePeptideSet();
		for(int i=0; i<filteredTargetPSM.size(); i++){
			MSGappedDictionaryPSM psm = filteredTargetPSM.get(i);
			String protein = psm.getProtein().toLowerCase();
			if(protein.contains("keratin") || 
			   protein.contains("trypsin")||
					
			(!protein.contains("collagen") && psm.getPeptideStr().contains("p"))){
				filteredTargetPSM.remove(i--);
			}
		}
		
		Tolerance t = new Tolerance(50, true);
		
		for(MSGappedDictionaryPSM psm : filteredTargetPSM){
			String fileName = psm.getSpecFileName();
			fileName = fileName.substring(fileName.lastIndexOf('/')+1, fileName.length());
			
			if(t.getToleranceAsDa(psm.getPeptide().getParentMass()) < psm.getParentMassError()) continue;
			
		//	System.out.println(fileName + "\t" + psm.getScanNum() + "\t" + psm.getProtein() + "\t" + psm.getPeptideStr() + "\t" + psm.getProbScore() + "\t" + psm.getParentMassError());
		}
	}
}
