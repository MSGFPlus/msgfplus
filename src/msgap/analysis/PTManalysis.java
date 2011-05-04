package msgap.analysis;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

import parser.MSGappedDictionaryPSM;
import parser.MSGappedDictionaryParser;
import parser.MgfSpectrumParser;
import parser.PSMList;
import msutil.AminoAcidSet;
import msutil.SpectraIterator;
import msutil.Spectrum;

public class PTManalysis {
	public static void main(String[] argv) throws IOException{
		String targetOutFileName = "/home/kwj/workspace/outputs/PTM/withPTM/out.txt";
		String aaFileName = "/home/kwj/workspace/inputs/AA_shewanella.txt";
		String spectrumFIleName = "/home/kwj/workspace/inputs/shewLengthAll.mgf";
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(aaFileName);
		
		PSMList<MSGappedDictionaryPSM> targetPSM = MSGappedDictionaryParser.parse(targetOutFileName, aaSet).getDistinctiveSpectralSet().getDistinctivePeptideSet();
		
		int numModified = 0;
		
		HashMap<Integer, String> peptideMap = new HashMap<Integer, String>();
		
		for(MSGappedDictionaryPSM psm : targetPSM){
			if(psm.isPeptideModified()) numModified++;
			peptideMap.put(psm.getScanNum(), psm.getPeptideStr());
		}
		
		System.out.println(numModified);
		
		
		Iterator<Spectrum> iterator = new SpectraIterator(spectrumFIleName, new MgfSpectrumParser());
		
		int numSpec = 0, numCorrectSpec = 0;
		
		while(iterator.hasNext())
		{
			Spectrum spec = iterator.next();
			
			if(!peptideMap.containsKey(spec.getScanNum())) continue;
			numSpec++;
			if(peptideMap.get(spec.getScanNum()).equals(spec.getAnnotationStr()))numCorrectSpec++;
			
		}
		
		System.out.println("Identified spec: " + numSpec + "\tCorrect : " + numCorrectSpec);
	}
}
