package edu.ucsd.msjava.parser;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.util.Collections;

import edu.ucsd.msjava.msutil.AminoAcidSet;

public class MSGappedDictionaryParser {
	
	public static PSMList<MSGappedDictionaryPSM> parse(String fileName){
		return parse(fileName, AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
	}
	
	public static PSMList<MSGappedDictionaryPSM> parse(String fileName, AminoAcidSet aaSet)
	{
		PSMList<MSGappedDictionaryPSM> psmList = new PSMList<MSGappedDictionaryPSM>();
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#") || s.length() == 0)
				continue;
			String[] token = s.split("\t");
			if(token.length < 13)
				continue;
			String specFileName = token[0];
			int scanNum = Integer.parseInt(token[1]);
//			String method = token[2];
			float precursorMz = Float.parseFloat(token[3]);
			int charge = Integer.parseInt(token[4]);
			String peptideStr = token[5].substring(token[5].indexOf('.')+1, token[5].lastIndexOf('.'));
			String protein = token[6];
			
//			int msgfScore = Integer.parseInt(token[7]);
			int peptideScore = Integer.parseInt(token[7]);
			float specProb = Float.parseFloat(token[8]);
			float massDiff = Float.parseFloat(token[12]);
			
			MSGappedDictionaryPSM psm = new MSGappedDictionaryPSM();
			psm.aaSet(aaSet).peptide(peptideStr).specFileName(specFileName).scanNum(scanNum).precursorMz(precursorMz).charge(charge)
			.protein(protein).rawScore(peptideScore).probScore(specProb);
			psmList.add(psm);
			
			psm.setParentMassError(massDiff);
		}
		return psmList;
	}
	
	public static PSMList<MSGappedDictionaryPSM> getMergedResults(String dirName, final String suffix)
	{
		File dir = new File(dirName);
		if(!dir.isDirectory())
			return null;
		class SuffixFileFilter implements FileFilter {
			public boolean accept(final File pathname) {
				if(pathname.getName().endsWith(suffix))
					return true;
				else
					return false;
			}
		}
		PSMList<MSGappedDictionaryPSM> psmList = new PSMList<MSGappedDictionaryPSM>();
		for(File f : dir.listFiles(new SuffixFileFilter()))
		{
			psmList.addAll(parse(f.getPath()));
		}
		Collections.sort(psmList, new PSM.PSMSpecNumComparator());
		return psmList;
	}
}
