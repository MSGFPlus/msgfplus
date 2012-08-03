package edu.ucsd.msjava.parser;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.util.Collections;

import edu.ucsd.msjava.msutil.Peptide;

public class MSGFDBParser {
	public static PSMList<PSM> parse(String fileName)
	{
		PSMList<PSM> psmList = new PSMList<PSM>();
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#") || s.length() == 0 || !Character.isDigit(s.charAt(0)))
				continue;
			String[] token = s.split("\t");
			if(token.length != 10)
				continue;
			String specFileName = token[0];
			int scanNum = Integer.parseInt(token[1]);
//			String method = token[2];
			float precursorMz = Float.parseFloat(token[3]);
			int charge = Integer.parseInt(token[4]);
			String peptideStr = token[5];
			String protein = token[6];
//			int msgfScore = Integer.parseInt(token[7]);
			int peptideScore = Integer.parseInt(token[8]);
			float specProb = Float.parseFloat(token[9]);
			PSM psm = new PSM();
			psm.specFileName(specFileName).scanNum(scanNum).precursorMz(precursorMz).charge(charge).peptide(new Peptide(peptideStr))
			.protein(protein).rawScore(peptideScore).probScore(specProb);
			psmList.add(psm);
		}
		return psmList;
	}
	
	public static PSMList<PSM> getMergedResults(String dirName, final String suffix)
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
		PSMList<PSM> psmList = new PSMList<PSM>();
		for(File f : dir.listFiles(new SuffixFileFilter()))
		{
			psmList.addAll(parse(f.getPath()));
		}
		Collections.sort(psmList, new PSM.PSMSpecNumComparator());
		return psmList;
	}
}
