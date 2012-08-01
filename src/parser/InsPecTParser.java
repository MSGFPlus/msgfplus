package parser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Modification;
import msutil.Peptide;

public class InsPecTParser {

	// InsPecT labels
	public static final String SPEC_FILE = "#SpectrumFile";
	public static final String SCAN_NUM = "Scan#";
	public static final String SPEC_INDEX = "SpecIndex";
	public static final String ANNOTATION = "Annotation";
	public static final String PROTEIN = "Protein";
	public static final String CHARGE = "Charge";
	public static final String MQ_SCORE = "MQScore";
	public static final String FDR = "FDR";
	public static final String SPEC_PROB = "SpecProb";
	public static final String F_SCORE = "F-Score";
	public static final String SPEC_FILE_POS = "SpecFilePos";
	
	private AminoAcidSet baseAASet;
	private String header;
	private PSMList<InsPecTPSM> psmList;
	
	public InsPecTParser(AminoAcidSet baseAASet)
	{
		this.baseAASet = baseAASet;
		header = null;
		psmList = null;
	}
	
	public String getHeader()	{ return header; }
	public PSMList<InsPecTPSM> getPSMList()	{ return psmList; }
	
	public PSMList<InsPecTPSM> getPSMList(String scoreName, float threshold, boolean isBiggerBetter)
	{
		PSMList<InsPecTPSM> filteredList = new PSMList<InsPecTPSM>();
		
		PSMList<InsPecTPSM> list = psmList.getDistinctiveSpectralSet();
		for(InsPecTPSM psm : list)
		{
			float score = psm.getScore(scoreName);
			if(isBiggerBetter)
			{
				if(score < threshold)
					continue;
			}
			else
			{
				if(score > threshold)
					continue;
			}
			filteredList.add(psm);
		}
		return filteredList;
	}
	
	public void parse(String fileName)
	{
		PSMList<InsPecTPSM> psmList = new PSMList<InsPecTPSM>();
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		String labelRow = in.readLine();
		
		if(labelRow == null || !labelRow.startsWith("#"))
		{
			return;	// illegal file format
		}
		
		int specFileColumn = -1;
		int specIndexColumn = -1;
		int scanNumColumn = -1;
		int annotationColumn = -1;
		int proteinColumn = -1;
		int chargeColumn = -1;
		int mqScoreColumn = -1;
		int fdrColumn = -1;
		int fScoreColumn = -1;
		int specProbColumn = -1;
		int specFilePosColumn = -1;
		String[] label = labelRow.split("\t");
		for(int i=0; i<label.length; i++)
		{
			if(label[i].equalsIgnoreCase(SPEC_FILE) || label[i].equalsIgnoreCase("#SpecFile"))
				specFileColumn = i;
			else if(label[i].equalsIgnoreCase(SPEC_INDEX))
				specIndexColumn = i;
			else if(label[i].equalsIgnoreCase(SCAN_NUM))
				scanNumColumn = i;
			else if(label[i].equalsIgnoreCase(ANNOTATION) || label[i].equalsIgnoreCase("Peptide"))
				annotationColumn = i;
			else if(label[i].equalsIgnoreCase(PROTEIN))
				proteinColumn = i;
			else if(label[i].equalsIgnoreCase(CHARGE))
				chargeColumn = i;
			else if(label[i].equalsIgnoreCase(MQ_SCORE))
				mqScoreColumn = i;
			else if(label[i].equalsIgnoreCase(FDR))
				fdrColumn = i;
			else if(label[i].equalsIgnoreCase(SPEC_FILE_POS))
				specFilePosColumn = i;
			else if(label[i].equalsIgnoreCase(F_SCORE))
				fScoreColumn = i;
			else if(label[i].equalsIgnoreCase(SPEC_PROB))
				specProbColumn = i;
		}
		
		// if there is specIndex column, use it instead of scanNum
		if(specIndexColumn >= 0)
			scanNumColumn = specIndexColumn;
		
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length < specFileColumn || token.length < scanNumColumn || token.length < annotationColumn)
				continue;
			String specFileName = token[specFileColumn].trim();

			// parse scan number(s)
			ArrayList<Integer> scanNumList = new ArrayList<Integer>();
			int scanNum;
			if(token[scanNumColumn].equals("-1"))
				scanNum = -1;
			else
			{
				String[] scanNumToken = token[scanNumColumn].split("-");
				scanNum = Integer.parseInt(scanNumToken[0]);
				for(int i=0; i<scanNumToken.length; i++)
					scanNumList.add(Integer.parseInt(scanNumToken[i]));
			}
			
			String annotationStr = token[annotationColumn];
			String proteinStr;
			if(proteinColumn >= 0)
				proteinStr = token[proteinColumn];
			else
				proteinStr = "";
			int charge = Integer.parseInt(token[chargeColumn]);
			float mqScore = Float.NaN;
			float fScore = Float.NaN;
			float fdr = Float.NaN;
			float specProb = Float.NaN;
			int specIndex = -1;
			
			long specFilePos = -1;
			if(mqScoreColumn >= 0)
				mqScore = Float.parseFloat(token[mqScoreColumn]);
			if(fScoreColumn >= 0)
				fScore = Float.parseFloat(token[fScoreColumn]);
			if(fdrColumn >= 0)
				fdr = Float.parseFloat(token[fdrColumn]);
			if(specFilePosColumn >= 0)
				specFilePos = Long.parseLong(token[specFilePosColumn]);
			if(specProbColumn >= 0)
			{
				if(token[specProbColumn].startsWith("N/A"))
					continue;
				specProb = Float.parseFloat(token[specProbColumn]);
			}
			if(specIndexColumn >= 0)
				specIndex = Integer.parseInt(token[specIndexColumn]);
			// parse specIndex

			// process specFileName
			if(specFileName.contains("/"))
				specFileName = specFileName.substring(specFileName.lastIndexOf('/')+1);
			else if(specFileName.contains("\\"))
				specFileName = specFileName.substring(specFileName.lastIndexOf('\\')+1);
			else if(specFileName.contains(File.separator))
				specFileName = specFileName.substring(specFileName.lastIndexOf(File.separatorChar)+1);
			
			// process annotation
			AminoAcid preAA = null;
			AminoAcid nextAA = null;
			int firstDotPos = annotationStr.indexOf('.');
			int lastDotPos = annotationStr.lastIndexOf('.');
			if(firstDotPos < lastDotPos)	// there are two dots in annotationStr
			{
				String preAAStr = annotationStr.substring(0, firstDotPos);
				assert(preAAStr.length()<=1);
				if(preAAStr.length() == 0)
					preAA = null;
				else
					preAA = baseAASet.getAminoAcid(preAAStr.charAt(0));
				String sucAAStr = annotationStr.substring(lastDotPos+1);
				assert(sucAAStr.length()<=1);
				if(sucAAStr.length() == 0)
					nextAA = null;
				else
					nextAA = baseAASet.getAminoAcid(sucAAStr.charAt(0));
			}
			else
			{
				firstDotPos = -1;
				lastDotPos = annotationStr.length();
			}
			String pepStr = annotationStr.substring(firstDotPos+1, lastDotPos);
			Peptide peptide = new Peptide(pepStr, baseAASet);
			if(peptide.isInvalid())
				peptide = null;
			
			InsPecTPSM psm = new InsPecTPSM();
			psm.specIndex(specIndex).scanNum(scanNum).peptide(peptide).protein(proteinStr).charge(charge).probScore(fdr).rawScore(mqScore);
			
			if(fScore != Float.NaN)
				psm.score(F_SCORE, fScore);
			if(mqScore != Float.NaN)
				psm.score(MQ_SCORE, mqScore);
			if(fdr != Float.NaN)
				psm.score(FDR, fdr);
			if(specProb != Float.NaN)
				psm.score(SPEC_PROB, specProb);
			
			if(s.endsWith("\t"))
				s = s.substring(0, s.length()-1);
			psm.setInsPecTString(s);
			psm.setPrecedingAA(preAA);
			psm.setSucceedingAA(nextAA);
			psm.setScanNumList(scanNumList);
			psm.specFileName(specFileName);
			psm.setSpecFilePos(specFilePos);
			psmList.add(psm);
		}
		
		try {
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		this.header = labelRow.trim();
		this.psmList = psmList;
	}
	
	public static void main(String argv[])
	{
		String fileName = System.getProperty("user.home")+"/Research/ToolDistribution/TestForNatalie/test.txt";
		InsPecTParser parser = new InsPecTParser(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
		parser.parse(fileName);
	}
}
