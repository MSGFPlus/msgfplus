package parser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Peptide;

public class InsPecTParser {

	// InsPecT labels
	public static final String SPEC_FILE = "#SpectrumFile";
	public static final String SCAN_NUM = "Scan#";
	public static final String ANNOTATION = "Annotation";
	public static final String PROTEIN = "Protein";
	public static final String CHARGE = "Charge";
	public static final String MQ_SCORE = "MQScore";
	public static final String FDR = "p-Value";
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
		int scanNumColumn = -1;
		int annotationColumn = -1;
		int proteinColumn = -1;
		int chargeColumn = -1;
		int mqScoreColumn = -1;
		int fdrColumn = -1;
		@SuppressWarnings("unused")
		int fScoreColumn = -1;
		int specFilePosColumn = -1;
		String[] label = labelRow.split("\t");
		for(int i=0; i<label.length; i++)
		{
			if(label[i].equalsIgnoreCase(SPEC_FILE) || label[i].equalsIgnoreCase("#SpecFile"))
				specFileColumn = i;
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
		}
		
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
			String[] scanNumToken = token[scanNumColumn].split("-");
			int scanNum = Integer.parseInt(scanNumToken[0]);
			ArrayList<Integer> scanNumList = new ArrayList<Integer>();
			for(int i=0; i<scanNumToken.length; i++)
				scanNumList.add(Integer.parseInt(scanNumToken[i]));
			
			String annotationStr = token[annotationColumn];
			String proteinStr;
			if(proteinColumn >= 0)
				proteinStr = token[proteinColumn];
			else
				proteinStr = "";
			int charge = Integer.parseInt(token[chargeColumn]);
			float mqScore = Float.NaN;
			float fdr = Float.NaN;
			long specFilePos = -1;
			if(mqScoreColumn >= 0)
					mqScore = Float.parseFloat(token[mqScoreColumn]);
			if(fdrColumn >= 0)
				fdr = Float.parseFloat(token[fdrColumn]);
			if(specFilePosColumn >= 0)
				specFilePos = Long.parseLong(token[specFilePosColumn]);
//			float fScore = Float.parseFloat(token[fScoreColumn]);

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
			Peptide peptide = baseAASet.getPeptide(pepStr);
				
			
			InsPecTPSM psm = new InsPecTPSM();
			psm.scanNum(scanNum).peptide(peptide).protein(proteinStr).charge(charge).probScore(fdr).rawScore(mqScore);
			if(s.endsWith("\t"))
				s = s.substring(0, s.length()-1);
			psm.setInsPecTString(s);
			psm.setPrecedingAA(preAA);
			psm.setSucceedingAA(nextAA);
			psm.setScanNumList(scanNumList);
			
//			AminoAcidSet unmodAASet = baseAASet.getUnmodifiedAminoAcidSet();
			//TODO: revise this
			AminoAcidSet unmodAASet = baseAASet;
			HashSet<AminoAcid> modifiedAASet = new HashSet<AminoAcid>();
			if(peptide != null && peptide.isModified())	// TODO: to be modified
			{
				for(AminoAcid aa : peptide)
				{
					if(aa.isModified())	// modified residue
						modifiedAASet.add(aa);
				}
				//TODO: revise this
//				AminoAcidSet newAASet = unmodAASet.getAminoAcidSet(modifiedAASet);
				AminoAcidSet newAASet = baseAASet;
				psm.setAASet(newAASet);
			}
			else 
				psm.setAASet(unmodAASet);
			
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
