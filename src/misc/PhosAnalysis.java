package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import parser.InsPecTPSM;
import parser.InsPecTParser;
import parser.MzXMLSpectraMap;

import msutil.AminoAcidSet;
import msutil.Peptide;
import msutil.Spectrum;

public class PhosAnalysis {
	public static void main(String argv[]) throws Exception
	{
		makeAnnotatedMgf();
	}
	
	public static void makeAnnotatedMgf() throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/Phospho");
		File resultDir = new File(dir.getPath()+"/inspect");
		File specRootDir = new File(dir.getPath()+"/spectra");
		File annotatedSpectra = new File(dir.getPath()+"/annotatedPhospho.mgf");
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(annotatedSpectra)));
		
		AminoAcidSet baseAASet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarboxymethylatedCys();
		
		int numSpec = 0;
		for(File f : resultDir.listFiles())
		{
			if(f.getName().endsWith(".txt"))
			{
				System.out.println(f.getName());
				InsPecTParser parser = new InsPecTParser(baseAASet);
				parser.parse(f.getPath());
				ArrayList<InsPecTPSM> psmList = parser.getPSMList();
				if(psmList == null || psmList.size() == 0)
					continue;
				// acquiring the spectrum
				InsPecTPSM firstPSM = psmList.get(0);
				String specFileName = firstPSM.getSpecFileName();
				String specSubDirName = specFileName.substring(0, specFileName.lastIndexOf('-'));
				File specDir = new File(specRootDir.getPath()+"/"+specSubDirName);
				assert(specDir.exists() && specDir.isDirectory()): specDir.getName() + " is missing or not a directory!";
				File specFile = new File(specDir.getPath()+"/"+specFileName);
				assert(specFile.exists()): specFile.getName() + " does not exist!";
				MzXMLSpectraMap map = new MzXMLSpectraMap(specFile.getPath());
				for(InsPecTPSM psm : psmList)
				{
					if(psm.getProbScore() > 0.01f)
						continue;
					Peptide pep = psm.getPeptide();
					if(pep.isModified())
					{
						String pepString = pep.toString();
						if(pepString.contains("s") || pepString.contains("t"))
						{
							Spectrum spec = map.getSpectrumByScanNum(psm.getScanNum());
							spec.setCharge(psm.getCharge());
							float expPM = spec.getParentMass();
							float theoPM = pep.getParentMass();
							float massDiff = expPM - theoPM;
							if(Math.abs(massDiff) > 5f)
							{
								System.out.println("Error! mass mismatch: " + 
										psm.getSpecFileName()+"\t"+psm.getScanNum()+"\t"+pep+"\t"+expPM+"!="+theoPM);
								spec.setCharge(psm.getCharge());
								continue;
							}
							numSpec++;
							spec.setAnnotation(pep);
							spec.outputMgf(out);
						}
					}
				}
			}
		}
		out.close();
	}
}
