package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import parser.BufferedLineReader;
import parser.InsPecTPSM;
import parser.InsPecTParser;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraMap;

import msgf.Histogram;
import msutil.AminoAcidSet;
import msutil.Peptide;
import msutil.SpectraMap;
import msutil.Spectrum;
import msutil.SpectrumAccessorBySpecIndex;

public class PhosAnalysis {
	public static void main(String argv[]) throws Exception
	{
//		makeAnnotatedMgf();
		compareMSGFDBAndInsPecT();
	}

	public static void compareMSGFDBAndInsPecT() throws Exception
	{
		String inspectResultFile = "/home/sangtaekim/Test/JunePhospho/Inspect/finalResult/dee75f91bf354b62aed3fc3436645826";
		String msgfdbResultFile = "/home/sangtaekim/Test/JunePhospho/MSGFDB/finalResult/5ab63ac652d840f58ab9e598cf796975";
		countIDByNumPhospho(inspectResultFile);
		countIDByNumPhospho(msgfdbResultFile);
	}

	public static void countIDByNumPhospho(String fileName) throws Exception
	{
		System.out.println(fileName);
		int idCol = -1;
		
		BufferedLineReader in = new BufferedLineReader(fileName);
		String header = in.readLine();
		String[] field = header.split("\t");
		
		Histogram<Integer> hist = new Histogram<Integer>();
		for(int i=0; i<field.length; i++)
		{
			if(field[i].equals("Peptide") || field[i].equals("Annotation"))
				idCol = i;
		}
		
		String s;
		int numID = 0;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length <= idCol)
				continue;
			String annotation = token[idCol];
			if(!annotation.contains("phos") && !annotation.contains("79.9"))
				continue;
			
			int numPhospho = 0;
			for(int i=0; i<annotation.length(); i++)
			{
				String substr = annotation.substring(i);
				if(substr.startsWith("phos") || substr.startsWith("+79.966"))
					numPhospho++;
			}
			hist.add(numPhospho);
			if(numPhospho > 0)
				numID++;
		}
		System.out.println("NumID: " + numID);
		hist.printSorted();
	}	
	
	public static void countIDByMSLevel(String fileName) throws Exception
	{
		System.out.println(fileName);
		File specDir = new File("/home/sangtaekim/Test/JunePhospho/Inspect/spectrum");
		int specFileCol = -1;
		int specIndexCol = -1;
		int idCol = -1;
		
		BufferedLineReader in = new BufferedLineReader(fileName);
		String header = in.readLine();
		String[] field = header.split("\t");
		
		Histogram<String> hist = new Histogram<String>();
		for(int i=0; i<field.length; i++)
		{
			if(field[i].equals("#SpecFile") || field[i].equals("#SpectrumFile"))
				specFileCol = i;
			else if(field[i].equals("SpecIndex") || field[i].equals("Scan#"))
				specIndexCol = i;
			else if(field[i].equals("Peptide") || field[i].equals("Annotation"))
				idCol = i;
		}
		
		HashMap<String,SpectrumAccessorBySpecIndex> specAccessorMap = new HashMap<String,SpectrumAccessorBySpecIndex>(); 
		String s;
		int numID = 0;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length <= specFileCol || token.length <= specIndexCol || token.length <= idCol)
				continue;
			String specFileName = token[specFileCol];
			specFileName = new File(specFileName).getName();
			int specIndex = Integer.parseInt(token[specIndexCol]);
			String annotation = token[idCol];
			if(!annotation.contains("phos") && !annotation.contains("79.9"))
				continue;
			
			SpectrumAccessorBySpecIndex specMap = specAccessorMap.get(specFileName);
			if(specMap == null)
			{
				String ext = specFileName.substring(specFileName.lastIndexOf('.'));
				if(ext.equalsIgnoreCase(".mzXML"))
					specMap = new MzXMLSpectraMap(specDir.getPath()+File.separator+specFileName);
				else if(ext.equalsIgnoreCase(".mgf"))
					specMap = new SpectraMap(specDir.getPath()+File.separator+specFileName, new MgfSpectrumParser());
				else
				{
					System.out.println("Unrecognized spectrum format: " + specFileName);
					System.exit(-1);
				}
				specAccessorMap.put(specFileName, specMap);
			}
			Spectrum spec = specMap.getSpectrumBySpecIndex(specIndex);
			hist.add(specFileName);
			numID++;
		}
		System.out.println("NumID: " + numID);
		hist.printSorted();
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
							Spectrum spec = map.getSpectrumBySpecIndex(psm.getScanNum());
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
