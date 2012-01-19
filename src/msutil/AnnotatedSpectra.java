package msutil;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import parser.BufferedLineReader;

public class AnnotatedSpectra {
	private File[] resultFiles;
	private File specDir;
	private AminoAcidSet aaSet;
	private float fdrThreshold = 0.01f;
	
	private SpectraContainer annotatedSpectra;
	
	public AnnotatedSpectra(File[] resultFiles, File specDir, AminoAcidSet aaSet)
	{
		this.resultFiles = resultFiles;
		this.specDir = specDir;
		this.aaSet = aaSet;
	}
	
	public AnnotatedSpectra fdrThreshold(float fdrThreshold)
	{
		this.fdrThreshold = fdrThreshold;
		return this;
	}
	
	public SpectraContainer getAnnotatedSpecContainer()	{ return annotatedSpectra; }
	
	public String parse()
	{
		annotatedSpectra = new SpectraContainer();
		
		for(File resultFile : resultFiles)
		{
			String errMsg = parseFile(resultFile);
			if(errMsg != null)
				return "Error while parsing " + resultFile.getName() + ": " + errMsg;
		}
		return null;
	}
	
	public String parseFile(File resultFile)
	{
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(resultFile.getPath());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		String s = in.readLine();
		
		if(!s.startsWith("#"))
		{
			return "Not a valid tsv result file";
		}
		
		int specIndexCol = -1;
		int specFileCol = -1;
		int pepCol = -1;
		int fdrCol = -1;
		int chargeCol = -1;
		
		String[] label = s.split("\t");
		for(int i=0; i<label.length; i++)
		{
			if(label[i].equalsIgnoreCase("#SpecFile"))
				specFileCol = i;
			else if(label[i].equalsIgnoreCase("SpecIndex"))
				specIndexCol = i;
			else if(label[i].equalsIgnoreCase("Peptide"))
				pepCol = i;
			else if(label[i].equalsIgnoreCase("FDR") || label[i].equalsIgnoreCase("EFDR"))
				fdrCol = i;
			else if(label[i].equalsIgnoreCase("Charge"))
				chargeCol = i;
		}
		if(specIndexCol < 0 || specFileCol < 0 || pepCol < 0 || fdrCol < 0)
			return "Not a valid tsv result file";
		
		ArrayList<String> resultList = new ArrayList<String>();
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length <= specIndexCol || token.length <= specFileCol || token.length <= pepCol || token.length <= fdrCol)
				continue;
			
			float fdr = Float.parseFloat(token[fdrCol]);
			
			if(fdr <= fdrThreshold)
			{
				resultList.add(s);
			}
		}
		
		Iterator<String> itr = resultList.iterator();

		HashMap<String,SpectrumAccessorBySpecIndex> specAccessorMap = new HashMap<String,SpectrumAccessorBySpecIndex>(); 
		while(itr.hasNext())
		{
			String str = itr.next();
			String[] token = str.split("\t");
			
			String pep = token[pepCol];
			if(pep.matches(".\\..+\\.."))
				pep = pep.substring(pep.indexOf('.')+1, pep.lastIndexOf('.'));
			
			String specFileName = token[specFileCol];
			specFileName = new File(specFileName).getName();
			
			int charge = Integer.parseInt(token[chargeCol]);
			
			SpectrumAccessorBySpecIndex specMap = specAccessorMap.get(specFileName);
			if(specMap == null)
			{
				File specFile = new File(specDir.getPath()+File.separator+specFileName);
				specMap = SpecFileFormat.getSpecMap(specFile);
				if(specMap == null)
					return "Unrecognized spectrum format";
				specAccessorMap.put(specFileName, specMap);
			}
			
			int specIndex = Integer.parseInt(token[specIndexCol]);
			Spectrum spec = specMap.getSpectrumBySpecIndex(specIndex);
			
			if(spec == null)
				return specFileName+":"+specIndex+" is not available!";
			else
			{
				Peptide peptide = new Peptide(pep, aaSet);
				spec.setCharge(charge);
				
				if(Math.abs(spec.getPeptideMass()-peptide.getMass()) < 5)
				{
					spec.setAnnotation(peptide);
					annotatedSpectra.add(spec);
				}
				else
				{
					return "parent mass doesn't match " + specFileName+":"+specIndex + " " + peptide.toString() + " " + spec.getPeptideMass() + " != " + peptide.getMass();
				}
			}
		}
		
		try {
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}	
}
