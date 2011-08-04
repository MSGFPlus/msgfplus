package parser;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Hashtable;

import msutil.SpectraMap;
import msutil.Spectrum;
import msutil.SpectrumAccessorBySpecIndex;

public class PSMList<T extends PSM> extends ArrayList<T> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	public PSMList()	{ super(); }
	public PSMList(Collection<T> c) { super(c); }
	
	/**
	 * Select the best PSM among multiple PSMs belong to the same spectrum. 
	 * @return list of PSMs, one PSM per spectrum.
	 */
	public PSMList<T> getDistinctiveSpectralSet()
	{
		Hashtable<Integer, T> table = new Hashtable<Integer, T>();
		for(T psm : this)
		{
			T bestPSM = table.get(psm.getScanNum());
			if(bestPSM == null || psm.getProbScore() < bestPSM.getProbScore())
				table.put(psm.getScanNum(), psm);
		}
		PSMList<T> distinctSpecList = new PSMList<T>(table.values());
		Collections.sort(distinctSpecList, new PSM.PSMSpecNumComparator());
		return distinctSpecList;
	}
	
	/**
	 * Select the best PSM among multiple PSMs belong to the same peptide. 
	 * @return list of PSMs, one PSM per peptide.
	 */
	public PSMList<T> getDistinctivePeptideSet()
	{
		Hashtable<String, T> table = new Hashtable<String, T>();
		for(T psm : this)
		{
			T bestPSM = table.get(psm.getPeptideStr());
			if(bestPSM == null || psm.getProbScore() < bestPSM.getProbScore())
				table.put(psm.getPeptideStr(), psm);
		}
		PSMList<T> distinctPepList = new PSMList<T>(table.values());
		Collections.sort(distinctPepList, new PSM.PSMSpecNumComparator());
		return distinctPepList;
	}
	
	public static<T extends PSM> PSMList<T> selectUsingFDR(PSMList<T> targetPSMList, PSMList<T> decoyPSMList, float fdrThreshold)
	{
		PSMList<T> selected = new PSMList<T>();
		
		Collections.sort(targetPSMList, new PSM.PSMProbScoreComparator());
		Collections.sort(decoyPSMList, new PSM.PSMProbScoreComparator());
		
		int indexTarget = -1;
		int indexDecoy = 0;
		float fdr = 0;
		for(T targetPSM : targetPSMList)
		{
			indexTarget++;
			if(indexDecoy >= decoyPSMList.size())
			{
				if(fdr <= fdrThreshold)
					selected.add(targetPSM);
			}
			else if(targetPSM.getProbScore() < decoyPSMList.get(indexDecoy).getProbScore())
				selected.add(targetPSM);
			else
			{
				while(indexDecoy < decoyPSMList.size() && targetPSM.getProbScore() >= decoyPSMList.get(indexDecoy).getProbScore())
					indexDecoy++;
				fdr = indexDecoy / (float)(indexTarget+1);
				if(fdr <= fdrThreshold)
					selected.add(targetPSM);
				else
					break;
			}
		}
		
		return selected;
	}
	
	public boolean outputMgf(File specDir, File outputFile, String scoreName, float threshold, boolean isGreaterBetter)
	{
		PrintStream out = null;
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		int numSpecs = 0;
		Hashtable<String, SpectrumAccessorBySpecIndex> fileList = new Hashtable<String, SpectrumAccessorBySpecIndex>();
		for(PSM psm : this)
		{
//			if(psm.getPeptide().isModified())	// modified peptides are ignored
//				continue;
			if(psm.getProtein().contains("DECOY"))
				continue;
			Float score = psm.getScore(scoreName);
			if(score == null)
			{
				System.err.println(scoreName + " is not supported.");
				return false;
			}
			if(isGreaterBetter && score < threshold ||
					!isGreaterBetter && score > threshold)
				continue;
				
			SpectrumAccessorBySpecIndex specAccessor = null;
			String fileName = psm.getSpecFileName();
			if((specAccessor = fileList.get(fileName)) == null)
			{
				File specFile = new File(specDir.getPath()+File.separator+psm.getSpecFileName());
				if(!specFile.exists())
				{
					System.err.println(psm.getSpecFileName() + " doesn't exist in " + specDir.getPath());
					return false;
				}
				String ext = fileName.substring(fileName.lastIndexOf('.'));
				if(ext.equalsIgnoreCase(".mzxml") || ext.equalsIgnoreCase(".xml"))	// mzXML
				{
					specAccessor = new MzXMLSpectraMap(specFile.getPath());
				}
				else if(ext.equalsIgnoreCase(".mgf"))
				{
					specAccessor = new SpectraMap(specFile.getPath(), new MgfSpectrumParser());
				}
				else if(ext.equalsIgnoreCase(".pkl"))
				{
					specAccessor = new SpectraMap(specFile.getPath(), new PklSpectrumParser());
				}
				else if(ext.equalsIgnoreCase(".ms2"))
				{
					specAccessor = new SpectraMap(specFile.getPath(), new MS2SpectrumParser());
				}
				else
				{
					System.err.println("Unrecognizable format: " + psm.getSpecFileName());
					return false;
				}
				fileList.put(fileName, specAccessor);
			}			
			
			Spectrum spec = specAccessor.getSpectrumBySpecIndex(psm.getScanNum());
			if(spec == null)
			{
				System.err.println("Spectrum doesn't exist: " + psm.getSpecFileName() + ":" + psm.getScanNum());
				return false;
			}
			spec.setTitle(psm.getSpecFileName()+":"+psm.getScanNum()+" "+scoreName+"="+score);
			spec.setAnnotation(psm.getPeptide());
			spec.outputMgf(out);
			numSpecs++;
		}
		out.close();
		System.out.println(numSpecs + " spectra are stored.");
		return true;
	}
	
	// added by kyowon
	public static<T extends PSM> PSM getBestPSM(PSMList<T> targetPSMList){
		if(targetPSMList.isEmpty()) return null;
		Collections.sort(targetPSMList, new PSM.PSMProbScoreComparator());
		return targetPSMList.get(0);
	}
	
	// added by kyowon
	public static<T extends PSM> PSM getWorstPSM(PSMList<T> targetPSMList){
		if(targetPSMList.isEmpty()) return null;
		Collections.sort(targetPSMList, new PSM.PSMProbScoreComparator());
		return targetPSMList.get(targetPSMList.size() - 1);
	}
	
}
