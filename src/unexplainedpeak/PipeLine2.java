package unexplainedpeak;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import parser.MgfSpectrumParser;

import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;

public class PipeLine2 {
	static int rankLimit = 1;
	static private HashMap<String, ArrayList<Spectrum>> specPairMap = new HashMap<String, ArrayList<Spectrum>>();
	
	static public void main(String[] argv) throws IOException{
		
		boolean writemgf = false;
		
		String specFileName = "/home/kwj/workspace/inputs/OffsetTest/merged_charge2_inconsist_Unexp.mgf";
		String mgfOutFolder = "/home/kwj/workspace/inputs/OffsetTest/";
		String outFileName = specFileName.substring(0,specFileName.lastIndexOf('.'))+ "_fragPeptides_"+ rankLimit+".m";
		if(writemgf){
			gatherPeptides(specFileName);
			outputTripletsAndPairs(mgfOutFolder);
		}
		ArrayList<float[]> diffs = new ArrayList<float[]> ();
		diffs.addAll(getPeakPrecursorMzDiffs(mgfOutFolder + File.separatorChar + "pairs" + File.separatorChar + rankLimit));
		diffs.addAll(getPeakPrecursorMzDiffs(mgfOutFolder + File.separatorChar + "duals" + File.separatorChar + rankLimit));
		diffs.addAll(getPeakPrecursorMzDiffs(mgfOutFolder + File.separatorChar + "triplets" + File.separatorChar + rankLimit));
	
		writePeakPrecursorMzDiffs(outFileName, diffs);
		
		
	}
	
	static void writePeakPrecursorMzDiffs(String file, ArrayList<float[]> diffs) throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(file));
		out.write("out = [\n");
		for(float[] diff: diffs){
			out.write(diff[0]+"\t"+diff[1]+"\n");
		}
		out.write("];");
		
		out.close();
	}
	
	static ArrayList<float[]> getPeakPrecursorMzDiffs(String mgfOutFolder) throws IOException{
		File folder = new File(mgfOutFolder);
		File[] specfiles = folder.listFiles();
		
		ArrayList<float[]> peakPrecursorMzDiffs = new ArrayList<float[]>();
		
		for(File specfile : specfiles){
			Iterator<Spectrum> iterator = new SpectraIterator(specfile.getAbsolutePath(), new MgfSpectrumParser());
			ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
			while(iterator.hasNext()){
				spectra.add(iterator.next());
			}
			
			for(int i=0; i<spectra.size(); i++){
				Spectrum s1 = spectra.get(i);
				for(int j=i+1; j<spectra.size(); j++){
					Spectrum s2 = spectra.get(j);
					float precursorMzDiff = getPepMassDiff(s1,s2);
					for(float peakMz:getPeakMzDifferences(s1,s2)){
						float[] v = {precursorMzDiff, peakMz};
						peakPrecursorMzDiffs.add(v);
					}
				}
			}
			
		}
		
		return peakPrecursorMzDiffs;
	}
	
	
	static ArrayList<Float> getPeakMzDifferences(Spectrum s1, Spectrum s2){
		ArrayList<Float> diffs = new ArrayList<Float>();
		
		for(Peak p1 : s1)
			for(Peak p2 : s2)
				diffs.add(p1.getMz() - p2.getMz());
		
		return diffs;
	}
	
	static float getPepMassDiff(Spectrum s1, Spectrum s2){
		return s1.getAnnotation().getMass() - s2.getAnnotation().getMass();
	}
	
	static void gatherPeptides(String filename) throws IOException{
		Iterator<Spectrum> iterator = new SpectraIterator(filename, new MgfSpectrumParser());
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			String pep = spec.getAnnotationStr();
			ArrayList<Spectrum> value;
			
			if(specPairMap.containsKey(pep))
				value = specPairMap.get(pep);
			else
				value = new ArrayList<Spectrum>();
			
			value.add(spec);
			
			specPairMap.put(pep, value);
		}
	}
	
	
	
	static void outputTripletsAndPairs(String outFolder) throws IOException{
		ArrayList<String> toRemove = new ArrayList<String>();
		for(String pep : specPairMap.keySet()){
			if(toRemove.contains(pep)) continue;
			ArrayList<String> candidate = new ArrayList<String>();
			candidate.add(pep);
			int l = pep.length();
			for(String pep2 : specPairMap.keySet()){
				if(pep2.length() != l-1) continue;
				if(toRemove.contains(pep2) || pep.equals(pep2)) continue;
				if(pep.substring(1).equals(pep2) || pep.substring(0, pep.length()-1).equals(pep2))
					candidate.add(pep2);				
			}
			if(candidate.size() > 1){
				boolean isPair = true;
				if(candidate.size() == 3) isPair = false;					
				//System.out.println(candidate);
				toRemove.addAll(candidate);
				String filename = outFolder + File.separatorChar + (isPair? "pairs" : "triplets") + File.separatorChar  + rankLimit + File.separatorChar + candidate.get(0)+".mgf";
				PrintStream ps = new PrintStream(filename);
				
				for(String annotation : candidate){
					for(Spectrum spec : specPairMap.get(annotation)){
						spec.outputMgf(ps);
					}
				}
				ps.flush();ps.close();
				
			}else if(specPairMap.get(pep).size() > 1){
				String filename = outFolder + File.separatorChar + "duals" + File.separatorChar  + rankLimit + File.separatorChar + candidate.get(0)+".mgf";
				PrintStream ps = new PrintStream(filename);
				toRemove.addAll(candidate);
				for(String annotation : candidate){
					for(Spectrum spec : specPairMap.get(annotation)){
						spec.outputMgf(ps);
					}
				}
				
				ps.flush();ps.close();
				
			}
		}
		
		for(String pep : toRemove)
			specPairMap.remove(pep);
	}
}
