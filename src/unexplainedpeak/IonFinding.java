package unexplainedpeak;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.Peak;
import msutil.Spectrum;

public class IonFinding {
	
	static HashMap<IonType, Integer> numExplainedPeaks;
	static HashMap<IonType, Integer> iterationNum; 
	static IonDistribution ionDist;
	static int ionDistBinNum = 100;
	static AminoAcidSet aaSet;
	static int numTotalPeaks;
	
	static  public IonDistribution getIonDistribution() {return ionDist;}
	static void setIonDistributionBinNum(int n) {ionDistBinNum = n;}
	
	static public void run(String filename, int rankLimit, float probabiltiy, Tolerance tolerance, AminoAcidSet aas, boolean writemgf, boolean calIonDists, boolean considerAAMassDiff) throws IOException{
		numExplainedPeaks = new HashMap<IonType, Integer>();
		iterationNum =  new HashMap<IonType, Integer>(); 
		numTotalPeaks = 0;
		
		String offsetFilenamePrefix = filename.substring(0,filename.lastIndexOf('.')) + "_off";
		String sigIonFileName = filename.substring(0,filename.lastIndexOf('.')) + ".sig";
		aaSet = aas;
		
		ArrayList<IonType> consideredIonTypes = new ArrayList<IonType>();
		if(calIonDists) ionDist = new IonDistribution(ionDistBinNum);
		
		int numIteration = 0, numSpec, numIons = 0;
		
		ArrayList<Spectrum> spectra = Util.getSpectra(filename, 2);
		numSpec = spectra.size();
	
		
		while(true){
			int numIonsBefore = numIons;
			OffsetFrequencyFunction offset = new OffsetFrequencyFunction(spectra, rankLimit, tolerance);
			
			for(IonType io : offset.getSignificantIons(probabiltiy))
			{
				if(!consideredIonTypes.contains(io)){
					consideredIonTypes.add(io);
					iterationNum.put(io, numIteration);
				}
			}
			numIons = consideredIonTypes.size();
			
			if(numIonsBefore == numIons) break;
			
			AnnotatedPeak.setConsideredIons(consideredIonTypes);
			spectra = getInconsistentSpectraOnlyWithUnexplainedPeaks(spectra, rankLimit, tolerance, calIonDists, considerAAMassDiff);
		
			offset.normalizeOffsetPeaks(numSpec * rankLimit);
			offset.printAll(offsetFilenamePrefix+ "_"+numIteration+".m");
			System.out.println((numIteration+1) + " - Num inconsistent spec : " + spectra.size());
			numIteration++;
		}
		
		writeSigIonFile(sigIonFileName, numSpec, spectra.size(), rankLimit);
			
		if(writemgf) writeSpecFiles(filename, spectra);
		if(calIonDists) ionDist.write(filename);
	}
	
	static void writeSpecFiles(String filename, ArrayList<Spectrum> spectra) throws IOException{
		String inconsistSpecFilename =  filename.substring(0,filename.lastIndexOf('.')) + "_inconsist.mgf";
		String inconsistSpecFilenameOnlyWithUnexplainedPeaks =  filename.substring(0,filename.lastIndexOf('.')) + "_inconsist_Unexp.mgf";
		String inconsistSpecFilenameOnlyWithExplainedPeaks =  filename.substring(0,filename.lastIndexOf('.')) + "_inconsist_Exp.mgf";
		
		ArrayList<Spectrum> spectraToWrite = Util.getSpectra(filename, spectra, 2);
		writeSpectra(inconsistSpecFilename, spectraToWrite);
		writeSpectra(inconsistSpecFilenameOnlyWithUnexplainedPeaks, spectra);
		writeSpectra(inconsistSpecFilenameOnlyWithExplainedPeaks, getSpecWithExplainedPeaks(spectraToWrite, spectra));
	
	}

	
	static HashMap<Character, Float> getAminoAcidRatio(ArrayList<Spectrum> spectra){
		HashMap<Character, Float> ratio = new HashMap<Character, Float>();
		for(Spectrum spec : spectra){
			String annotation = spec.getAnnotationStr();
			HashSet<Character> aaCharSet = new HashSet<Character>();
			for(int i=0; i<annotation.length();i++){
				aaCharSet.add(annotation.charAt(i));
			}
			for(char aa : aaCharSet){
				float num = 0;
				if(ratio.containsKey(aa)) num = ratio.get(aa);
				num++;
				ratio.put(aa, num);
			}
				
		}
		
		for(char aa : ratio.keySet()){
			float num = ratio.get(aa)/spectra.size();
			ratio.put(aa, num);
		}
		
		return ratio;
	}
	
	static void writeSpectra(String filename, ArrayList<Spectrum> spectra) throws IOException{
		PrintStream mgfout = new PrintStream(filename);
		for(Spectrum spec : spectra)
			spec.outputMgf(mgfout);
		mgfout.close();
	}
	
	
	static ArrayList<Spectrum> getSpecWithExplainedPeaks(ArrayList<Spectrum> origspectra, ArrayList<Spectrum> spectraWithUnexplainedPeaks){
		
		for(int i=0; i<origspectra.size();i++){
			Spectrum spec = spectraWithUnexplainedPeaks.get(i);
			Spectrum origSpec = origspectra.get(i);
	
			HashSet<Float> mzs = new HashSet<Float>();
			
			for(Peak p:spec)
				mzs.add(p.getMz());
			
			for(int j=0; j<origSpec.size(); j++){
				Peak p = origSpec.get(j);
				if(mzs.contains(p.getMz())) origSpec.remove(j--);
			}
		}
		
		return origspectra;
		
	}
	
	public static ArrayList<IonType> getOrderedIonTypes(){
		ArrayList<Integer> numPeaks = new ArrayList<Integer>();
		for(IonType ion : numExplainedPeaks.keySet()){
			numPeaks.add(numExplainedPeaks.get(ion));
		}
		
		Collections.sort(numPeaks);
		ArrayList<IonType> ret = new ArrayList<IonType>();
		while(!numPeaks.isEmpty()){
			for(IonType ion : numExplainedPeaks.keySet()){
				if(ret.contains(ion)) continue;
				if(numExplainedPeaks.get(ion) == numPeaks.get(numPeaks.size()-1)){
					numPeaks.remove(numPeaks.size()-1);
					ret.add(ion);
					break;
				}	
			}
		}
		
		return ret;
		
	}
	
	//TODO use getOrderedIonTypes()
	static void writeSigIonFile(String filename, int numSpec, int numInconsistentSpec, int rankLimit) throws IOException{
		ArrayList<Integer> numPeaks = new ArrayList<Integer>();
		BufferedWriter sigOut = new BufferedWriter(new FileWriter(filename));
		System.out.println("Num Spec: " + numSpec + "\tNum Inconsistent Spec: " + numInconsistentSpec +"\tNum Peaks considered: " + numTotalPeaks);
		sigOut.write("Num Spec: " + numSpec + "\tNum Inconsistent Spec: " + numInconsistentSpec +"\tNum Peaks considered: " + numTotalPeaks	+"\n");
		
		for(IonType ion : numExplainedPeaks.keySet()){
			numPeaks.add(numExplainedPeaks.get(ion));
		}
		
		Collections.sort(numPeaks);
		HashSet<IonType> written = new HashSet<IonType>();
		while(!numPeaks.isEmpty()){
			for(IonType ion : numExplainedPeaks.keySet()){
				if(written.contains(ion)) continue;
				if(numExplainedPeaks.get(ion) == numPeaks.get(numPeaks.size()-1)){
					numPeaks.remove(numPeaks.size()-1);
					String s = (iterationNum.get(ion)+1) +"&$"+Util.getIonString(ion) +"$&" + numExplainedPeaks.get(ion) + "&" + String.format("%.1f", (float)numExplainedPeaks.get(ion)/(float)numTotalPeaks*100) + "\\%\\\\";
					System.out.println(s);
					written.add(ion);
					sigOut.write(s+"\n");
					break;
				}	
			}
		}
		
		sigOut.close();
	}
	
	static ArrayList<Spectrum> getInconsistentSpectraOnlyWithUnexplainedPeaks(ArrayList<Spectrum> spectra, int rankLimit, Tolerance tolerance, boolean calIonDists,boolean considerAAMassDiff){
		ArrayList<Spectrum> unexpSpectra = new ArrayList<Spectrum>();
		
		boolean calNumTotalpeaks = numTotalPeaks == 0;
		int[] sum = new int[2];
		for(Spectrum spec : spectra){
			Spectrum unexpSpec = spec.getCloneWithoutPeakList();
			boolean isThisSpecConsistent = true;
			ArrayList<Peak> peaks = null;
			ArrayList<IonType> ions = new ArrayList<IonType>();
			
			peaks = Util.pickPeaks(spec, rankLimit, considerAAMassDiff, aaSet, tolerance);
			if(calNumTotalpeaks) numTotalPeaks += peaks.size();
			for(Peak peak : peaks){
				AnnotatedPeak annotatedPeak = new AnnotatedPeak(peak.getMz(), peak.getIntensity(), peak.getCharge(), spec.getAnnotation(), tolerance);
				if(!annotatedPeak.isExplained()){
					isThisSpecConsistent = false;
					unexpSpec.add(peak);
				}else{
					for(int i=0; i<annotatedPeak.getExplainingIonTypes().size(); i++){
						IonType ion = annotatedPeak.getExplainingIonTypes().get(i);
						ions.add(ion);
						int n;
						if(numExplainedPeaks.containsKey(ion)){
							n = numExplainedPeaks.get(ion)+1;
						}else{
							n = 1;
						}
						numExplainedPeaks.put(ion, n);
						if(calIonDists) ionDist.add(ion, spec.getAnnotation(), annotatedPeak.getMz());
						
					}
				}
			}
			if(ions.size() > 1){
				if(ions.get(0).equals(ions.get(1))){
					sum[0]++;
					//System.out.println("same");
				}
				else{
					sum[1]++; 
				//	System.out.println("****false" + ions.get(0) + "\t" + ions.get(1));
				}
			}
			if(!isThisSpecConsistent) unexpSpectra.add(unexpSpec);
		}
		System.out.println(sum[0] + "\t" + sum[1]);
		return unexpSpectra;
	}
	
	
}
