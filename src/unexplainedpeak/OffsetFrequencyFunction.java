package unexplainedpeak;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import parser.MgfSpectrumParser;
import msgf.Tolerance;
import msscorer.PrecursorOffsetFrequency;
import msutil.Composition;
import msutil.Constants;
import msutil.IonType;
import msutil.Peak;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;

public class OffsetFrequencyFunction {
	// prefix, suffix, internal, precursor, charge 1, 2 for now.
    class offsetPeak  implements Comparable<offsetPeak>{
    	float x;
    	float y;
    	offsetPeak(float x, float y){
    		this.x = x; this.y = y;
    	}
    	
    	public int compareTo(offsetPeak o) {
    		return new Float(this.x).compareTo(new Float(o.x));
    	}
    }
	private int rankLimit;
	private int maxCharge = 2;
	private Tolerance tolerance;
	private int min = -60, max = 60;
	private int sigmin = -56, sigmax = 56;
	
	private HashMap<String, ArrayList<offsetPeak>> offsetPeakMap = new HashMap<String, ArrayList<offsetPeak>>();
	private HashMap<String, ArrayList<Float>> offsetMap = new HashMap<String, ArrayList<Float>>();
	
	private void addOffsetMap(Spectrum spec, String type, int charge){
		String key = type + charge;
		ArrayList<Float> toput;
		
		if(!offsetMap.containsKey(key))
			toput = new ArrayList<Float>();
		else
			toput = offsetMap.get(key);
		
		if(type.equals("P"))
			toput.addAll(getPrefixOffset(spec, charge));
		else if(type.equals("S"))
			toput.addAll(getSuffixOffset(spec, charge));
		else if(type.equals("I"))
			toput.addAll(getInternalOffset(spec, charge));
		else if(type.equals("C"))
			toput.addAll(getCyclicOffset(spec, charge));
		else if(type.equals("R"))
			toput.addAll(getPrecursorOffset(spec, charge));
		
		offsetMap.put(key, toput);
	}
	
	private void addAllTypeOffset(Spectrum spec){
		for(int charge=1; charge<=spec.getCharge(); charge++){
			addOffsetMap(spec, "P", charge);
			addOffsetMap(spec, "S", charge);
			addOffsetMap(spec, "I", charge);
			addOffsetMap(spec, "R", charge);
			addOffsetMap(spec, "C", charge);
		}
	}
	
	
	public OffsetFrequencyFunction(String specFileName, int rankLimit, Tolerance tolerance) throws IOException{
		this.rankLimit = rankLimit;
		this.tolerance = tolerance;
		Iterator<Spectrum> iterator = new SpectraIterator(specFileName, new MgfSpectrumParser());
		
		System.out.println("Generating OFF");
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			if(spec.getAnnotation().isModified()) continue;
			System.out.println(spec.getScanNum());
			spec.setRanksOfPeaks();
			addAllTypeOffset(spec);			
		}
		
		for(String key : offsetMap.keySet()){
			offsetPeakMap.put(key, getOffSetFrequencyFunction(offsetMap.get(key), key));
		}
		System.out.println("Done");
		
	}
	
	public OffsetFrequencyFunction(ArrayList<Spectrum> spectra, int rankLimit, Tolerance tolerance){
		this.rankLimit = rankLimit;
		this.tolerance = tolerance;
		System.out.println("Generating OFF");
		for(Spectrum spec : spectra){
			if(spec.getAnnotation().isModified()) continue;
			//System.out.println(spec.getScanNum());
			spec.setRanksOfPeaks();
			addAllTypeOffset(spec);			
		}
		
		for(String key : offsetMap.keySet()){
			offsetPeakMap.put(key, getOffSetFrequencyFunction(offsetMap.get(key), key));
		}
		System.out.println("Done");
	}
	
	
	

	private ArrayList<offsetPeak> getOffSetFrequencyFunction(ArrayList<Float> offsets, String key){
		ArrayList<offsetPeak> offsetPeaks = new  ArrayList<offsetPeak>();
		int charge = Integer.parseInt(key.substring(1));
		
		float resolution = tolerance.getToleranceAsDa(57)*2/charge;
		
		int[] peaks = new int[(int)((max - min) / resolution)+1];
		for(float offset : offsets){
			peaks[Math.round((float)(offset - min) / resolution)]++; 	
		}
		
		for(int i=0; i<peaks.length; i++){
			if(peaks[i] > 0)
				offsetPeaks.add(new offsetPeak(i*resolution + min, peaks[i]));
		}
	
		return offsetPeaks;
		
	}
	
	private ArrayList<Float> getPrefixOffset(Spectrum spec, int charge){
		ArrayList<Float>  prefixOffset = new ArrayList<Float> ();
		Peptide annotation = spec.getAnnotation();
		
		float[] prefixMasses = annotation.getPRMMasses(true, 0);
		for(Peak peak : spec){
			if(peak.getRank() <= rankLimit){
				for(float prefixMass : prefixMasses){
					float offset = (float) (peak.getMz() - prefixMass / charge);
					if(offset >= min && offset <= max)
						prefixOffset.add(offset);
				}
			}
		}
		
		return prefixOffset;
	}
	
	private ArrayList<Float> getSuffixOffset(Spectrum spec, int charge){
		ArrayList<Float>  suffixOffset = new ArrayList<Float> ();
		Peptide annotation = spec.getAnnotation();
		
		float[] suffixMasses = annotation.getPRMMasses(false, 0);
		for(Peak peak : spec){
			if(peak.getRank() <= rankLimit){
				for(float suffixMass : suffixMasses){
					float offset = (float) (peak.getMz() - suffixMass / charge);
					if(offset >= min && offset <= max)
						suffixOffset.add(offset);
				}
			}
		}
		
		return suffixOffset;
	}
	
	private ArrayList<Float> getInternalOffset(Spectrum spec, int charge){
		ArrayList<Float>  internalOffset = new ArrayList<Float> ();
		Peptide annotation = spec.getAnnotation();
		
		ArrayList<Float> internalMasses = new ArrayList<Float>();
		
		for(int i=1; i<annotation.size()-1; i++){
			for(int j=i+1; j<annotation.size()-1;j++){
				internalMasses.add(annotation.getMass(i, j));
			}
		}
		
		for(Peak peak : spec){
			if(peak.getRank() <= rankLimit){
				for(float internalMass : internalMasses){
					float offset = (float) (peak.getMz() - internalMass / charge);
					if(offset >= min && offset <= max)
						internalOffset.add(offset);
				}
			}
		}
		return internalOffset;
	}
	
	private ArrayList<Float> getCyclicOffset(Spectrum spec, int charge){
		ArrayList<Float>  cyclicOffset = new ArrayList<Float> ();
		Peptide annotation = spec.getAnnotation();
		
		ArrayList<Float> cyclicMasses = new ArrayList<Float>();
		
		for(int i=2; i<annotation.size();i++){
			for(int j=annotation.size()+1;;j++){
				if(j-i >= annotation.size()) break;
				float cyclicMass = annotation.getMass(i, annotation.size());
				cyclicMass += annotation.getMass(0, j-annotation.size());
				cyclicMasses.add(cyclicMass);
			}
		}
		
		for(Peak peak : spec){
			if(peak.getRank() <= rankLimit){
				for(float cyclicMass : cyclicMasses){
					float offset = (float) (peak.getMz() - cyclicMass / charge);
					if(offset >= min && offset <= max)
						cyclicOffset.add(offset);
				}
			}
		}
		return cyclicOffset;
	}
	
	private ArrayList<Float> getPrecursorOffset(Spectrum spec, int charge){
		ArrayList<Float>  precursorOffset = new ArrayList<Float> ();
		
		float precrusorMass = spec.getAnnotation().getParentMass();
		
		for(Peak peak : spec){
			if(peak.getRank() <= rankLimit){
				float offset = (float) (peak.getMz() - precrusorMass / charge);
				if(offset >= min && offset <= max){
					precursorOffset.add(offset);
				}
			}
		}
		
		return precursorOffset;
	}
	

	//TODO make chebyshev thing in Util
	public ArrayList<IonType> getSignificantIons(float probability){
		ArrayList<IonType> ionTypes = new  ArrayList<IonType>();
		for(String key : offsetPeakMap.keySet()){
			
			float avg=0, var = 0;
			for(offsetPeak op : offsetPeakMap.get(key))
				avg += op.y;
			
			avg /= offsetPeakMap.get(key).size();
			
			for(offsetPeak op : offsetPeakMap.get(key))
				var += Math.pow((op.y - avg), 2);
			
			var /= offsetPeakMap.get(key).size()-1;
			
			float threshold = (float)(avg + Math.sqrt(var / probability));
			
			for(offsetPeak op : offsetPeakMap.get(key))
				if(op.y >= threshold){
					int charge = Integer.parseInt(key.substring(1));
					if(!(key.startsWith("R") && charge>1) && (op.x < sigmin/charge || op.x > sigmax/charge)) continue;
					if(key.startsWith("S"))
						ionTypes.add(new IonType.SuffixIon(key, charge, op.x / Constants.INTEGER_MASS_SCALER));
					else if(key.startsWith("I"))
						ionTypes.add(new IonType.InternalIon(key, charge, op.x/ Constants.INTEGER_MASS_SCALER));
					else if(key.startsWith("C"))
						ionTypes.add(new IonType.CyclicIon(key, charge, op.x/ Constants.INTEGER_MASS_SCALER));
					else if(key.startsWith("P"))
						ionTypes.add(new IonType.PrefixIon(key, charge, op.x/ Constants.INTEGER_MASS_SCALER));
					else
						ionTypes.add(new IonType.PrecursorIon(key, charge, op.x/ Constants.INTEGER_MASS_SCALER));
				}
		}
		
		
		return ionTypes;
	}
	
	
	public void normalizeOffsetPeaks(int specNum){
		for(String key : offsetPeakMap.keySet()){
			ArrayList<offsetPeak> offsetPeaks = offsetPeakMap.get(key);
			for(offsetPeak op : offsetPeaks){
				op.y /= (float)specNum;
			}
		}
	}
	
	public void printAll(String fileName) throws IOException{
		new File(fileName).delete();
		for(String key : offsetPeakMap.keySet()){
			print(key, fileName);
		}
	}
	
	public void print(String type, String fileName) throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(fileName, true));
		if(offsetPeakMap.containsKey(type)){
			int charge = Integer.parseInt(type.substring(1));
			ArrayList<offsetPeak> offsetPeaks = offsetPeakMap.get(type);
			if(offsetPeaks.isEmpty()) return;
			Collections.sort(offsetPeaks);
			out.write(type+"=[\n");
			for(offsetPeak op : offsetPeaks)
				out.write((op.x + 1) * charge+"\t"+op.y+"\n"); //TODO mass diff, not mz diff
			out.write("];");
		}
		out.close();
	}
}
