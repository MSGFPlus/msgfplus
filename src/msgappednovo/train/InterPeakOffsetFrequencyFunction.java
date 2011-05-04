package msgappednovo.train;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import msgf.Tolerance;
import msutil.Peak;

public class InterPeakOffsetFrequencyFunction {
	static public class NewGeneralizedOffsetPeak  implements Comparable<NewGeneralizedOffsetPeak>{
    	private InterPeakOffset gof;
    	private float y;
    	
    	NewGeneralizedOffsetPeak(InterPeakOffset gof, float y){
    		this.gof = gof; this.y = y;
    	}
    	
    	public int compareTo(NewGeneralizedOffsetPeak o) {
    		return new Float(this.gof.getOffset()).compareTo(new Float(o.gof.getOffset()));
    	}
    	
    	public float getProbability() { return y; }
    	public InterPeakOffset getGeneralizedOffset() { return gof; }
    }
	
	static public final float MAX = 38;
	static public final float MIN = -38;
	
/*	static public ArrayList<NewGeneralizedOffsetPeak> getOffSetFrequencyFunction(HashMap<InterPeakOffset, Integer> offsetnums, float normalizer){
		return getOffSetFrequencyFunction(offsetnums, normalizer, 0);
	}
	
	static public ArrayList<NewGeneralizedOffsetPeak> getOffSetFrequencyFunction(HashMap<InterPeakOffset, Integer> offsetnums, float normalizer, float threshold){
		return getOffSetFrequencyFunction(offsetnums, normalizer, threshold, null);
	}
	*/
	static public ArrayList<NewGeneralizedOffsetPeak> getOffSetFrequencyFunction(HashMap<InterPeakOffset, Integer> offsetnums, float normalizer, float threshold, String filename){
		ArrayList<NewGeneralizedOffsetPeak> offsetPeaks = new  ArrayList<NewGeneralizedOffsetPeak>();
		ArrayList<NewGeneralizedOffsetPeak> offsetPeaksforOutput = new  ArrayList<NewGeneralizedOffsetPeak>();
		if(offsetnums == null || offsetnums.isEmpty()) return offsetPeaks;
		
		for(InterPeakOffset gof : offsetnums.keySet()){
			float prob = (float)offsetnums.get(gof)/normalizer;
			
			if(prob > threshold)
				offsetPeaks.add(new NewGeneralizedOffsetPeak(gof, prob));
			
			if(filename != null) offsetPeaksforOutput.add(new NewGeneralizedOffsetPeak(gof, prob));
		}
		//for(NewGeneralizedOffsetPeak d : offsetPeaks)
		//	System.out.println(d.getGeneralizedOffset());
		Collections.sort(offsetPeaks);
		if(filename != null) Collections.sort(offsetPeaksforOutput);
		
		if(filename != null && !offsetPeaksforOutput.isEmpty()){
			try {
				HashSet<Integer> cs = new HashSet<Integer>();
				HashSet<Integer> cos = new HashSet<Integer>();
				HashSet<Boolean> comps = new HashSet<Boolean>();
				for(NewGeneralizedOffsetPeak p : offsetPeaksforOutput){
					cs.add(p.gof.getBaseCharge());
					cos.add(p.gof.getChargeOffset());
					comps.add(p.gof.isComplementary());
				}
				
				for(int c : cs){
					for(int co : cos){
						for(boolean com : comps){
							boolean towrite = false;
							for(NewGeneralizedOffsetPeak p : offsetPeaksforOutput){
								if(p.gof.getBaseCharge() == c && p.gof.getChargeOffset() == co && p.gof.isComplementary() == com)
									if(p.y > threshold){
										towrite = true;
										break;
									}
							}
							
							if(!towrite) continue;
							
							PrintStream out = new PrintStream(filename + "charge: " + c +  " chargeOff: " + co + " comp: " + com);
							out.println("off=[");
							for(NewGeneralizedOffsetPeak p : offsetPeaksforOutput){
								if(p.gof.getBaseCharge() == c && p.gof.getChargeOffset() == co && p.gof.isComplementary() == com)
									out.println(p.gof.getOffset() + "\t" + p.y);
							}
							out.println("];");
							out.close();
						}
					}
				}
				
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		return offsetPeaks;
	}
	
	
	
	
		
	static public ArrayList<InterPeakOffset> getGeneralizedOffsets(Peak bp, ArrayList<Peak> cps, boolean isComplementary, int chargeOffset, Tolerance tol){
		ArrayList<InterPeakOffset> offs = new ArrayList<InterPeakOffset>();
		for(Peak cp : cps){
			InterPeakOffset off = getGeneralizedOffset(bp, cp, isComplementary, chargeOffset, tol);
			if(off != null) offs.add(off);
		}
		
		return offs;
	}
	
	static public ArrayList<InterPeakOffset> getGeneralizedOffsets(ArrayList<Peak> bps, ArrayList<Peak> cps, boolean isComplementary, int chargeOffset, Tolerance tol){
		ArrayList<InterPeakOffset> offs = new ArrayList<InterPeakOffset>();
		for(Peak bp : bps){
			for(Peak cp : cps){
				InterPeakOffset off = getGeneralizedOffset(bp, cp, isComplementary, chargeOffset, tol);
				if(off != null) offs.add(off);
			}
		}
		return offs;
	}
	
	static public InterPeakOffset getGeneralizedOffset(Peak bp, Peak cp, boolean isComplementary, int chargeOffset, Tolerance tol){
		// base peaks already have charge offset	
		float offset = (cp.getMz() - bp.getMz()) * bp.getCharge();
	//	if(bp.getCharge() == 2) System.out.println(offset);
		if(offset > MAX || offset < MIN) return null;
		return new InterPeakOffset(offset, isComplementary, chargeOffset, bp.getCharge() - chargeOffset, tol.getToleranceAsDa(500)*2);
	}
	
}
