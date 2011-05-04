package unexplainedpeak;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import msutil.Peak;
import msutil.Spectrum;

public class PeakDistribution {
	static int minPM = 300, maxPM = 3000;
	static int pmResolution = 100;
	static int numPMBins = (maxPM-minPM)/pmResolution;
	static int numBins = 100;
	
	public static void write(String filename) throws IOException{
		String outfilename = filename.substring(0,filename.lastIndexOf('.')) + "_peakDist.m";
		BufferedWriter out = new BufferedWriter(new FileWriter(outfilename));
		
		ArrayList<Spectrum> spectra = Util.getSpectra(filename, 2);
		float[][] dist = new float[numPMBins+1][numBins+1];
	//	System.out.println(spectra.size());
		for(Spectrum spec : spectra){
			float pm = spec.getParentMass();
			if(pm < minPM || pm > maxPM) continue;
			
			for(Peak peak : spec){
				dist[Math.round((pm - minPM)/pmResolution)][Math.round((peak.getMz()/pm) * numBins)]++;
			}
		}
		
		String s="";
		int sum=0;
		for(int i=0;i<=numPMBins;i++){
			s+="pm(:," + (i+1) + ")=[";
			for(int j=0; j<=numBins;j++){
				s+= dist[i][j] + "\t";
				sum += dist[i][j];
			}
			s+= "];\n";
		}
		out.write(s);
		//System.out.println(sum);
		out.close();
	}
}
