package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Hashtable;

import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.MzXMLSpectraMap;


public class HeckPercolator {

	private static Hashtable<String, Float[]> mascotThresholds;
	private static Hashtable<String, Float[]> percolatorThresholds;
	private static Hashtable<String, Float[]> msgfThresholds;

	public static void main(String argv[]) throws Exception
	{
		makeAnnotatedSpectraWithMascot();
	}
	
	public static void makeAnnotatedSpectraWithMascot() throws Exception
	{
		File specDir = new File("/home/sangtaekim/Research/Data/HeckWhole/Spectra");
		HashMap<String, MzXMLSpectraMap> specMap = new HashMap<String, MzXMLSpectraMap>();
		for(File f : specDir.listFiles())
		{
			String fileName = f.getName();
			if(!fileName.endsWith(".mzXML"))
				continue;
			String name = fileName.substring(0, fileName.lastIndexOf('.')).toLowerCase();
//			System.out.println(name);
			specMap.put(name, new MzXMLSpectraMap(f.getPath()));
		}
		
		File mascotDir = new File("/home/sangtaekim/Research/Data/HeckRevision/Mascot23");
		File annotatedSpecDir = new File("/home/sangtaekim/Research/Data/HeckRevision/AnnotatedSpectra");
		
		
		for(File f : mascotDir.listFiles())
		{
			String fileName = f.getName();
			if(!fileName.endsWith("Target.txt"))
				continue;
			String method = fileName.substring(0, 7);
			String outputFileName = annotatedSpecDir.getPath() + File.separator + "MascotAnnotated_" + method + ".mgf";
			PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
			
			System.out.print(fileName+": ");
			
			HashMap<String, Spectrum> pepSpecMap = new HashMap<String, Spectrum>();	// peptide -> spectrum
			HashMap<String, Float> pepScoreMap = new HashMap<String, Float>();	// peptide -> score
			String s;
			BufferedLineReader in = new BufferedLineReader(f.getPath());
			String prevTitle = "";
			while((s=in.readLine()) != null)
			{
				if(s.startsWith("#"))
					continue;
				String[] token = s.split("\t");
				if(token.length < 4)
					continue;
				if(prevTitle.equalsIgnoreCase(token[0]))
					continue;
				else
					prevTitle = token[0];
				float score = Float.parseFloat(token[3]);
				int chargeIndex = Integer.parseInt(token[1]) - 2;
				if(chargeIndex > 2)
					chargeIndex = 2;
				float threshold = mascotThresholds.get(method)[chargeIndex]; 
				
				if(score > threshold)
				{
					String annotation = token[2];
					String peptide = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
					if(pepSpecMap.get(peptide+":"+chargeIndex) == null || pepScoreMap.get(peptide+":"+chargeIndex) < score)
					{
						// get the spectrum
						String specFileID = token[0].substring(token[0].indexOf("090"), token[0].lastIndexOf('_'));
						String[] titleToken = token[0].split("\\s+");
						int scanNum = Integer.parseInt(titleToken[titleToken.length-1]);
						Spectrum spec = specMap.get(specFileID.toLowerCase()).getSpectrumBySpecIndex(scanNum);
						Peptide pep = new Peptide(peptide);
						spec.setAnnotation(pep);
						assert(spec != null);
						assert(Math.abs(spec.getParentMass() - pep.getParentMass()) < 5f): s + " " + spec.getParentMass() + " " + pep.getParentMass();
//						if(Math.abs(spec.getParentMass() - pep.getParentMass()) > 0.1f)
//							System.out.println(s + " " + spec.getParentMass() + " " + pep.getParentMass());
						pepSpecMap.put(peptide+":"+chargeIndex, spec);
						pepScoreMap.put(peptide+":"+chargeIndex, score);
					}
				}
			}
			
			int numSpecs = 0;
			for(Spectrum spec : pepSpecMap.values())
			{
				spec.outputMgf(out);
				numSpecs++;
			}
			
			in.close();
			out.close();
			System.out.println(numSpecs);
		}
	}
	
	static {
		Float[]	mascotThresholdCIDTryp = {37.08f, 23.35f, 24.18f};
		Float[]	mascotThresholdETDTryp = {46.79f, 74.49f, 44.2f};
		Float[]	mascotThresholdCIDLysN = {33.53f, 22.31f, 23.01f};
		Float[]	mascotThresholdETDLysN = {35.8f, 30.48f, 21.07f};
		mascotThresholds = new Hashtable<String, Float[]>();		
		mascotThresholds.put("CIDTryp", mascotThresholdCIDTryp);
		mascotThresholds.put("ETDTryp", mascotThresholdETDTryp);
		mascotThresholds.put("CIDLysN", mascotThresholdCIDLysN);
		mascotThresholds.put("ETDLysN", mascotThresholdETDLysN);
		
		Float[]	percolatorThresholdCIDTryp = {0.0239881f, 0.0467094f, 0.0178952f};
		Float[]	percolatorThresholdETDTryp = {0.0147269f, 0.050837f, 0.0568241f};
		Float[]	percolatorThresholdCIDLysN = {0.0295369f, 0.0437894f, 0.0125539f};
		Float[]	percolatorThresholdETDLysN = {0.0207123f, 0.0456574f, 0.0456574f};
		percolatorThresholds = new Hashtable<String, Float[]>();		
		percolatorThresholds.put("CIDTryp", percolatorThresholdCIDTryp);
		percolatorThresholds.put("ETDTryp", percolatorThresholdETDTryp);
		percolatorThresholds.put("CIDLysN", percolatorThresholdCIDLysN);
		percolatorThresholds.put("ETDLysN", percolatorThresholdETDLysN);
		
		Float[]	thresholdCIDTryp = {4.9982934E-11f, 2.5958782E-11f, 1.0835081E-11f};
		Float[]	thresholdETDTryp = {2.3154301E-11f, 7.760064E-13f, 2.299933E-13f};
		Float[]	thresholdSumTryp = {5.6320247E-11f, 2.3743533E-11f, 1.8667672E-11f};
		Float[]	thresholdCIDLysN = {3.2987557E-11f, 1.5884185E-11f, 1.615089E-12f};
		Float[]	thresholdETDLysN = {2.2711777E-11f, 2.1367989E-11f, 9.767584E-12f};
		Float[]	thresholdSumLysN = {2.4446337E-11f, 2.8512514E-11f, 9.074788E-12f};
		msgfThresholds = new Hashtable<String, Float[]>();		
		msgfThresholds.put("CIDTryp", thresholdCIDTryp);
		msgfThresholds.put("ETDTryp", thresholdETDTryp);
		msgfThresholds.put("SumTryp", thresholdSumTryp);
		msgfThresholds.put("CIDLysN", thresholdCIDLysN);
		msgfThresholds.put("ETDLysN", thresholdETDLysN);
		msgfThresholds.put("SumLysN", thresholdSumLysN);
	}
	
}
