package trex.analysis;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import msgf.AminoAcidGraph;
import msgf.GeneratingFunction;
import msgf.NominalMass;
import msgf.ScoredSpectrum;
import msgf.Tolerance;
import msscorer.NewScorerFactory;
import msutil.ActivationMethod;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.SpectrumRecalibrator;
import parser.MgfSpectrumParser;

public class Test_found {
	public static void main(String[] argv) throws IOException
	{

		//float dictionaryProb;
	
		Tolerance tolerance = new Tolerance(0.5f, false);
		boolean isTolerancePPM = false;
		boolean rescaling = true;
		boolean recalibration = true;
		boolean trypticOnly = true;
		
		String fileName = null;

		String rankScorerFileName = null;
		
		String dirName = "/home/kwj/workspace/";
		
		//fileName = dirName+"inputs/Had/Hadrosaur_14217_Orbitrap_mz_fixed.mgf";
		fileName = dirName+"inputs/TRex/TRex48216.mgf";

		
		SpectrumRecalibrator recalibrator = new SpectrumRecalibrator(tolerance);

		//String[] m = { "GATGApGIAGApGFpGAR",};// "GSNGEpGSAGPpGPAGLR", "GLPGESGAVGPAGPpGSR"};
		String[] m = { "GVVGLpGQR", "GVQGPpGPQGPR", "GLPGESGAVGPAGPIGSR", "GAPGPQGPSGApGPK"};
		
		ArrayList<String> matches = new ArrayList<String>();
		
		for(String s: m)
			matches.add(s);
		
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet("/home/kwj/workspace/inputs/cavebear/AASetWithHydroProline.txt");
		
		SpectraIterator iterator = null;
		
		try {
			iterator = new SpectraIterator(fileName, new MgfSpectrumParser());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		while(iterator.hasNext())
		{
			Spectrum origSpec = iterator.next();
			if(origSpec.getCharge() != 2)
				continue;

		//	Peptide annotation = origSpec.getAnnotation();
		//	if(trypticOnly && !annotation.hasTrypticCTerm()) continue;
				
			origSpec.filterPrecursorPeaks(tolerance);
			
			origSpec.correctParentMass();
			Spectrum spec = null;
			
			if(recalibration){
				spec = recalibrator.recalibrateUsingRescaling(origSpec);
			}else{
				spec = origSpec;
			}
			
			spec.correctParentMass();
			
		//	float mz = spec.getPrecursorPeak().getMz();
		//	if(mz )
			
			
			ScoredSpectrum<NominalMass> scoredSpec = NewScorerFactory.get(ActivationMethod.CID, Enzyme.TRYPSIN).getScoredSpectrum(spec);
			
			AminoAcidGraph graph = new AminoAcidGraph(spec.getParentMass(), aaSet, Enzyme.TRYPSIN);
			
			//GappedGeneratingFunction<IntMass> gap = new GappedGeneratingFunction<IntMass>(scoredSpec, graph, 1e-10f, trypticOnly, 20);
			GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(scoredSpec, graph);
			gf.computeGeneratingFunction();
			//System.out.println(spec.getScanNum());
			for(String s: matches){
				ArrayList<AminoAcid> aaArray = new ArrayList<AminoAcid>();
				for(int i=0; i<s.length();i++){
					char c = s.charAt(i);
					aaArray.add(aaSet.getAminoAcid(c));
				}
				
				Peptide pep = new Peptide(aaArray);
				//System.out.println(pep);
				if(Math.abs(pep.getParentMass() - spec.getParentMass()) > 5||  gf.getScore(pep) < 0) continue;
				
				System.out.print(s+"\t" + pep.getParentMass() + "\t" + spec.getParentMass()+"\t");
				System.out.println(spec.getScanNum() + "\t" + gf.getScore(pep) + "\t" + gf.getSpectralProbability(pep));
			}
				
		}
			

	}
}
