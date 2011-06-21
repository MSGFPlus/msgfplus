package misc;

import parser.MgfSpectrumParser;
import msgf.DeNovoGraph;
import msgf.FlexAminoAcidGraph;
import msgf.GeneratingFunction;
import msgf.NominalMass;
import msscorer.NewRankScorer;
import msscorer.NewScoredSpectrum;
import msscorer.NewScorerFactory;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.SpectraIterator;
import msutil.Spectrum;

public class ChargePrediction {
	public static void main(String argv[]) throws Exception
	{
		testChargePrediction();
	}
	
	public static void testChargePrediction() throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		ActivationMethod method = ActivationMethod.CID;
		Enzyme enzyme = Enzyme.TRYPSIN;
		NewRankScorer scorer = NewScorerFactory.get(method, enzyme);
		
		String fileName = "/home/sangtaekim/Research/Data/HeckRevision/AnnotatedSpectra/CID_Tryp_Confident.mgf";
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		
		int numPredictions = 0;
		int numCorrect = 0;
		int minCharge = 2;
		int maxCharge = 4;
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int charge = spec.getCharge();
			if(charge < minCharge || charge > maxCharge)
				continue;

			numPredictions++;
			System.out.print(spec.getAnnotationStr()+"\t"+charge);
			int bestScore = Integer.MIN_VALUE;
			int bestCharge = -1;
			for(int c=minCharge; c<=maxCharge; c++)
			{
				spec.setCharge(c);
				int nominalPeptideMass = NominalMass.toNominalMass(spec.getPeptideMass());
				NewScoredSpectrum<NominalMass> scoredSpec = scorer.getScoredSpectrum(spec);
				DeNovoGraph<NominalMass> graph = new FlexAminoAcidGraph(
						aaSet, 
						nominalPeptideMass,
						enzyme,
						scoredSpec
						);
				
				GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(graph)
				.doNotBacktrack()
				.doNotCalcNumber()
				.doNotCalcProb();
				gf.computeGeneratingFunction();
				System.out.print("\t"+(gf.getMaxScore()-1));
				if(bestScore < gf.getMaxScore()-1)
				{
					bestScore = gf.getMaxScore()-1;
					bestCharge = c;
				}
			}
			if(charge == bestCharge)
			{
				System.out.println("\t1");
				numCorrect++;
			}
			else
			{
				System.out.println("\t0");
			}
		}
		System.out.println(numCorrect/(float)numPredictions);
	}
}
