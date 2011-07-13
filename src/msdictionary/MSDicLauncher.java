package msdictionary;

import java.util.ArrayList;
import java.util.Iterator;

import msgf.AminoAcidGraph;
import msgf.GeneratingFunction;
import msgf.NominalMassFactory;
import msgf.ScoredSpectrum;
import msgf.ToolLauncher;
import msgf.NominalMass;
import msscorer.NewAdditiveScorer;
import msscorer.NewScorerFactory;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.Peptide;
import msutil.Spectrum;
import suffixarray.SuffixArray;

/**
 * MSDicLauncher class for launching MS-Dictinoary. 
 * @author sangtaekim
 */
public class MSDicLauncher extends ToolLauncher {
	private final SuffixArray sa;
	
	// Determine dictionary size
	private float numRecs = 1e6f;
	private boolean isNumInclusive = false;
	
	/**
	 * A constructor specifies spectral file name and database file name. Database must be "fasta" format.
	 * @param specIterator spectra iterator.
	 * @param sa  a suffix array indexing protein database.
	 * @param scorer a scorer object.
	 */
	public MSDicLauncher(Iterator<Spectrum> specIterator, NewAdditiveScorer scorer, SuffixArray sa)
	{
		super(specIterator, scorer);
		this.sa = sa;
	}

	public MSDicLauncher(Iterator<Spectrum> specIterator, NewAdditiveScorer scorer)
	{
		this(specIterator, scorer, null);
	}
	
	/**
	 * A builder method to set the number of reconstructions searched. 
	 * @param numRecs the number of reconstructions.
	 * @return this object.
	 */
	public MSDicLauncher numRecs(float numRecs)		{ this.numRecs = numRecs; return this; }
	
	/**
	 * If this method is called, it is guaranteed that #Reconstructions is at least numRecs.
	 * @return this object.
	 */
	public MSDicLauncher setNumInclusive()	{ this.isNumInclusive = true; return this; }
	
	public void runMSDictionary()
	{
	    int numSpecs = 0;
	    int numProcessed = 0;
	    
		// TODO: aaSet must be dynamically defined. 
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		// TODO: enzyme must be dynamically defined.
		Enzyme enzyme = Enzyme.TRYPSIN;
		
		int maxLength = 30;
		
		int dataStructure = 0;	// 0: amino acid graph with 1Da bins, 1: composition graph
		// Currently, composition graph is not supported.
//		if(this.fragTolerance.getToleranceAsDa(1000f) < 0.05f)
//			dataStructure = 1;
			
		if(this.sa != null)
			out.println("Title\tScanNum\tPeptide\tProtein\tPrecursorMass\tCharge\tMSGFScore\tPeptideScore\tSpecProb");
		
		NominalMassFactory factory = new NominalMassFactory(aaSet, enzyme, maxLength);
	    while(specIterator.hasNext())
	    {
	    	Spectrum spec = specIterator.next();
	    	numSpecs++;
	    	
	    	int charge = spec.getCharge();

	    	// for charge 0, we assume it is charge 2
	    	if(charge == 0)
	    		spec.setCharge(2);

	    	if(spec.getParentMass() > maxParentMass || spec.getParentMass() < minParentMass)
	    		continue;

			GeneratingFunction<?> gf = null;
			GeneratingFunction<?> gfWellCleaved = null;
			
	    	if(dataStructure == 0)
	    	{
	    		ScoredSpectrum<NominalMass> scoredSpec = NewScorerFactory.get(ActivationMethod.CID, Enzyme.TRYPSIN).getScoredSpectrum(spec);
	    		AminoAcidGraph graph = new AminoAcidGraph(factory, spec.getParentMass(), scoredSpec);
	    		
				gf = new GeneratingFunction<NominalMass>(graph);
				gfWellCleaved = new GeneratingFunction<NominalMass>(graph).enzyme(Enzyme.TRYPSIN);
				gfWellCleaved.doNotBacktrack().doNotCalcNumber();
	    	}
			
	    	numProcessed++;
			
			gf.computeGeneratingFunction();
			
			if(gf.getMaxScore() <= this.msgfScoreThreshold)
				continue;
			
			ArrayList<String> dictionary = gf.getReconstructions(specProb, numRecs, isNumInclusive, sa);

			if(dictionary.size() == 0)
				continue;

	    	if(sa == null)
	    	{
	    		String title = spec.getTitle();
	    		if(title == null)
	    			title = "Spec:"+spec.getSpecIndex();
	    		out.println("#"+spec.getTitle()+"\t"+spec.getSpecIndex()+"\t"+spec.getPrecursorPeak().getMz()+"\t"+charge+"\t"+(gf.getMaxScore()-1));
	    		for(String s : dictionary)
	    		{
	    			Peptide p = new Peptide(s, aaSet);
	    			int score = gf.getGraph().getScore(p);
	    			out.println(s+"\t"+score+"\t"+gf.getSpectralProbability(score));
//	    			out.println(s+"\t"+score);
	    		}
	    	}
	    	else
	    	{
				for(String s : dictionary)	// db search
				{
					ArrayList<String> matchedPeptides = sa.getAllMatchedStrings(s);
					ArrayList<String> matchedProteinAnnotations = sa.getAllMatchingAnnotations(s);
					assert(matchedPeptides.size() == matchedProteinAnnotations.size());
					for(int i=0; i<matchedPeptides.size(); i++)
					{
						Peptide pep = new Peptide(s, aaSet);
						
						int peptideScore = gf.getGraph().getScore(pep);
						
						if(this.trypticOnly && !Enzyme.TRYPSIN.isCleaved(pep))
							continue;
						
						double specProb = gf.getSpectralProbability(peptideScore);
						String id = matchedPeptides.get(i);
						
						out.println(spec.getTitle()+"\t"+spec.getSpecIndex()+"\t"+id+"\t"+matchedProteinAnnotations.get(i)+"\t"+spec.getPrecursorPeak().getMz()+"\t"+spec.getCharge()+"\t"+(gf.getMaxScore()-1)+
								"\t"+peptideScore+"\t"+specProb);
					}
				}	    		
	    	}
	    }
	    out.close();
	}
}
