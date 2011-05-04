package msgap;

import java.io.PrintWriter;

import msgf.AminoAcidGraph;
import msgf.DeNovoGraph;
import msgf.NominalMassFactory.NominalMass;
import msgf.ScoreDist;
import msgf.ScoredSpectrum;
import msutil.Spectrum;

public class ScoringParameter {
	private int specID;
	private ScoredSpectrum<NominalMass> scoredSpec;
	private AminoAcidGraph graph;
	private ScoreDist finalDist, finalDistWellCleaved;
	private float correctedParentMass;
	private float originalParentMass;
	private int scanNum;
	private String specfilename;
	private String annotation = null;
	private int charge;
	
	public ScoringParameter(int specID, Spectrum s, String specfilename, GappedGeneratingFunction<NominalMass> gap, float originalParentMass){
		this.specID = specID;
		this.scoredSpec = gap.getScoredSpectrum();
		this.graph = (AminoAcidGraph) gap.getGraph();
		this.correctedParentMass = s.getParentMass();
		this.scanNum = s.getScanNum();
		this.finalDist = gap.getScoreDist();
		this.finalDistWellCleaved = gap.getDistWellCleaved();
		this.specfilename = specfilename;
		this.originalParentMass = originalParentMass;
		this.annotation = s.getAnnotationStr();
		this.charge = s.getCharge();
	}
	
	ScoringParameter(int specID, ScoredSpectrum<NominalMass> scoredSpec, AminoAcidGraph graph, 
			float originalParentMass,
			float correctedParentMass, 
			int scanNum,
			int charge,
			ScoreDist finalDist, 
			ScoreDist finalDistWellCleaved, 
			String specfilename,
			String annotation){
		this.specID = specID;
		this.scoredSpec = scoredSpec;
		this.graph = graph;
		this.correctedParentMass = correctedParentMass;
		this.scanNum = scanNum;
		this.finalDist = finalDist;
		this.finalDistWellCleaved = finalDistWellCleaved;
		this.specfilename = specfilename;
		this.originalParentMass = originalParentMass;
		this.annotation = annotation;
		this.charge = charge;
	}
	
	public int getSpecID() { return specID; }
	public ScoredSpectrum<NominalMass> getScoredSpec() { return scoredSpec; }
	public DeNovoGraph<NominalMass> getGraph() { return graph; }
	public float getCorrectedParentMass() { return correctedParentMass; }
	public float getOriginalParentMass() { return originalParentMass; }
	public int getScanNum() { return scanNum; }
	public ScoreDist getFinalDist() { return finalDist; }
	public ScoreDist getFinalDistWellCleaved() { 
		if(finalDistWellCleaved!=null)
			return finalDistWellCleaved; 
		else 
			return finalDist;
	}
	public String getSpecFileName() { return specfilename; }
	public String getAnnotation() { return annotation; }
	public int getCharge() { return charge; }
	
	/**
	 * Setter method that changes the corrected and original parent masses
	 * @param mass
	 */
	public void setParentMass(float mass) {
	  this.correctedParentMass = this.originalParentMass = mass;
	}
	
	public void outputFIle(PrintWriter out){ // only distributions and corrected parent masses are written
		out.println("#"+ specID + "\t" +specfilename+"\t"+scanNum+"\t" + charge + "\t"+ originalParentMass + "\t" + correctedParentMass + "\t" + annotation);
		out.println("DIST1\t"+finalDist.getMinScore() + "\t" + finalDist.getMaxScore());
		for(int score = finalDist.getMinScore() ; score < finalDist.getMaxScore(); score ++){
			out.println(score+"\t"+finalDist.getProbability(score));
		}
		
		if(finalDistWellCleaved != null){
			out.println("DIST2\t"+finalDistWellCleaved.getMinScore() + "\t" + finalDistWellCleaved.getMaxScore());
			for(int score = finalDistWellCleaved.getMinScore() ; score < finalDistWellCleaved.getMaxScore(); score ++){
				out.println(score+"\t"+finalDistWellCleaved.getProbability(score));
			}
		}
	}
}
