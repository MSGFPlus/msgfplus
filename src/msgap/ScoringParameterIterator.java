package msgap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;


import msgf.AminoAcidGraph;
import msgf.NominalMass;
import msgf.ScoreDist;
import msgf.ScoreDistFactory;
import msgf.ScoredSpectrum;
import msscorer.NewRankScorer;
import msscorer.NewScorerFactory;
import msutil.AminoAcidSet;
import msutil.Spectrum;

public class ScoringParameterIterator implements  Iterator<ScoringParameter>,Iterable<ScoringParameter> {
	private String sprFileName;
	private Parameters par;

  private boolean hasNext = true;
  private BufferedReader sprReader;
  private Iterator<Spectrum> specIterator = null;
  private String currentSpecFile = new String();
  private String nextLine = new String();
    
	public ScoringParameterIterator(Parameters par) { // TODO: first parse add
		this.sprFileName = par.getSPRPath();
		this.par = par;

	    try {
	    	  sprReader = new BufferedReader(new FileReader(sprFileName));
	      }
	      catch (IOException ioe) {
	        System.err.println(ioe);
	        System.exit(-9);
	      }
	      
	}  
	

	    @Override
	    public boolean hasNext() {
	      return hasNext;
	    }

	    @Override
	    public ScoringParameter next() {
	      if (!hasNext) return null;
	      
	      if(nextLine.isEmpty()){// init
	    	  try {
				nextLine = sprReader.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	      }
	      
	      if(!nextLine.startsWith("#")) return null;
	      
    	  nextLine = nextLine.substring(1); // elminimate '#'
    	  String[] token = nextLine.split("\t");
    	  int specID = Integer.parseInt(token[0]);
    	  String specFileName = token[1];
    	  int scanNum = Integer.parseInt(token[2]);
    	  int charge = Integer.parseInt(token[3]);
    	  float originalParentMass = Float.parseFloat(token[4]);
    	  float correctedParentMass = Float.parseFloat(token[5]);
    	  String annotation = token[6];
    	  
    	  Spectrum spec = null;
    	  
    	  if(specIterator == null || !currentSpecFile.equals(specFileName)){
    		  currentSpecFile = specFileName;
    		  specIterator = NewMSGappedDictionary.getSpectralIterator(currentSpecFile);
    	  }
    	  
    	  while(specIterator.hasNext()){
    		  spec = specIterator.next();
    		  if(spec.getScanNum() == scanNum) break;
          // cannot assume that scanNumbers are of increasing order.
    		  //if(spec.getScanNum() > scanNum) {
    		    //return null;
    		  //}
    	  }
	      
    	  if(spec == null) return null;
    	  
    	  spec.setPrecursorCharge(charge);
    	  spec.correctParentMass(correctedParentMass);
    	  
		  NewRankScorer scorer = NewScorerFactory.get(par.getActivationMethod(spec), par.enzyme());
		  if(par.getMSGFParaFile()!= null) scorer = new NewRankScorer(par.getMSGFParaFile());
		  
    	  ScoredSpectrum<NominalMass> scoredSpec =  scorer.getScoredSpectrum(spec);
    	  AminoAcidGraph graph =  new AminoAcidGraph(spec.getParentMass(), par.aaSet(), par.enzyme());
    	  if(par.allowNonEnzymaticCleavage()) graph.allowNonEnzymaticCleavage();
    	  
    	  String s;
    	  ScoreDist finalDist = null;
    	  ScoreDist finalDistWellCleaved = null;
    	  ScoreDistFactory factory = new ScoreDistFactory(false, true);
    	  
	      try {
	    	int mod = 0;
			while((s=sprReader.readLine())!=null){
			  	if(s.startsWith("#")){
			  		nextLine = s;
			  		break;
			  	}
			  	
			  	String[] t = s.split("\t");
		  		
			  	if(s.startsWith("DIST1")){
			  		mod = 0;
			  		finalDist = factory.getInstance(Integer.parseInt(t[1]), Integer.parseInt(t[2]));
			  		continue;
			  	}else if(s.startsWith("DIST2")){
			  		mod = 1;
			  		finalDistWellCleaved = factory.getInstance(Integer.parseInt(t[1]), Integer.parseInt(t[2]));
			  		continue;
			  	}
				
			  	if(mod == 0 && finalDist!=null){
			  		finalDist.setProb(Integer.parseInt(t[0]), Float.parseFloat(t[1]));
			  	}else if(mod == 1 && finalDistWellCleaved!=null){
			  		finalDistWellCleaved.setProb(Integer.parseInt(t[0]), Float.parseFloat(t[1]));
			  	}
			  }
			if(s == null) hasNext = false;
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//System.out.println(scanNum + "\t" + currentSpecFile + "\t" + spec.getParentMass() + "\t" + finalDist.getProbability(0));
		
		return new ScoringParameter(
				specID,
				scoredSpec, graph, 
	  			originalParentMass,
				correctedParentMass, 
	  			scanNum, 
	  			charge,
	  			finalDist, 
	  			finalDistWellCleaved, 
	  			currentSpecFile,
	  			annotation);
 
	    }

	    @Override
	    public void remove() {
	      throw new RuntimeException("Cannot remove entry from the MSGDResultFileParser iterator. This is not supported");
	    }
	    
	  
		@Override
		public Iterator<ScoringParameter> iterator() {
			return this;
		
		}

		


/*	static public void main(String[] argv){
		ScoringParameterIterator tmp = new ScoringParameterIterator(new Parameters("/home/kwj/workspace/confs/conf_standard", true));
		
		while(tmp.hasNext()){
			tmp.next();
		}
	}
*/

}
