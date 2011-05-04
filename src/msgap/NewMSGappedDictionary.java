package msgap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import msgf.AminoAcidGraph;
import msgf.GeneratingFunction;
import msgf.NominalMass;
import msgf.ScoredSpectrum;
import msscorer.NewRankScorer;
import msscorer.NewScorerFactory;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;

import parser.MS2SpectrumParser;
import parser.MSGappedDictionaryPSM;
import parser.MSGappedDictionaryParser;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraIterator;
import parser.PSMList;
import suffixtree.actions.ExactMatching;



public class NewMSGappedDictionary {
	static private boolean verbose = true;
	static private int numSpecLimit = Integer.MAX_VALUE;
	static private int minScanNum = -1;
	static boolean test = false; // TODO make it false
	
	/*generate reversed decoy and return if decoy db is generated or not. 
	 * */
	/*static private boolean generateReversedDecoyDataBase(String decoyDBFileName, String targetDBFileName) throws IOException{
		if(new File(decoyDBFileName).exists()) return false;
		if(verbose) System.out.println("Generating Decoy DataBase from " + targetDBFileName);
		BufferedLineReader in = new BufferedLineReader(targetDBFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(decoyDBFileName));
		int dbSize = 0;
		
		String s, seqTitle = "", seq = "";
		char[] reversedSeq = null;
		
		while((s=in.readLine()) != null){
			if(s.startsWith(">")){
				if(!seqTitle.isEmpty()){
					reversedSeq = new char[seq.length()];
					dbSize += seq.length();
					for(int i=0; i<reversedSeq.length; i++)
						reversedSeq[i] = seq.charAt(seq.length() - i - 1);
					out.write(seqTitle + "_Reversed\n" + new String(reversedSeq) + "\n");
				}
				seqTitle = s;
				seq = "";
			}else{
				seq += s;
			}
		}
		
		reversedSeq = new char[seq.length()];
		dbSize += seq.length();
		for(int i=0; i<reversedSeq.length; i++)
			reversedSeq[i] = seq.charAt(seq.length() - i - 1);
		out.write(seqTitle + "_Reversed\n" + new String(reversedSeq) + "\n");
		
		out.close();
		in.close();
		
		return true;
	}
	*/
	
	static public Iterator<Spectrum> getSpectralIterator(String spectrumFileName){
		Iterator<Spectrum> iterator = null;
		try {
			if(spectrumFileName.endsWith("mgf"))
				iterator = new SpectraIterator(spectrumFileName, new MgfSpectrumParser());
			else if (spectrumFileName.endsWith("mzXML")) 
				iterator = new MzXMLSpectraIterator(spectrumFileName);
			else if (spectrumFileName.endsWith("ms2")) 
			  iterator = new SpectraIterator(spectrumFileName, new MS2SpectrumParser());
		}catch (IOException e) {
			 System.err.println("IOException no spectrum file found named " + spectrumFileName);
			 e.printStackTrace();
			 System.exit(-1);
		}
		return iterator;
	}
	
	
	static private boolean isQualifiedAfterPMCorrection(Spectrum s, Parameters par, boolean correct){
	    boolean isQualified = false;
		int maxThresholdScore = Integer.MIN_VALUE;
	    
	    HashSet<Float> parentMasses = new HashSet<Float>();
	    
	    if(par.correctPM){
	    	s.correctParentMass();
	    	correct = false;
	    }
	    
	    float originalParentMass = s.getParentMass();

	    float correctedParentMass = originalParentMass;
	    
	    NewRankScorer scorer = NewScorerFactory.get(par.getActivationMethod(s), par.enzyme());
	    if(par.getMSGFParaFile()!= null) scorer = new NewRankScorer(par.getMSGFParaFile());
	    	
	    if(correct){
	    	ArrayList<Integer> nominalPMs = new ArrayList<Integer>();
	        for(float pmOffset = 0; pmOffset<=par.pmTolerance().getToleranceAsDa(s.getParentMass()); pmOffset+=0.5f) {
	        	float offsetpm = pmOffset + originalParentMass;
	        	int nm = NominalMass.toNominalMass(offsetpm);
	        	if(!nominalPMs.contains(nm)){
	        		nominalPMs.add(nm);
	        		parentMasses.add(offsetpm);
	        	}
	        	
	        	if(pmOffset > 0) {
	        		offsetpm = originalParentMass - pmOffset;
	        		nm = NominalMass.toNominalMass(offsetpm);
	        		if(!nominalPMs.contains(nm)){
		        		nominalPMs.add(nm);
		        		parentMasses.add(offsetpm);
	        		}
        		}
	        }
	    }else{
	    	parentMasses.add(originalParentMass);
	   }
	   // for(float p : parentMasses){
	   // 	System.out.println(originalParentMass + " " + p + " " +  NominalMass.toNominalMass(p));
	   // }
	 
        for(float pm : parentMasses){
	      GeneratingFunction<NominalMass> gf = null;
	      
	      Spectrum s2 = s.getCloneWithoutPeakList();
	        
	      for(Peak p : s){
        	s2.add(p.clone());
	      }
	      
	      s2.correctParentMass(pm);
	      
	      ScoredSpectrum<NominalMass> ss =  scorer.getScoredSpectrum(s2);
		  AminoAcidGraph g =  new AminoAcidGraph(s2.getParentMass(), par.aaSet(), par.enzyme());
	      
	      if(par.allowNonEnzymaticCleavage()){
	        g.allowNonEnzymaticCleavage();
	      }
	      
	      gf = new GeneratingFunction<NominalMass>(ss, g).enzyme(par.enzyme());
	      gf.doNotBacktrack().doNotCalcProb().doNotCalcNumber().computeGeneratingFunction();
	  		
	      if(!gf.isGFComputed() || gf.getMaxScore() <= Math.min(par.msgfScore(), par.matchScore())) continue;
	     // if(Float.isInfinite(gf.getDictionarySize(par.specProbThresholdBeforeConsideringFlankingAAs()) )) continue;
	      
	      int tmpThreshold = gf.getMaxScore();
	      maxThresholdScore = Math.max(maxThresholdScore, tmpThreshold);
	  		
	    //  if(maxThresholdScore < par.matchScore()) continue;
	  		 
	      if(maxThresholdScore == tmpThreshold){
	        correctedParentMass = pm;
	        isQualified = true;
		    	
	      }	 
	    }
	    
		s.correctParentMass(correctedParentMass);	  
		
		return isQualified;
		
	}
	
	
	//returns the number of spectra
	static private int writeGappedSpectralDictionary(String spectrumFileName, PrintWriter grcFile, PrintWriter scoringParaFile, int cumCount, Parameters par){
	
		int qualifiedSpecNum = 0, correctNum = 0, scanNum = 0;
		float averageDictionarySize = 0;
		
		Iterator<Spectrum> iterator = getSpectralIterator(spectrumFileName);
		
		while(iterator.hasNext())
		{
			Spectrum spec = iterator.next();
			cumCount++; // should be incremented everytime
			scanNum = spec.getScanNum();
			
			if(cumCount > numSpecLimit) break;
			if(minScanNum > scanNum) continue;
			
			if(spec.getCharge() == 0) spec.setPrecursorCharge(2); // charge correction
			else if(par.minSpecCharge() > spec.getCharge() || par.maxSpecCharge() < spec.getCharge())
		    	continue;
			
			float originalParentMass = spec.getParentMass();

			//spec.correctParentMass();
			//originalParentMass = spec.getParentMass();
			
			boolean isCorrect = false;
			boolean isQualified = isQualifiedAfterPMCorrection(spec, par, true);
			
			/*if(test) {
				if(spec.getAnnotation().getNominalMass() == NominalMass.toNominalMass(spec.getParentMass()))
					qualifiedSpecNum++;
				continue;
			}*/
		    if(!isQualified) continue;		
		    
			qualifiedSpecNum++;
		    //spec.correctParentMass(spec.getParentMass()+1);
		    
		    NewRankScorer scorer = NewScorerFactory.get(par.getActivationMethod(spec), par.enzyme());
		    if(par.getMSGFParaFile()!= null) scorer = new NewRankScorer(par.getMSGFParaFile());
		    
			ScoredSpectrum<NominalMass> scoreSpec =  scorer.getScoredSpectrum(spec);
		    AminoAcidGraph graph =  new AminoAcidGraph(spec.getParentMass(), par.aaSet(), par.enzyme());
		    scoreSpec.precomputeNodeScores(graph.getIntermediateNodeList());
			GappedGeneratingFunction<NominalMass> gap = new GappedGeneratingFunction<NominalMass>(scoreSpec, graph, par);
			//System.out.println(spec.getParentMass() + "\t" + scoreSpec.getScore(new Peptide(spec.getAnnotationStr(), par.aaSet()), graph));

			gap.computeGappedGeneratingFunction();
					
			ScoringParameter scoringPara = new ScoringParameter(cumCount, spec, spectrumFileName, gap, originalParentMass);
			
			scoringPara.outputFIle(scoringParaFile);
			
			//int delta = Math.min((int)((spec.getParentMass()-Composition.H2O)/121.6), par.delta());
			
			gap.setBackTrackPara(true, par.delta(), par.dictionarySize(), par.maxGapMass());
			String toWriteInGRCFile = "";
			
			int scoreThreshold = gap.getThresholdScore(par.specProbThresholdBeforeConsideringFlankingAAs());
			while(gap.getNumEqualOrBetterPeptides(scoreThreshold) == Float.POSITIVE_INFINITY){
				scoreThreshold++;
			}
			
			for(GappedReconstruction<NominalMass> gappedReconstruction : gap.generateGappedDictionary(Math.max(scoreThreshold, par.matchScore()))){
				toWriteInGRCFile += gappedReconstruction.getGapMassRepresentation() + "\n";
				averageDictionarySize++;
				
				if(spec.getAnnotation() != null)
					if(spec.getAnnotation().isGappedPeptideTrue(gappedReconstruction.getPrefixMassRepresentation())) isCorrect = true;
	
			}
			
			String headerLine = String.format("#%d\t%d\t%s\t%.3f\t%d\t%s\t%s\n", cumCount, spec.getScanNum(), spectrumFileName, spec.getPrecursorPeak().getMz(), spec.getCharge(), spec.getActivationMethod(), spec.getAnnotationStr());
			if(!toWriteInGRCFile.isEmpty()){
			  // the header contains extra information about the spectrum
				grcFile.printf(headerLine);
				grcFile.print(toWriteInGRCFile);
			}
			
			if(isCorrect) correctNum++;
			
			if(verbose){
				//spec.correctParentMass();
				System.out.print(headerLine);	
				//if(spec.getAnnotation() != null) System.out.println(isCorrect + "\t" + spec.getParentMass() + "\t" + spec.getAnnotation().getParentMass() + "\t" + originalParentMass + "\t" + gap.getSpectralProbability(spec.getAnnotation()));
			}
			
			gap = null;
			grcFile.flush();
  			scoringParaFile.flush();
		}
		iterator = null;
		
		if(test) {
			System.out.println("\t# Total Spectrum :"+cumCount+"\t# Qualified Spectrum :"+qualifiedSpecNum);
			if(correctNum>0) System.out.println("\t# Correct Spectrum :"+correctNum);
			System.out.println("\tAverage Dictionary Size : " + (averageDictionarySize/qualifiedSpecNum));
		}

		return cumCount;
	}
	
	static private void writeSpectrumMatchResult(Parameters par){
		
		ExactMatching.run(par);
		
		//PSMList<MSGappedDictionaryPSM> list = MSGappedDictionaryParser.parse(par.getOutFileName()+".txt", par.aaSet());
		//System.out.println("Distinctive pep number : " + list.getDistinctivePeptideSet().size());
		//System.out.println("Matched spectra number : " + list.getDistinctiveSpectralSet().size());
		
	}
	
	
	
	
	static void initialize(Parameters par){
		
		ModifiedAAinGap.initialize(par.aaSet(), par.maxModNum(), par.maxModNumPerOneModification());

	/*	if(par.decoydbFileName() != null){
			try {
				generateReversedDecoyDataBase(par.decoydbFileName(), par.dbFileName());
			} catch (IOException e) {
				 System.err.println("IOException caught when running decoy db()");
				 e.printStackTrace();
				 System.exit(-1);
			}
		}*/

	}
	
	
	
	public static void run(Parameters par, boolean ver){
		verbose = ver;
		long starttime = System.currentTimeMillis();
		
		initialize(par);		
		
		String outFileName = par.getOutFileName();
		//String decoyOutFileName = outFileName+".decoy_out";
		
		if(new File(outFileName).exists()) {
		  new File(outFileName).delete();
		}
	//	if(new File(decoyOutFileName).exists()) {
	//	  new File(decoyOutFileName).delete();
	//	}
	
		String grcFileName = par.getGRCPath();
		String scoringParaFileName = par.getSPRPath();
		
		if (ver) {
		  System.out.println("GRC path " + grcFileName);
		  System.out.println("ScoreParams path " + scoringParaFileName);
		}
		
		//System.out.println(scoringParaFileName);
		/*
		if(!par.useGeneratedGrcFile()) {
		  if (new File(grcFileName).exists()) {
			  new File(grcFileName).delete();
		  }
		  if (new File(scoringParaFileName).exists()) {
			  new File(scoringParaFileName).delete();
		  }
		}
		 */
		
		int cumCount = 0;

		try {
  		PrintWriter grcFile = new PrintWriter(new FileWriter(grcFileName, par.useGeneratedGrcFile()));
  		PrintWriter scoringParaFile = new PrintWriter(new FileWriter(scoringParaFileName, par.useGeneratedGrcFile()));
  		
  		for(String specFile : par.specFiles()){
  			if(verbose) System.out.println("** Processing "+specFile);
  			
  			if(!par.useGeneratedGrcFile())
  				cumCount = writeGappedSpectralDictionary(specFile, grcFile, scoringParaFile, cumCount, par);
  			else
  				System.out.println("[" + grcFileName + " ] is used for [" + specFile +"]");
  		}
  		grcFile.close();
  		scoringParaFile.close();
		}
		catch (IOException ioe) {
		  System.err.println("Error creating the grc & spr file");
		  System.err.println(ioe);
		  System.exit(-1);
		}
		
		if(par.dbFileName() != null){
			writeSpectrumMatchResult(par);
		}
		

		if(verbose){ // report!
			long elapsedTime = System.currentTimeMillis() - starttime;
			if(cumCount > 0) System.out.println("\nTotal number of specta: " + cumCount);
			System.out.println("Totoal amont of time elapsed: " + (float)elapsedTime/1000 + " sec");
			if(cumCount > 0) System.out.println("Per spectrum: " + (float)elapsedTime / cumCount / 1000 + " sec");
		}
	}
	
	private static void GenerateSampleConfigFile(String filename){
		String s = "";
		s+="# expression after '#' character is ignored.\n";
		s+= "-s /home/user/spectrum.mzXML # input spectrum file\n";
		s+= "-o /home/user/output # output prefix\n";
		s+= "#-d /home/user/protein.fasta # db file\n";
		s+= "#-t 2.5Da # parent mass tolerance\n";
		s+= "#-c 2:3 # charge range\n";
		s+= "#-e 3 # enzyme selection 0: No enzyme, 1: Trypsin (default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N\n";
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(filename));
			out.write(s);
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	
	}
	
	public static void main(String[] argv){
		if (argv.length < 1) printUsageAndExit();
		
		if(argv[0].equals("-h")){
			if(argv.length < 2){
				System.err.print("Please specify the path for sample config file (e.g., -h /home/user/sample.txt).");
				System.exit(0);
			}
			
			GenerateSampleConfigFile(argv[1]);
			System.exit(0);
		}
		Parameters par;
		
		if(argv.length == 1)
			par = new Parameters(argv[0], true);
		else 
			par = new Parameters(argv, 0);

	    //if(argv.length > 1) numSpecLimit = Integer.parseInt(argv[1]);
	    //if(argv.length > 2) minScanNum = Integer.parseInt(argv[2]);
	    
		run(par, true);
	}
	
	
	  private static void printUsageAndExit() {
		    System.out.println("MS-Java: MS-GappedDictionary ver. 102010\n");
		    System.out.println("Usage 1: java -jar -Xmx3500m MS-GappedDictionary.jar [configFilePath]");
		   // System.out.println("	or");
		    System.out.print("Usage 2: java -jar -Xmx3500m MS-GappedDictionary.jar\n"
					+ "\t-s spectrumFile (*.mzXML, *.mgf, *.ms2, or directory name containing spectrum files)\n" //, *.mgf, *.pkl, *.ms2)\n"
					+ "\t-o outputFilePrefix ([prefix].grc and [prefix].spr will be generated. If DB file is specified, [prefix].txt will be generated.)\n"
					+ "\t[-d DBFile (*.fasta or directory name containing *.fasta files)]\n"
					
					+ "\t[-t parentMassTolerance (ex: 2.5Da, 50ppm, no space allowed Default: 2.0Da)]\n"
					+ "\t[-c chargeRange] (ex: -c 2, -c 2:4 (2 to 4), -c :4 (1 to 4), -c 2: (from 2 to inf) Default: all charges. Charge 0 spectra will be assumed to be charge 2)]\n"
					+ "\t[-e enzymeID] (0: No enzyme, 1: Trypsin (default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N)\n"
					+ "\t[-m fragmentationMethodID] (1: CID (default) , 2: ETD)\n"//, 3: CID/ETD pair)\n"
					+ "\t[-fixMod 0/1/2] (0: NoCysteineProtection, 1: Carbamidomethyl-C (default), 2: Carboxymethyl-C)\n"
					+ "\t[-p spectralProbabilityThreshold] (Default: 1e-9)]\n"
					+ "\t[-filter msgfScoreThreshold (default: 0)]\n"
					+ "\t[-ps peptideMatchScoreThreshold (default: 0)]\n"
					+ "\t[-l minimumGappedPeptideLength] (Default: 5)]\n"
					+ "\t[-u] (use previously built output file (*.grc, *.spr) for further DB search.)\n"//, 3: CID/ETD pair)\n"
		    );
		    System.exit(0);
	  }
	  
}
