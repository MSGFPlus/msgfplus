package msgappednovo.train;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import java.util.Iterator;

import parser.MgfSpectrumParser;
import msgappednovo.IPE;
import msgappednovo.parameters.PeakParameter;
import msgappednovo.parameters.SpectrumParameter;
import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;

public class Trainer {
	private String inputmgf;
	private Tolerance tol;
	private Tolerance pmtol;
	private int minCaseNum = 100;
	private int minSpecNum = 1000;
	private int maxConditionNum = 5000;
	private int specMzNum = 5;
	private int peakGroupNum = 10;
	private int peakRatioNum = 10;
	private int peakPartitionNum = 4;
	
	private AminoAcidSet aaSet;
	private boolean discard = true;
	private int maxIonNum = 8;
	private boolean considerIonIntensity = true;
	private WindowFilter filter = null;
	
	public Trainer(String inputmgf, Tolerance tol, Tolerance pmtol, AminoAcidSet aaSet){
		this.inputmgf = inputmgf;
		this.tol = tol;
		this.pmtol = pmtol;
		this.aaSet = aaSet;
	}
	
	public Trainer filter(WindowFilter b) {filter = b; return this;}
	
	public Trainer doNotDiscard() {discard = false; return this;}
	
	public Trainer setMaxIonNum(int n) {maxIonNum = n;  return this;}
	
	private ArrayList<Integer> getChargeRange(){
		Iterator<Spectrum> iterator;
		int[] numPerCharge = new int[100];
		int min = 100, max = 0;
		
		try {
			iterator = new SpectraIterator(inputmgf, new MgfSpectrumParser());
			while(iterator.hasNext()){
				Spectrum s = iterator.next();
				int charge = s.getCharge();
				
				numPerCharge[charge]++;

			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
		
		for(int i = 0 ; i<numPerCharge.length; i++){
			if(numPerCharge[i] > minSpecNum){
				max = max < i ? i  : max;
				min = min > i ? i : min;
			}
		}
		
		ArrayList<Integer> range = new ArrayList<Integer>();
		range.add(min); range.add(max);
		System.out.println("Charge range : " + min + "-" + max);
		return range;
		
		
	}
	
	public void train(String para, int charge, int maxRank, int maxIterationNum){
		
		SpectrumParameter.setSpecMzRangeNum(specMzNum);
		PeakParameter.setGroupNum(peakGroupNum);
		PeakParameter.setPartitionNum(peakPartitionNum);
		PeakParameter.setPeakIntensityRatioNum(peakRatioNum);
		
		new File(para).delete();
		
		ArrayList<Integer> range = getChargeRange();
		
		for(int c=range.get(0); c<=Math.min(range.get(1), charge);c++){
			long time = System.currentTimeMillis();
			String mgf = new String(inputmgf);
			
			IonSelector ionfinder = new IonSelector(mgf, new Tolerance(0.5f, false), filter, maxRank, maxIonNum);
			int sn = ionfinder.train(c, considerIonIntensity);
			
			ArrayList<IonType> sigIons = new ArrayList<IonType>();
			for(IonType ion : ionfinder.getSigIons().keySet()){
				if(!ion.equals(IonType.NOISE)) 
					sigIons.add(ion);
			}
			
			if(sigIons.isEmpty()) continue;//
			
			for(int iterationNum = 0; iterationNum < maxIterationNum; iterationNum++){
				FeatureSelector confinder = new FeatureSelector(mgf, sigIons, tol, pmtol).filter(filter);
				confinder.train(c, maxRank, iterationNum, aaSet);
				
				FeatureProbTrainer contrainer = new FeatureProbTrainer(mgf, ionfinder.getSigIons(), confinder.getSignificantConditions(), tol, pmtol).filter(filter);
				contrainer.train(c, maxRank, iterationNum);
				contrainer.discardConditionsWithLessThan(minCaseNum, maxConditionNum, discard);
				contrainer.writeInFile(para, c, iterationNum);
				
				String nextinputmgf = mgf.substring(0, mgf.length()-1) + iterationNum;
				
				if(iterationNum < maxIterationNum - 1){
					IPE gen = new IPE(para, tol, pmtol, aaSet).filter(filter);
					gen.run(mgf, nextinputmgf, c, maxRank, iterationNum+1, iterationNum, true);
				}
				
				if(iterationNum>0) new File(mgf).delete();
				mgf =  nextinputmgf;
			}
			System.out.println("Charge " + c + " training Done : " + (float)(System.currentTimeMillis() - time)/sn/1000 + " sec/spec");
			new File(mgf).delete();
		}
		
	    /*	for(int c=range.get(0); c<=Math.min(range.get(1), charge);c++){
			MassAccuracyGraphTrainer mtrainer = new MassAccuracyGraphTrainer(inputmgf, para, tol, pmtol,  aaSet,  c);
			mtrainer.train(maxRank, maxIterationNum);
		}*/
		for(int c=range.get(0); c<=Math.min(range.get(1), charge);c++){
			ASGTrainer mtrainer = new ASGTrainer(inputmgf, para, tol, pmtol,  aaSet).filter(filter);
			mtrainer.train(c, maxRank, maxIterationNum);
		}
	}
	
	
	static public void main(String[] args) throws IOException{

		int maxCharge = 2;
		Tolerance tol = new Tolerance(0.5f, false);
		Tolerance pmtol = new Tolerance(0.5f, false);
		
	//	tol = new Tolerance(10f, true);
	//	pmtol = new Tolerance(20f, true);
		
		String trainmgf = "/home/kwj/workspace/inputs/Training/shewLengthAll.mgf";
		
		int maxRank = 100;
		int maxIterationNum = 5;
		
		if(args.length > 0) trainmgf = args[0];
		if(args.length > 1) maxCharge = Integer.parseInt(args[1]);
		if(args.length > 2) maxRank = Integer.parseInt(args[2]);
		if(args.length > 3) maxIterationNum = Integer.parseInt(args[3]);
		
		String para = trainmgf.substring(0, trainmgf.lastIndexOf(".")) + ".par";
		
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		Trainer trainer = new Trainer(trainmgf, tol, pmtol, aaSet).filter(new WindowFilter(8, 50));//.setMaxIonNum(ionnum);
		trainer.train(para, maxCharge, maxRank, maxIterationNum);

	}
}
