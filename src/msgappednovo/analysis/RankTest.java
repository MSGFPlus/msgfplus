package msgappednovo.analysis;

import java.io.IOException;
import java.util.ArrayList;

import unexplainedpeak.GenerateIonDistGraph;
import msgappednovo.IPE;
import msgappednovo.train.Trainer;
import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.WindowFilter;

public class RankTest {
	static public void main(String[] args) throws IOException{
		WindowFilter filter = new WindowFilter(8, 50);
		int charge =2 ; 
		//int maxCharge = 2;
	//	Tolerance tol = new Tolerance(10f, true);
	//	Tolerance pmtol = new Tolerance(20f, true);
		Tolerance tol = new Tolerance(0.5f, false);
		Tolerance pmtol = new Tolerance(0.5f, false);
		
		String trainmgf = "/home/kwj/workspace/inputs/Training/CID_Tryp_Confident.mgf";//shewLengthAll.mgf
		String para = trainmgf.substring(0, trainmgf.lastIndexOf(".")) + ".par";
		String inputmgf ="/home/kwj/workspace/inputs/Training/CID_Tryp_Confident.mgf";// trainmgf;//Zubarev_HCD_Annotated.mgf
		String outmgf = inputmgf.substring(0, inputmgf.lastIndexOf(".")) + "_ranked.mgf";
		String pepoutmgf = "/home/kwj/workspace/inputs/Training/CID_Tryp_Confident_PepNovo+"+charge+".mgf";
		
		int maxRank =100;
		int maxIterationNum = 10;
		
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();

		IPE gen = new IPE(para, tol, pmtol, aaSet).filter(filter);
		gen.run(inputmgf, outmgf, charge, maxRank, maxIterationNum, 0, false);
		
		charge =2;
		ArrayList<IonType> sigIons=gen.getSigIonsOrderedByIntensityWithOutNoiseIon(charge);
		sigIons.remove(IonType.NOISE);
		//System.out.println(sigIons.size());
		float probabiltiy = 0.1f;
	
		//String filename, int rankLimit, ArrayList<IonType> sigIons, int specCharge, float probabiltiy
	//	sigIons =  GenerateIonDistGraph.gen(pepoutmgf, maxRank, sigIons, charge, tol, probabiltiy);
		
		sigIons =  GenerateIonDistGraph.gen(inputmgf, maxRank, sigIons, charge, tol, probabiltiy);

		sigIons=GenerateIonDistGraph.gen(outmgf, maxRank, sigIons, charge, tol, probabiltiy);
		

	}
}
