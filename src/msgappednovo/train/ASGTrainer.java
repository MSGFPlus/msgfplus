package msgappednovo.train;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import msgappednovo.IPE;
import msgappednovo.AugmentedSpectrumGraph;
import msgappednovo.AugmentedSpectrumGraph.Edge;
import msgappednovo.AugmentedSpectrumGraph.Node;
import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.MgfSpectrumParser;

public class ASGTrainer {
	
	final int SPECNUM = 50000;
	
	private boolean highPrecision = false;
	private String specfilename;
	private String para;
	private static AminoAcidSet aaSet;
	private WindowFilter filter;
	private Tolerance tol, pmtol;
	
	private float[] ionOffset;
	private float[][] ionWeights;
	private float [][][] covariances;
	private float [][][] edgeAccuracies;
	
	private void trainNodeCovariance(int specCharge, int maxRank, int maxIterationNum){
		Iterator<Spectrum> iterator;
		
		int length = 20; // interpolation needed
		
		int sn = 0;
		HashMap<Node, Integer> minAANumTable = null;
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			IPE gen = new IPE(para, tol, pmtol, aaSet).filter(filter);
			int maxAANum = 5;
			float[][][] divider = new float[maxAANum+1][length+1][length+1];
			float[][][] meanL = new  float[maxAANum+1][length+1][length+1];
			float[][][] meanR = new  float[maxAANum+1][length+1][length+1];
			
			while(iterator.hasNext()){ // mean estimation
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				sn++;
				if(sn > SPECNUM) break;
				
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training covarances (mean estimation) : " + sn);
				}
				
				AugmentedSpectrumGraph mag = new AugmentedSpectrumGraph(spec, gen, para, maxRank);
				minAANumTable = AugmentedSpectrumGraph.getMinAANumTable();

				ArrayList<Node> nodes = mag.getNodes();
				ArrayList<Integer> correctNodeNumbers = getCorrectNodeBinNumbers(spec);
				
				for(int i=1; i<nodes.size()-1;i++){ // exclude sink
					Node rnode = nodes.get(i);
					int i3 = Math.round((meanR[0][0].length-1)*rnode.getAccuracy());
					for(Edge e : mag.getEdges(rnode)){
						Node lnode = e.getLeftNode();
			
						int i1 = Math.min(meanR.length-1, minAANumTable.get(Node.getNode(rnode.getBinNumber()-lnode.getBinNumber())));
						int i2 = Math.round((meanR[0].length-1)*lnode.getAccuracy()); 

						divider[i1][i2][i3] ++;
						if(correctNodeNumbers.contains(rnode.getBinNumber())) meanR[i1][i2][i3] ++;
						if(correctNodeNumbers.contains(lnode.getBinNumber())) meanL[i1][i2][i3] ++;
					}
				}
			}
			for(int i=0;i<meanR.length;i++){
				for(int j=0;j<meanR[i].length;j++){
					for(int k=0;k<meanR[i][j].length;k++){
						if(divider[i][j][k] > 0){
							meanR[i][j][k] /= divider[i][j][k];
							meanL[i][j][k] /= divider[i][j][k];	
						}
					}
				}
			}		
			
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			prevsn = 0; sn = 0; 
			
			covariances = new float[maxAANum+1][length+1][length+1];
			divider = new float[maxAANum+1][length+1][length+1];
			
			while(iterator.hasNext()){ // cov estimation
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				sn++;
				if(sn > SPECNUM) break;
				
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training covarances (cov estimation) : " + sn);
				}
				
				AugmentedSpectrumGraph mag = new AugmentedSpectrumGraph(spec, gen, para, maxRank);
				minAANumTable = AugmentedSpectrumGraph.getMinAANumTable();

				ArrayList<Node> nodes = mag.getNodes();
				ArrayList<Integer> correctNodeNumbers = getCorrectNodeBinNumbers(spec);
				
				for(int i=1; i<nodes.size()-1;i++){ // exclude sink
					Node rnode = nodes.get(i);
					int i3 = Math.round((meanR[0][0].length-1)*rnode.getAccuracy());
					for(Edge e : mag.getEdges(rnode)){
						Node lnode = e.getLeftNode();
			
						int i1 = Math.min(meanR.length-1, minAANumTable.get(Node.getNode(rnode.getBinNumber()-lnode.getBinNumber())));
						int i2 = Math.round((meanL[0].length-1)*lnode.getAccuracy()); 

						divider[i1][i2][i3] ++;
						
						float lc = 0, rc = 0;
						if(correctNodeNumbers.contains(rnode.getBinNumber())) rc = 1;
						if(correctNodeNumbers.contains(lnode.getBinNumber())) lc = 1;
						covariances[i1][i2][i3] += (rc - meanR[i1][i2][i3]) * (lc - meanL[i1][i2][i3]);
					}
				}
			}
			for(int i=0;i<meanR.length;i++){
				for(int j=0;j<meanR[i].length;j++){
					for(int k=0;k<meanR[i][j].length;k++){
						if(divider[i][j][k] > 0)
							covariances[i][j][k] /= divider[i][j][k];
					}
				}
			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
	}
	
	private void trainEdgeAccuracy(int specCharge, int maxRank, int maxIterationNum){
		Iterator<Spectrum> iterator;
		
		int length = 20; // interpolation needed
		
		int sn = 0;
		HashMap<Node, Integer> minAANumTable = null;
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			IPE gen = new IPE(para, tol, pmtol, aaSet).filter(filter);
			int maxAANum = 5;
			float[][][] divider = new float[maxAANum+1][length+1][length+1];
			edgeAccuracies = new float[maxAANum+1][length+1][length+1];
			
			while(iterator.hasNext()){ 
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				sn++;
				if(sn > SPECNUM) break;
				
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training edge accuracy : " + sn);
				}
				
				AugmentedSpectrumGraph mag = new AugmentedSpectrumGraph(spec, gen, para, maxRank);
				minAANumTable = AugmentedSpectrumGraph.getMinAANumTable();

				ArrayList<Node> nodes = mag.getNodes();
				ArrayList<Integer> correctNodeNumbers = getCorrectNodeBinNumbers(spec);
			
				for(int i=1; i<nodes.size()-1;i++){ // exclude sink
					Node rnode = nodes.get(i);
	
					for(Edge e : mag.getEdges(rnode)){
						Node lnode = e.getLeftNode();
			
						if(!correctNodeNumbers.contains(lnode.getBinNumber())) continue;
						
						int i1 = Math.min(divider.length-1, minAANumTable.get(Node.getNode(rnode.getBinNumber()-lnode.getBinNumber())));
						int i2 = Math.round((divider[0].length-1)*rnode.getAccuracy()); 
						int i3 = Math.round((divider[0][0].length-1)*lnode.getAccuracy()); 
						
						divider[i1][i2][i3] ++;
						if(correctNodeNumbers.contains(rnode.getBinNumber())) edgeAccuracies[i1][i2][i3] ++;
					}
				}
			}
			for(int i=0;i<edgeAccuracies.length;i++){
				for(int j=0;j<edgeAccuracies[i].length;j++){
					for(int k=0;k<edgeAccuracies[i][j].length;k++){
						if(divider[i][j][k] > 0){
							edgeAccuracies[i][j][k] /= divider[i][j][k];
						}
					}
				}
			}		
			
		
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
	}
	
	private void trainIonWeights(int specCharge, int maxRank, int maxIterationNum){
		Iterator<Spectrum> iterator;
		double[][][] meanY;
		double[][][] meanYY;
		double[][][] meanXY;
		double[] meanX;
		double[] divider;
		
		int sn = 0;
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			IPE gen = new IPE(para, tol, pmtol, aaSet).filter(filter);
			ArrayList<IonType> sigIons = gen.getSigIonsOrderedByIntensityWithOutNoiseIon(specCharge);

			meanY = new double[1<<(sigIons.size())][][];
			meanYY = new double[1<<(sigIons.size())][][];
			meanXY = new double[1<<(sigIons.size())][][];
			meanX = new double[1<<(sigIons.size())];
			divider = new double[1<<(sigIons.size())];
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				sn++;
				if(sn > SPECNUM) break;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training ion weights for NewMassAccuracyGraph : " + sn);
				}
				
				ArrayList<Integer> correctNodeNumbers = getCorrectNodeBinNumbers(spec);
				
				spec.setRanksOfPeaks();
		
				HashMap<Node, float[]> nodeIonProbs = AugmentedSpectrumGraph.getNodeIonProbs(spec, gen, maxRank);
				
				/*
				for(int i=0; i<500; i++) if(nodeIonProbs.containsKey(NewNode.getNewNode(i)))
					if(nodeIonProbs.get(NewNode.getNewNode(i))[0]>0)
						System.out.println(nodeIonProbs.get(NewNode.getNewNode(i))[0]);
				*/
				
				for(Node node : nodeIonProbs.keySet()){ // node
					int caseIndex=0;
					int cardinality = 0;
					float[] ionprobs = nodeIonProbs.get(node);
					
					for(int j=0; j<sigIons.size(); j++){ // ion
						if(ionprobs[j] > 0){
							caseIndex += 1<<j;
							cardinality++;
						}
					}
					if(cardinality==0) continue;
					int[] indices = new int[cardinality];
					int t = 0;
					for(int j=0; j<sigIons.size(); j++){ // ion
						if(ionprobs[j] > 0){
							indices[t++] = j;
						}
					}
					
					int d=0;
					if(correctNodeNumbers.contains(node.getBinNumber())){
						d=1;
					//	System.out.print(node.getBinNumber());
					//	for(float p : ionprobs) System.out.print("\t" + p);
					//	System.out.println();
					}
					divider[caseIndex]++;
					meanX[caseIndex] += d;
					
					if(meanY[caseIndex] == null){
						meanY[caseIndex] = new double[1][cardinality];
						meanYY[caseIndex] = new double[cardinality][cardinality];
						meanXY[caseIndex]  = new double[1][cardinality];
					}
					
					for(int j =0; j< cardinality;j++){ // ion
						float prob = ionprobs[indices[j]];
						meanY[caseIndex][0][j] += 	prob;
						meanXY[caseIndex][0][j] += d * prob;
						for(int k=0; k< cardinality; k++){
							float prob2 = ionprobs[indices[k]];
							meanYY[caseIndex][k][j] += prob2 * prob;
						}
					}
				}
			}
			
			ionOffset= new float[1<<(sigIons.size())];
			ionWeights = new float[1<<(sigIons.size())][];
			
			for(int caseIndex = 0; caseIndex < 1<<(sigIons.size()); caseIndex++){
				if(divider[caseIndex] < 100) continue;
				
				meanX[caseIndex]/=divider[caseIndex];
				
				//System.out.println(meanX[caseIndex]);
				
				for(int i=0; i<meanY[caseIndex][0].length; i++){
					meanY[caseIndex][0][i]/=divider[caseIndex];
				}
				for(int i=0; i<meanXY[caseIndex][0].length; i++) meanXY[caseIndex][0][i] = 
					(meanXY[caseIndex][0][i]-divider[caseIndex]*meanX[caseIndex]*meanY[caseIndex][0][i])/(divider[caseIndex]-1);
				for(int i=0; i<meanYY[caseIndex].length; i++){
					for(int j=0; j<meanYY[caseIndex][i].length; j++) meanYY[caseIndex][i][j] = (meanYY[caseIndex][i][j] - divider[caseIndex]*meanY[caseIndex][0][i]*meanY[caseIndex][0][j])/(divider[caseIndex]-1);
				}
				
				meanYY[caseIndex] = MatrixCalculus.sum(meanYY[caseIndex], MatrixCalculus.multiply(MatrixCalculus.transpose(meanY[caseIndex]), meanY[caseIndex]), -1);
				meanYY[caseIndex] = MatrixCalculus.invert(meanYY[caseIndex]);
	
				double[][] coeff  = MatrixCalculus.multiply(meanYY[caseIndex], MatrixCalculus.transpose(MatrixCalculus.sum(meanXY[caseIndex], meanY[caseIndex], -meanX[caseIndex])));

				ionOffset[caseIndex] = (float) (meanX[caseIndex]-MatrixCalculus.multiply(meanY[caseIndex], coeff)[0][0]);
				ionWeights[caseIndex] = new float[coeff.length];
				for(int i=0; i<ionWeights[caseIndex].length;i++){
					ionWeights[caseIndex][i] = (float) coeff[i][0];
				}
			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
	}
	
	private ArrayList<Integer> getCorrectNodeBinNumbers(Spectrum spec){
		ArrayList<Integer> correctNodeNumbers = new ArrayList<Integer>();
		correctNodeNumbers.add(0);
		for(float m : spec.getAnnotation().getPRMMasses(true, 0)){
			if(highPrecision){ // TODO 		+ sink node!			
				float delta = tol.getToleranceAsDa(Math.max(m, spec.getParentMass() - m));
				int min = Node.getNode(m-delta).getBinNumber();
				int max = Node.getNode(m+delta).getBinNumber();
				
				for(int k=min;k<=max;k++)
					correctNodeNumbers.add(k);
			}else{
				correctNodeNumbers.add(Node.getNode(m).getBinNumber());
			}
		}
		
		if(highPrecision){
			
		}else{
			correctNodeNumbers.add(Node.getNode(spec.getAnnotation().getMass()).getBinNumber());
		}
		
		return correctNodeNumbers;
	}
	
	protected void train(int charge, int maxRank, int maxIterationNum){
		
		IPE.setMaxIterationNum(0);

		//trainIonWeights(charge, maxRank, maxIterationNum);
		//write(para, charge);
		
		IPE.setMaxIterationNum(0);
		trainNodeCovariance(charge, maxRank, maxIterationNum);
		write(para, charge);
		
		IPE.setMaxIterationNum(0);
		trainEdgeAccuracy(charge, maxRank, maxIterationNum);
		write(para, charge);
	
	}
	
	private void write(String para, int charge){
		PrintWriter out;
		try {
			out = new PrintWriter(new FileWriter(para, true));
			
			if(ionWeights != null){
				int written = 0;
				
				
				out.println("#IONWEIGHT\t" + charge);
				for(int i=0; i< ionWeights.length;i++){
					boolean towrite = false;
					if(ionWeights[i] == null) continue;
					for(int j=0; j< ionWeights[i].length; j++){
						if(ionWeights[i][j] != 0) towrite = true;
					}
					if(!towrite) continue;
					if(ionOffset[i] > 0.3 || ionOffset[i] < -0.3) continue;
					written++;
					
					out.println("#OFF\t"+i + "\t" + ionOffset.length);
					out.println(ionOffset[i]);
					
					out.println("#WEIGHTS\t"+ i + "\t" +ionWeights[i].length+"\t" + ionWeights.length);
					for(float w : ionWeights[i]) out.println(w);
				}
				
			
				out.println("#ENDIONWEIGHT\t" + written);
			}
			
			if(covariances != null){
				
				out.println("#COVARIANCE\t"+covariances.length+"\t"+covariances[0].length+"\t"+covariances[0][0].length+"\t"+  charge);
				
				for(int i=0;i<covariances.length;i++){
					boolean toWrite = false;
					for(int j=0;j<covariances[i].length;j++){
						for(int k=0;k<covariances[i][j].length;k++){
							if(covariances[i][j][k] != 0){
								toWrite = true;
								break;
							}
						}
						if(toWrite == true) break;
					}
					
					if(toWrite){
						out.println("##NUM\t"+i);
						for(int j=0;j<covariances[i].length;j++){
							for(int k=0;k<covariances[i][j].length;k++){
								out.print(covariances[i][j][k]+" ");
							}
							out.println();
						}
					}
				}
				
				out.println("#ENDCOVARIANCE");
			}
			
			if(edgeAccuracies != null){
				
				out.println("#EDGEACCURACY\t"+edgeAccuracies.length+"\t"+edgeAccuracies[0].length+"\t"+edgeAccuracies[0][0].length+"\t"+  charge);
				
				for(int i=0;i<edgeAccuracies.length;i++){
					boolean toWrite = false;
					for(int j=0;j<edgeAccuracies[i].length;j++){
						for(int k=0;k<edgeAccuracies[i][j].length;k++){
							if(edgeAccuracies[i][j][k] != 0){
								toWrite = true;
								break;
							}
						}
						if(toWrite == true) break;
					}
					
					if(toWrite){
						out.println("##NUM\t"+i);
						for(int j=0;j<edgeAccuracies[i].length;j++){
							for(int k=0;k<edgeAccuracies[i][j].length;k++){
								out.print(edgeAccuracies[i][j][k]+" ");
							}
							out.println();
						}
					}
				}
				
				out.println("#ENDEDGEACCURACY");
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	protected ASGTrainer(String specfilename, String para, Tolerance tol, Tolerance pmtol, AminoAcidSet aaSet){
		this.specfilename = specfilename;
		this.para = para;
		this.tol = tol;
		this.pmtol = pmtol;
		ASGTrainer.aaSet = aaSet;
	}
	
	protected ASGTrainer filter(WindowFilter b) {filter = b; return this;}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String trainmgf = "/home/kwj/workspace/inputs/Training/shewLengthAll.mgf"; //CID_Tryp_Confident.mgf Zubarev_HCD_Annotated.mgf shewLengthAll.mgf
		int charge = 2;
		
		//Tolerance tol = new Tolerance(10f, true);
	//	Tolerance pmtol = new Tolerance(20f, true);
		Tolerance tol = new Tolerance(0.5f, false);
		Tolerance pmtol =  new Tolerance(0.5f, false);
		
		if(args.length>0){
			trainmgf = args[0];
			charge = Integer.parseInt(args[1]);
		}
		
		String para = trainmgf.substring(0, trainmgf.lastIndexOf(".")) + ".par";
		aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		for(int c=2;c<=charge;c++){
			ASGTrainer mtrainer = new ASGTrainer(trainmgf, para, tol, pmtol,  aaSet).filter(new WindowFilter(8, 50));
			mtrainer.train(c, 100, 100);
		}
	}

}
