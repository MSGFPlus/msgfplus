package msgappednovo.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import msgappednovo.AugmentedSpectrumGraph;
import msgappednovo.IPE;
import msgappednovo.GappedPeptideForMSGappedNovo;
import msgappednovo.MSGappedNovo;
import msgappednovo.AugmentedSpectrumGraph.Edge;
import msgappednovo.AugmentedSpectrumGraph.Node;
import msgappednovo.train.PeakGenerator;
import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.Peak;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.MgfSpectrumParser;

public class AccuracyTest {
	static int MaxSpecNum =10000;
	static float setAccuracyThreshold = 0.8f;
	static float minGappedPeptideAccuracy = 0.0f;
	static int numPerLength = 15, minGappedPeptideLength = 5;
	
	static ArrayList<Peptide> getFragmentedPeptidesFrom(Spectrum spec, ArrayList<IonType> ions, Tolerance tol, int maxRank){
		ArrayList<Peptide> prefixPeptides = new ArrayList<Peptide>();
		
		float[] pms = spec.getAnnotation().getPRMMasses(true, 0);
		for(int i=0; i<pms.length; i++){
			float pm = pms[i];
			boolean isFound = false;
			for(IonType ion : ions){
			//	if(ions.indexOf(ion) > n) continue;
				if(ion instanceof IonType.PrefixIon){
					float mz = ion.getMz(pm);
					for(Peak p : spec.getPeakListByMass(mz, tol)){
						if(p.getRank() <= maxRank){
						//	System.out.println(i+"\t"+p.getRank());
							prefixPeptides.add(new Peptide(spec.getAnnotationStr().substring(0, i+1)));
							isFound = true;
							break;
						}
					}
					break;
				}
				if(isFound) break;
			}
		}
		
		pms = spec.getAnnotation().getPRMMasses(false, 0);
		for(int i=0; i<pms.length; i++){
			float pm = pms[i];
			boolean isFound = false;
			for(IonType ion : ions){
			//	if(ions.indexOf(ion) > n) continue;
				if(ion instanceof IonType.SuffixIon){
					float mz = ion.getMz(pm);
					for(Peak p : spec.getPeakListByMass(mz, tol)){
						if(p.getRank() <= maxRank){
							Peptide sp = new Peptide(spec.getAnnotationStr().substring(0, spec.getAnnotation().size()-1-i));
							if(!prefixPeptides.contains(sp))
									prefixPeptides.add(sp);
							isFound = true;
							break;
						}
					}
					break;
				}
				if(isFound) break;
			}
		}
		
		return prefixPeptides;
	}
	
	
	public static void main(String[] args) throws IOException {
		String inputmgf ="/home/kwj/workspace/inputs/Training/CID_Tryp_Confident.mgf";// trainmgf;//CID_Tryp_Confident
		String trainmgf = "/home/kwj/workspace/inputs/Training/shewLengthAll.mgf";//shewLengthAll.mgf Zubarev_HCD_Annotated.mgf
		String para = trainmgf.substring(0, trainmgf.lastIndexOf(".")) + ".par";
		
		int charge = 2;
		int maxRank = 100;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		for(AminoAcid aa : aaSet){
			System.out.println(aa + " " + aa.getMass());
		}
	//	Tolerance tol = new Tolerance(10f, true);
	//	Tolerance pmtol = new Tolerance(20f, true);
		Tolerance tol = new Tolerance(0.5f, false);
		Tolerance pmtol = new Tolerance(0.5f, false);
		WindowFilter filter = new WindowFilter(8, 50);
		
		float[][] peakaccuracydivider;
		float[][] peakaccuracy;
		float[][] gpaccuracydivider = new float[50][10];
		float[][] gpcorrectaccuracy = new float[50][10];
		int[] numcorrect = new int[50];
		float[] avgsize = new float[50];
		float[] avglength = new float[50];
		float[] avgIdellength = new float[50];
		float[] avgcorrectlength = new float[50];
		int[] sn = new int[50];
		int[] tsn = new int[50];
		int specnum = 0;
		int qualifiedSpecnum = 0;
		int totalnumcorrect = 0;
		float[][] nodeaccuracydivider = new float[50][10];
		float[][] nodeaccuracy = new float[50][10];
		float[][] edgeaccuracydivider = new float[50][10];
		float[][] edgeaccuracy = new float[50][10];
		
		IPE ipe =  new IPE(para, tol, pmtol, aaSet).filter(filter);
		
		peakaccuracydivider = new float[ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(charge).size()][10];
		peakaccuracy = new float[ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(charge).size()][10];
		
		Iterator<Spectrum> iterator = new SpectraIterator(inputmgf, new MgfSpectrumParser());
		
		float time = 0;
		
		HashSet<String> peps = new HashSet<String>();
		
		while(iterator.hasNext()){
			float etime = System.nanoTime(); 
			///////////////////////////////
			
			Spectrum spec = iterator.next();
			if(charge > 0 && spec.getCharge() != charge) continue;
			if(peps.contains(spec.getAnnotationStr())) continue;
			if(spec.getAnnotationStr().contains("C")) continue; // TODO delete
			//if(spec.getScanNum()<312)continue;
			//spec.correctParentMass();
			peps.add(spec.getAnnotationStr());
			int pepLength = spec.getAnnotation().size();
		//	if(pepLength < 13) continue;
			tsn[pepLength]++; specnum++;
			
			
			AugmentedSpectrumGraph graph = new AugmentedSpectrumGraph(spec, ipe, para, maxRank);
			MSGappedNovo mg = new MSGappedNovo(graph);
			
			HashMap<Integer, ArrayList<GappedPeptideForMSGappedNovo>> can =  mg.getCandidateGappedPeptides(numPerLength, minGappedPeptideLength, minGappedPeptideAccuracy);
			ArrayList<GappedPeptideForMSGappedNovo> out = new ArrayList<GappedPeptideForMSGappedNovo>();
			
			float setaccuracy = mg.selectOutputFromCandidates(out, can, numPerLength, minGappedPeptideLength, setAccuracyThreshold);
			
			time += System.nanoTime() - etime; 
			///////////////////////////////
			
			HashMap<Peak, HashMap<IonType, Float>>  profile = ipe.getProfile(spec, maxRank, 100);
			PeakGenerator pg = new PeakGenerator(spec);
			
			for(Peak p : profile.keySet()){
				HashMap<IonType, Float> prof = profile.get(p);
				for(int i=0; i< ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(charge).size(); i++){
					IonType ion = ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(charge).get(i);
					
					int bin = Math.round(prof.get(ion) * peakaccuracydivider[i].length);
					peakaccuracydivider[i][Math.min(bin, peakaccuracydivider[i].length-1)] ++;
					if(pg.isExplainedBy(p, ion, tol, pmtol)){
						peakaccuracy[i][Math.min(bin, peakaccuracy[i].length-1)] ++;
					}
				}
			}
			
			
			ArrayList<Node> correctNodes = new ArrayList<Node>();
			for(int i=0; i<graph.getNodes().size();i++){ 
				Node m = graph.getNodes().get(i);
				for(float prm : spec.getAnnotation().getPRMMasses(true, 0)){
					float delta = tol.getToleranceAsDa(Math.max(prm, spec.getParentMass() - prm));
					if(Math.abs(m.getMass() - prm) <= delta){
						correctNodes.add(m);
					}
				}
			}
			
			for(int i=1; i<graph.getNodes().size()-1;i++){ // except sink source
				
				Node m = graph.getNodes().get(i);
				if(graph.getSinkNode().equals(m)) continue;
				
				float accuracy = m.getAccuracy();
				int bin = Math.round(accuracy * nodeaccuracydivider[pepLength].length);
				
				nodeaccuracydivider[pepLength][Math.min(bin, nodeaccuracydivider[pepLength].length-1)]++;

				int bin2 = 0 ;
				boolean nodeCorrect = correctNodes.contains(m);
				
				for(Edge e : graph.getEdges(m)){
					if(!correctNodes.contains(e.getLeftNode())) continue;
					Integer min = AugmentedSpectrumGraph.getMinAANumTable().get(Node.getGapNode(e));
					if(min != null){
						bin2 =Math.round(e.getAccuracy() * edgeaccuracydivider[pepLength].length);
						edgeaccuracydivider[pepLength][Math.min(bin2, edgeaccuracydivider[pepLength].length-1)]++;
						if(nodeCorrect){
							edgeaccuracy[pepLength][Math.min(bin2, edgeaccuracydivider[pepLength].length-1)]++;
						}
					}
				}
	
				if(nodeCorrect){
					nodeaccuracy[pepLength][Math.min(bin, nodeaccuracydivider[pepLength].length-1)]++;
				}else{
			//		if(bin == nodeaccuracydivider.length) System.out.println("***" + accuracy + "\t" + cgp.toBitSetString() + "\t" + m + "\t" +(spec.getAnnotation().getNominalMass()-m.getNominalMass() ));
				}
				
			}
			
			int maxcorrectgplength = 0;

			for(int len : can.keySet()){
				ArrayList<GappedPeptideForMSGappedNovo> gps = can.get(len);
				for(GappedPeptideForMSGappedNovo gp : gps){
	
					int bin = Math.round(gp.getAccuracy()*gpaccuracydivider[0].length);
					gpaccuracydivider[pepLength][Math.min(bin, gpaccuracydivider[0].length-1)]++;
					if(gp.isCorrect(spec.getAnnotation(), tol, pmtol)){
						maxcorrectgplength = maxcorrectgplength > gp.length()? maxcorrectgplength : gp.length();
						gpcorrectaccuracy[pepLength][Math.min(bin, gpaccuracydivider[0].length-1)]++;
					}
				}
			}

			if(setaccuracy < setAccuracyThreshold)continue;

			
			qualifiedSpecnum++;
			boolean isCorrect = false;
			sn[pepLength]++;
			if(specnum > MaxSpecNum) break;
			System.out.println(spec.getScanNum() + "\t" + spec.getAnnotationStr());
			int size = 0;
			
			ArrayList<Peptide> fragd = getFragmentedPeptidesFrom(spec, ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(charge), tol, maxRank);
			System.out.println(fragd);
			
			for(GappedPeptideForMSGappedNovo gp : out){
				System.out.print(gp + " : " + Math.round(gp.getAccuracy()*100)+"% ");
				if(gp.isCorrect(spec.getAnnotation(), tol, pmtol)){
					isCorrect = true;
					System.out.println(" true*\t" + gp.length() + "\t" + maxcorrectgplength + "\t" + ( fragd.size()+1));
				}else{
					System.out.println(" false\t" + gp.length() +"\t"+  maxcorrectgplength + "\t" + (fragd.size()+1));
				}
				avglength[pepLength] += gp.length();
				size++;
			}
			

			
			if(isCorrect){
				numcorrect[pepLength]++;
				totalnumcorrect++;
			}
			avgsize[pepLength] += size;
			
			avgIdellength[pepLength] += fragd.size() + 1;
			avgcorrectlength[pepLength] += maxcorrectgplength;
				
			System.out.println(size + " " + isCorrect);
		
		}
		
		
		System.out.println();
		for(int i=0;i<sn.length;i++){
			if(sn[i] > 0){
				System.out.println(i+" " +(float) numcorrect[i]/sn[i] + " " +(float) numcorrect[i]/tsn[i]);
			}
		}
		System.out.println();
		for(int i=0;i<sn.length;i++){
			if(sn[i] > 0){
				System.out.println(i+ " " + (float) avgsize[i]/sn[i] + " " +(float) avglength[i]/avgsize[i] + " " + (float)avgcorrectlength[i]/sn[i]+" "+(float) avgIdellength[i]/sn[i]);
			}
		}
		
		System.out.print("%");	
		for(int i=0;i<sn.length;i++){
			if(sn[i] > 0){
				System.out.print(i+" ");
			}
		}
		System.out.println();
		System.out.println("accuracy=[");
		for(int j=0; j<gpcorrectaccuracy[0].length;j++){
			System.out.println("%" + j*gpaccuracydivider[0].length+"-"+ ((j+1)*gpaccuracydivider[0].length) + "% : ");
			for(int i=0;i<sn.length;i++){
				if(sn[i] > 0){
					System.out.print(100 * gpcorrectaccuracy[i][j]/gpaccuracydivider[i][j]+" ");
				}else System.out.print(Float.NaN + " ");
			}
			System.out.println();
		}
		System.out.println("];");
		
		System.out.println();
		System.out.println("nodeaccuracy=[");
		for(int j=0;j<nodeaccuracydivider[0].length;j++){
			System.out.print(j + " ");
			for(int k=0; k<nodeaccuracydivider.length;k++){
				System.out.print(nodeaccuracy[k][j]/nodeaccuracydivider[k][j] + " ");
			}
			System.out.println();
		}
		System.out.println("];");
		
		System.out.println();
		System.out.println("edgeaccuracy=[");
		for(int j=0;j<edgeaccuracydivider[0].length;j++){
			System.out.print(j + " ");
			for(int k=0; k<edgeaccuracydivider.length;k++){
				System.out.print(edgeaccuracy[k][j]/edgeaccuracydivider[k][j] + " ");
			}
			System.out.println();
		}
		
		System.out.println("];");
		
		System.out.println();
		System.out.println("% " + ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(charge));
		System.out.println("peakaccuracy=[");
		for(int j=0;j<peakaccuracydivider[0].length;j++){
			System.out.print(j + " ");
			for(int k=0; k<peakaccuracydivider.length;k++){
				System.out.print(peakaccuracy[k][j]/peakaccuracydivider[k][j] + " ");
			}
			System.out.println();
		}
		System.out.println("];");
		
		
		System.out.println("spec: " + specnum + "\tqspec : " + qualifiedSpecnum);
		System.out.println("accuracy: " + ((float)totalnumcorrect/qualifiedSpecnum));
		
		System.out.println(time/specnum/1e9+" sec/spec ");
	}
}
