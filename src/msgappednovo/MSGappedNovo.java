package msgappednovo;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import parser.MgfSpectrumParser;

import msgappednovo.AugmentedSpectrumGraph.Edge;
import msgappednovo.AugmentedSpectrumGraph.Node;
import msgappednovo.train.MatrixCalculus;
import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;

public class MSGappedNovo {
	private final static float minAccuracyForBackTrack = 0.05f;
	
	private AugmentedSpectrumGraph graph;
	private void union(HashMap<Node, HashMap<Integer, HashMap<Integer, HashSet<Edge>>>> e, Node v, int l, int s, Edge edge){
		if(!e.containsKey(v)) e.put(v, new HashMap<Integer, HashMap<Integer, HashSet<Edge>>>());
		HashMap<Integer, HashMap<Integer, HashSet<Edge>>> k = e.get(v);
		if(!k.containsKey(l)) k.put(l, new HashMap<Integer, HashSet<Edge>>());
		HashMap<Integer, HashSet<Edge>> q = k.get(l);
		if(!q.containsKey(s)) q.put(s, new HashSet<Edge>());
		HashSet<Edge> p = q.get(s);
		p.add(edge);		
	}
	
	public MSGappedNovo(AugmentedSpectrumGraph graph){
		this.graph = graph;
	}
	
	public HashMap<Integer, ArrayList<GappedPeptideForMSGappedNovo>> getCandidateGappedPeptides(int numPerLength, int minGappedPeptideLength, float minGappedPeptideAccuracy){
		HashMap<Integer, ArrayList<GappedPeptideForMSGappedNovo>> candidateGappedPeptides = new HashMap<Integer, ArrayList<GappedPeptideForMSGappedNovo>>();
		HashMap<Node, HashMap<Integer, HashMap<Integer, HashSet<Edge>>>> es = new HashMap<Node, HashMap<Integer, HashMap<Integer, HashSet<Edge>>>>();
		
		union(es, Node.getNode(0f), 0, 0, Edge.getNullEdge());
		
		for(Node v : graph.getNodes()){
			for(Edge e: graph.getEdges(v)){
				// search l nodes 
				Node ln = e.getLeftNode();
				HashMap<Integer, HashMap<Integer, HashSet<Edge>>> subes = es.get(ln);
				if(subes == null) continue;
				for(int l : subes.keySet()){
					HashMap<Integer, HashSet<Edge>> ssubes = subes.get(l);
					for(int s : ssubes.keySet()){
						union(es, v, l+1, s+e.getRoundedLogAccuracy(), e);
					}
				}
			}
		}
		
		for(int l=50;l>=minGappedPeptideLength;l--){
			ArrayList<GappedPeptideForMSGappedNovo> candidates = new ArrayList<GappedPeptideForMSGappedNovo>();
			for(int score = 0; score >= AugmentedSpectrumGraph.getRoundedLogAccuracyFromAccuracy(minAccuracyForBackTrack); score--){
				backtrack(candidates, new ArrayList<Node> (), es, graph.getSinkNode(), l, score, numPerLength);	
			}
			//Collections.sort(candidates, Collections.reverseOrder());
			
			for(int i=0; i<candidates.size(); i++){
				if(candidates.get(i).getAccuracy() < minGappedPeptideAccuracy) candidates.remove(i--);
			}
			
			candidateGappedPeptides.put(l, candidates);
		}
		return candidateGappedPeptides;
	}
	
	public float selectOutputFromCandidates(ArrayList<GappedPeptideForMSGappedNovo> out, HashMap<Integer, ArrayList<GappedPeptideForMSGappedNovo>> candidateGappedPeptides, int numPerLength, int minGappedPeptideLength, float setAccuracyThreshold){
		float setAccuracy = 0;
		
		ArrayList<GappedPeptideForMSGappedNovo> total = new ArrayList<GappedPeptideForMSGappedNovo>();
		
		for(int l : candidateGappedPeptides.keySet()){
			total.addAll(candidateGappedPeptides.get(l));
		}
		
		for(int alpha = 5 ; alpha >=0; alpha -= 1){
			out.clear();
			GappedPeptideForMSGappedNovo.setAlpha(alpha);
			
			Collections.sort(total, Collections.reverseOrder());
			
			for(int i=0; i < Math.min(numPerLength, total.size()); i++){
				out.add(total.get(i));
				setAccuracy = getSetAccuracy(out);
				
				if(setAccuracy > setAccuracyThreshold){// conservative
					GappedPeptideForMSGappedNovo.setAlpha(0);
					Collections.sort(out, Collections.reverseOrder());
					return setAccuracy;
				}
			}
		}
		
		/*
		for(int l = 30; l>=minGappedPeptideLength; l--){
			ArrayList<GappedPeptideForMSGappedNovo> candidates = candidateGappedPeptides.get(l);
			if(candidates == null) continue;
			//total.clear();///////////////
			//total.addAll(candidates);
			total = candidates;
			GappedPeptideForMSGappedNovo.setAlpha(0);
			Collections.sort(total, Collections.reverseOrder());
			
			out.clear();
			
			for(int i=0; i < Math.min(numPerLength, total.size()); i++){
				out.add(total.get(i));
				setAccuracy = getSetAccuracy(out);
				
				if(setAccuracy > setAccuracyThreshold){// conservative
					return setAccuracy;
				}
			}
		}*/
		
		return setAccuracy;
	}
	
	public ArrayList<GappedPeptideForMSGappedNovo> getOutput(int numPerLength, int minGappedPeptideLength, float setAccuracyThreshold, float minGappedPeptideAccuracy){
		ArrayList<GappedPeptideForMSGappedNovo> out = new ArrayList<GappedPeptideForMSGappedNovo>();
		float setAccuracy = selectOutputFromCandidates(out, getCandidateGappedPeptides(numPerLength, minGappedPeptideLength, minGappedPeptideAccuracy), numPerLength, minGappedPeptideLength, setAccuracyThreshold);
		if(setAccuracy >= setAccuracyThreshold) return out;
		else return new ArrayList<GappedPeptideForMSGappedNovo>();
	}
	
	private float getSetAccuracy(ArrayList<GappedPeptideForMSGappedNovo> gappedPeptides){
		float accuracy = 0;
		int i=gappedPeptides.size()-1;
		
		float prob0 = gappedPeptides.get(0).getAccuracy();
		
		if(i>0){
			double[][] P = new double[1][i+1];
			double[][] Q = new double[i+1][i+1];
			
			Random ran = new Random();
			for(int k=0; k<=i; k++){
				GappedPeptideForMSGappedNovo g = gappedPeptides.get(k);
				P[0][k] = g.getAccuracy();
				for(int j=0; j<k; j++){
					
					GappedPeptideForMSGappedNovo u = g.getUnionGappedPeptide(gappedPeptides.get(j));
					float acc =u.getAccuracy();
					
					Q[k][j] = Math.min(g.getAccuracy(), acc);
					Q[k][j] = Math.min(Q[k][j], gappedPeptides.get(j).getAccuracy());
		
					double r = (ran.nextDouble()-0.5) * Double.MIN_VALUE;
					Q[k][j] += r;
					Q[k][j] = Math.max(Q[k][j] , 0);
					Q[k][j] = Math.min(Q[k][j] , 1);
				}
				
				Q[k][k] = g.getAccuracy();

			}
			for(int k=0; k<=i; k++){
				for(int j=k+1; j<=i; j++){
					Q[k][j] = Q[j][k];
				}
			}
			
			Q = MatrixCalculus.invert(Q);
			Q = MatrixCalculus.multiply(P, Q);
			
			
			prob0 = (float) MatrixCalculus.multiply(Q, MatrixCalculus.transpose(P))[0][0];
			if(Float.isNaN(prob0)){
				prob0 = 0;
			}
		}
		//System.out.println(prob0);
		//de Caen
		float prob1 = 0;
		/*
		for(int k=0;k<=i;k++){
			NewGappedPeptideForMSGappedNovo g2 = gappedPeptides.get(k);
			
			float nominator = g2.getAccuracy() * g2.getAccuracy();
			float denominator = g2.getAccuracy();
			float beta = 0;
			
			for(int j=0; j<=i; j++){
				if(k==j) continue;
				
				NewGappedPeptideForMSGappedNovo w = gappedPeptides.get(j);
				
				ArrayList<Integer> key = new ArrayList<Integer>();
				if(k<j){
					key.add(k); key.add(j);
				}else{
					key.add(j); key.add(k);
				}
				Float acc = unionGPAccuracy.get(key);
				
				if(acc == null){
					NewGappedPeptideForMSGappedNovo u = g2.getUnionGappedPeptide(gappedPeptides.get(j));
					acc =u.getAccuracy();
					unionGPAccuracy.put(key, acc);
				}
				
				float plus = Math.min(g2.getAccuracy(), w.getAccuracy());
				plus = Math.min(acc, plus);
				beta += plus;
				denominator += plus;
			}
			
			float theta = beta/g2.getAccuracy() - (int) (beta/g2.getAccuracy());

			prob1 += theta * nominator / (denominator + (1-theta) * g2.getAccuracy())
				+ (1-theta)*nominator / (denominator - theta * g2.getAccuracy());
		}
		
		if(Float.isNaN(prob1)){
			prob0 = 1;
		}
		
		*/
		float prob2 = 0;
		/*
		 if(i>0){ // Kounias
			NewGappedPeptideForMSGappedNovo g = gappedPeptides.get(0);
			for(int j=0;j<i;j++){
				NewGappedPeptideForMSGappedNovo gb = gappedPeptides.get(j);
				
				ArrayList<Integer> key = new ArrayList<Integer>();
				if(i<j){
					key.add(i); key.add(j);
				}else{
					key.add(j); key.add(i);
				}
				Float acc = unionGPAccuracy.get(key);
				
				if(acc == null){
					NewGappedPeptideForMSGappedNovo u = g.getUnionGappedPeptide(gappedPeptides.get(j));
					acc =u.getAccuracy();
					unionGPAccuracy.put(key, acc);
				}
				
				float w = Math.min(acc, g.getAccuracy());
				w = Math.min(w, gb.getAccuracy());

				prob2 = Math.max(prob2, g.getAccuracy() + gb.getAccuracy() - w);
			}
		}
		
		if(Float.isNaN(prob2)){
			prob2 = 0;
		}
		*/
		accuracy = Math.max(prob0, prob1);
		accuracy = Math.max(accuracy, prob2);
		//System.out.println(prob0 >= prob1);
		return accuracy;
	}
	
	private void backtrack(ArrayList<GappedPeptideForMSGappedNovo> rec, ArrayList<Node> suffix, HashMap<Node, HashMap<Integer, HashMap<Integer, HashSet<Edge>>>> es, Node v, int l, int s, int maxRecNum){
		if(rec.size() > maxRecNum) return;
		HashMap<Integer, HashMap<Integer, HashSet<Edge>>> subes = es.get(v);
		if(subes == null) return;
		HashMap<Integer, HashSet<Edge>> ssubes = subes.get(l);
		if(ssubes == null) return;
		HashSet<Edge> edges = ssubes.get(s);
		if(edges == null) return;
		
		if(v.equals(Node.getNode(0f))){// source
			suffix.add(v);
			Collections.sort(suffix);
			rec.add(new GappedPeptideForMSGappedNovo(suffix, graph.getNodes()));
			return;
		}
		
		for(Edge edge : edges){
			ArrayList<Node> nextSuffix = new ArrayList<Node>();
			nextSuffix.addAll(suffix);
			nextSuffix.add(v);
			backtrack(rec, nextSuffix, es, edge.getLeftNode(), l-1, s-edge.getRoundedLogAccuracy(), maxRecNum);
		}
	}
	
	public static void main(String[] args){
		String specfilename = "/home/kwj/workspace/inputs/Training/CID_Tryp_Confident.mgf";
		String trainfilename ="/home/kwj/workspace/inputs/Training/shewLengthAll.mgf";
		String para = trainfilename.substring(0, trainfilename.lastIndexOf(".")) + ".par";
		
		int charge  = 2;
		int numPerLength = 10;
		int minGappedPeptideLength = 5;
		float minGappedPeptideAccuracy = 0;
		
		Iterator<Spectrum> iterator;
		Tolerance tol = new Tolerance(0.5f, false);
		Tolerance pmtol = new Tolerance(0.5f, false);
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			IPE ipe = new IPE(para, tol, pmtol, aaSet).filter(new WindowFilter(8, 50));
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				if(spec.getCharge() != charge) continue;
				
				AugmentedSpectrumGraph graph = new AugmentedSpectrumGraph(spec, ipe, para, 100);
				/*
				for(NewNode node : graph.getNodes()){
					System.out.println("Node " + node.getMass() + " acc " + node.getAccuracy());
					for(Edge e : graph.getEdges(node)){
						System.out.println("\t Edge " + e.getLeftNode().getMass() + " : " + e.getRightNode().getMass() + " acc " + AugmentedSpectrumGraph.getAccuracyFromRoundedLogAccuracy(e.getAccuracy()));
					}
				}*/
				MSGappedNovo mg = new MSGappedNovo(graph);
				HashMap<Integer, ArrayList<GappedPeptideForMSGappedNovo>> can =  mg.getCandidateGappedPeptides(numPerLength, minGappedPeptideLength, minGappedPeptideAccuracy);
				
				/*for(int l : can.keySet()){
					for(NewGappedPeptideForMSGappedNovo gp : can.get(l)){
						System.out.println(spec.getAnnotationStr() + "\t" + l + "\t" + gp + "\t" + gp.getAccuracy());
					}
				}*/
				ArrayList<GappedPeptideForMSGappedNovo> out = new ArrayList<GappedPeptideForMSGappedNovo>();
				System.out.println(spec);
				mg.selectOutputFromCandidates(out, can, numPerLength, minGappedPeptideLength, 0.9f);
				for(GappedPeptideForMSGappedNovo gp : out){
					System.out.println(spec.getAnnotationStr() + "\t" + gp + "\t" + gp.length() + "\t" + gp.getAccuracy());
					float a = 1;
					for(Node n : gp.getNodes()){
						a *= n.getAccuracy();
					}
				//	System.out.println("\t" + a);
				}
				
				//break;
			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
	}
}
