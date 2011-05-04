package msgappednovo;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import parser.BufferedLineReader;
import parser.MgfSpectrumParser;

import msgappednovo.train.PeakGenerator;
import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Constants;
import msutil.IonType;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;

public class AugmentedSpectrumGraph {
		
	static private boolean useCovariance = false;
	
	final static private float minNodeAccuracy = 0.05f;
	final static private float maxEdgeMass = 3000f;
	
	static private boolean highPrecision = false;
	static private Tolerance tol, pmtol;// only used when highPrecision is on
	static private AminoAcidSet aaSet;
	static private HashMap<Node, Integer> minAANumTable;
	static private HashMap<Node, ArrayList<Integer>>hiddenNodeTable;
	static private HashMap<Integer, AminoAcid> aaBinDiffs;
	
	static public void turnHighPrecisionModeOn(Tolerance tol, Tolerance pmtol) {
		highPrecision = true;
		AugmentedSpectrumGraph.tol = tol;
		AugmentedSpectrumGraph.pmtol = pmtol;
	}
	static private float[][][] weights;
	static private float[][] offset;
	static private float[][][][] covariances;
	static private float [][][][] edgeAccuracies;
	
	static public int getRoundedLogAccuracyFromAccuracy(float accuracy){
		return (int) Math.round(Math.log10(accuracy) * 100);
	}
	
	static public HashMap<Node, float[]> getNodeIonProbs(Spectrum spec, IPE ipe, int maxRank){
		HashMap<Peak, HashMap<IonType, Float>> profile = ipe.getProfile(spec, maxRank, Integer.MAX_VALUE);
		ArrayList<IonType> sigIons = ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(spec.getCharge());
		HashMap<Node, float[]> nodeIonProbs = new HashMap<Node, float[]>();
		//new Node(spec.getParentMass() - (float)Composition.H2O).getBinNumber() + 1
		float sinkMass = spec.getParentMass() - (float)Composition.H2O;
		
		for(Peak peak : profile.keySet()){
			HashMap<IonType, Float> prof = profile.get(peak);
			for(int i=0; i<sigIons.size();i++){
				IonType ion = sigIons.get(i);
				float mass = PeakGenerator.getPrefixMass(peak, ion, spec);
				if(mass >= sinkMass || mass <= 0) continue;
				Node node = Node.getNode(mass);//TODO high precision mode?
				
				if(!nodeIonProbs.containsKey(node)) nodeIonProbs.put(node, new float[sigIons.size()]);
				float[] probs = nodeIonProbs.get(node);
				probs[i] += prof.get(ion);
			}
		}
		
		float[] sourceSinkProbs = new float[sigIons.size()];
		for(int i=0; i<sigIons.size();i++){sourceSinkProbs[i] = 1;}
		nodeIonProbs.put(Node.getNode(sinkMass), sourceSinkProbs);
		nodeIonProbs.put(Node.getNode(0), sourceSinkProbs);
		
		return nodeIonProbs;
	}
	static public ArrayList<Node> getNodes(HashMap<Node, float[]> nodeIonProbs, float[][][] weights, float[][] offset, int charge){
		
		ArrayList<Node> nodes = new ArrayList<Node>();
		
		for(Node node : nodeIonProbs.keySet()){
			float[] ionProbs = nodeIonProbs.get(node);
			int caseIndex=0;
			int cardinality = 0;
			int sigIonSize = ionProbs.length;
			
			while(sigIonSize > 0){
				caseIndex=0;
				cardinality = 0;
				for(int j=0; j<sigIonSize; j++){ // ion
					if(ionProbs[j] > 0){
						caseIndex += 1<<j;
						cardinality++;
					}
				}
				if(weights[charge][caseIndex] != null) break;
				sigIonSize--;
			}

			if(cardinality==0) continue;
			int[] indices = new int[cardinality];
			int t = 0;
			
			for(int j=0; j<sigIonSize; j++){ // ion
				if(ionProbs[j] > 0){
					indices[t++] = j;
				}
			}

			float weightedSum = offset[charge][caseIndex];
			for(int j=0; j<cardinality; j++){
				weightedSum+= ionProbs[indices[j]] * weights[charge][caseIndex][j];
			}
			
			weightedSum = Math.min(1, weightedSum);
			weightedSum = Math.max(0, weightedSum);
			
			if(weightedSum >= minNodeAccuracy){
				node.setAccuracy(weightedSum);
				nodes.add(node);
			}
		}
		Collections.sort(nodes);
		
		return nodes;
	}
	
	static private void read(String para){
		if(weights!=null) return;
		
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(para);
			String s;
			int mode = -100;
			int index = 0;
			int index2 = 0;
			int charge = 0;
			
			while((s=in.readLine())!=null){
				if(s.startsWith("#IONWEIGHT")){
					charge = Integer.parseInt(s.split("\t")[1]);
					continue;
				}else if(s.startsWith("#ENDIONWEIGHT")){
					mode = -100;
					continue;					
				}
				else if(s.startsWith("#OFF")){
					mode = 0;
					if(offset == null){
						offset = new float[100][];
					}
					if(offset[charge] == null)
						offset[charge] = new float[Integer.parseInt(s.split("\t")[2])];
					
					index = Integer.parseInt(s.split("\t")[1]);
					continue;
				}else if(s.startsWith("#WEIGHTS")){
					mode = 1;
					//	weightsForEachIon = new float[Integer.parseInt(s.split("\t")[1])][2];index = 0;
					if(weights == null){
						weights = (new float[100][][]);
					}
					if(weights[charge] == null)
						weights[charge] = new float[Integer.parseInt(s.split("\t")[3])][];
	
					index = 0;
					index2 = Integer.parseInt(s.split("\t")[1]);
					weights[charge][index2] = new float[Integer.parseInt(s.split("\t")[2])];
					continue;
				}else if(s.startsWith("#COVARIANCE")){
					mode = 2;
					index = 0; index2 = 0;
					charge = Integer.parseInt(s.split("\t")[4]);
					if(covariances == null)
						covariances = new float[100][][][];
					
					if(covariances[charge] == null)
						covariances[charge] = new float[Integer.parseInt(s.split("\t")[1])][Integer.parseInt(s.split("\t")[2])][Integer.parseInt(s.split("\t")[3])];
					continue;
				}else if(s.startsWith("#ENDCOVARIANCE")){
					mode = -100;
					continue;
				}else if(s.startsWith("#EDGEACCURACY")){
					mode = 3;
					index = 0; index2 = 0;
					charge = Integer.parseInt(s.split("\t")[4]);
					if(edgeAccuracies == null)
						edgeAccuracies = new float[100][][][];
					
					if(edgeAccuracies[charge] == null)
						edgeAccuracies[charge] = new float[Integer.parseInt(s.split("\t")[1])][Integer.parseInt(s.split("\t")[2])][Integer.parseInt(s.split("\t")[3])];
					continue;
				}else if(s.startsWith("#ENDEDGEACCURACY")){
					mode = -100;
					continue;
				}
				
				
	
				if(mode == -1){
					continue;
				}else if(mode == 0){
					offset[charge][index] = Float.parseFloat(s);
					continue;
				}else if(mode == 1){
					weights[charge][index2][index] = Float.parseFloat(s);
					index++;
				}else if(mode == 2){
					if(s.startsWith("##NUM")){
						index = Integer.parseInt(s.split("\t")[1]);
						index2 = 0;
						continue;
					}else{
						int index3=0;
					
						for(String t:s.split(" ")){
							if(!t.isEmpty()){
								covariances[charge][index][index2][index3++] = Float.parseFloat(t);
							}
						}
						index2++;
					}
				}else if(mode == 3){
					if(s.startsWith("##NUM")){
						index = Integer.parseInt(s.split("\t")[1]);
						index2 = 0;
						continue;
					}else{
						int index3=0;
					
						for(String t:s.split(" ")){
							if(!t.isEmpty()){
								edgeAccuracies[charge][index][index2][index3++] = Float.parseFloat(t);
							}
						}
						index2++;
					}
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	
	static private void updateTables(AminoAcidSet aaSet){
		if(AugmentedSpectrumGraph.aaSet == aaSet) return;
		
		minAANumTable = new HashMap<Node, Integer>();
		minAANumTable.put(Node.getNode(0f), 0);
		setHiddenNodeTable(new  HashMap<Node, ArrayList<Integer>>());
		getHiddenNodeTable().put(Node.getNode(0f), new ArrayList<Integer>());
		
		AugmentedSpectrumGraph.aaBinDiffs = new HashMap<Integer, AminoAcid>();
		
		for(AminoAcid aa : aaSet){
			getAABinDiffs().put(Node.getNode(aa.getMass()).getBinNumber(), aa);
		}
		
		for(int i=0; i<=Node.getNode(maxEdgeMass).getBinNumber();i++){
			int m = Integer.MAX_VALUE;
			int numEdge = 0;
			int prev = 0, cprev = 0;
			
			for(int j : getAABinDiffs().keySet()){
				prev = i-j;
				
				if(minAANumTable.containsKey(Node.getNode(prev))){
					m = Math.min(m, minAANumTable.get(Node.getNode(prev)));
					numEdge++;
					cprev = prev;
				}
			}
			
			if(numEdge == 1) {
				ArrayList<Integer> indices = new ArrayList<Integer>();
				indices.addAll(getHiddenNodeTable().get(Node.getNode(cprev)));
				indices.add(i);
				getHiddenNodeTable().put(Node.getNode(i), indices);
			}
			
			if(m<Integer.MAX_VALUE){
				minAANumTable.put(Node.getNode(i), m+1);
			}
		}	
	}
	
	private ArrayList<Node> nodes;
	private Node sinkNode;
	private HashMap<Node, ArrayList<Edge>> adjList;
	
	
	private void setEdgeAccuracy(Edge e, int charge){
		float acc;
		
		Node gapNode = Node.getGapNode(e);
		Integer min = minAANumTable.get(gapNode);
		float accR = e.r.getAccuracy();
		float accL = e.l.getAccuracy();
		
		if(min == null || min == 0){
			acc = 0;
		}else if(accR >= 1 || e.getRightNode().equals(getSinkNode())){
			acc = 1;
		}else if(edgeAccuracies != null && !useCovariance){
			int i1 = Math.min(edgeAccuracies[charge].length-1, min);
			int i2_1 = (int) ((edgeAccuracies[charge][0].length-1)*accR);
			int i2_2 = i2_1 + 1;
			int i3 = Math.round((edgeAccuracies[charge][0].length-1)*accL);
			float acc1 = edgeAccuracies[charge][i1][i2_1][i3];
			float acc2 = i2_2 < edgeAccuracies[charge][i1].length ? edgeAccuracies[charge][i1][i2_2][i3] : 1;
			
			acc = acc1 + (acc2-acc1)*((edgeAccuracies[charge][0].length-1)*accR - i2_1);
		
		}else if(covariances != null && useCovariance){
			int i1 = Math.min(covariances[charge].length-1, min);
			int i2 = Math.round((covariances[charge][0].length-1)*accL);
			int i3 = Math.round((covariances[charge][0].length-1)*accR);
			float cov = covariances[charge][i1][i2][i3];
			
			acc = accR + cov/accL;
		}else{
			acc = e.getRightNode().getAccuracy();
		}
		
		acc = Math.max(acc, 0);
		acc = Math.min(acc,1);
		
		e.setAccuracy(acc);
	}
	
	private void updateAdjList(int charge){
		adjList = new HashMap<Node, ArrayList<Edge>>();
	
		for(Node node : nodes){
			ArrayList<Edge> connectingEdges = new ArrayList<Edge>();
			for(Node cnode : nodes){
				if(node.getMass() <= cnode.getMass()) continue;
				
				Edge connectingEdge = Edge.getEdge(cnode, node);
				setEdgeAccuracy(connectingEdge, charge);
				if(connectingEdge.getAccuracy() > 0 ){
					connectingEdges.add(connectingEdge);
				}
			}

			if(!connectingEdges.isEmpty()) adjList.put(node, connectingEdges);
	
		}
		
	
		if(highPrecision){//TODO sink edge update
			
		}
	}
	
	
	public ArrayList<Node> getNodes(){return nodes;}
	public ArrayList<Edge> getEdges(Node node){
		ArrayList<Edge> edges;
		edges = adjList.get(node);
		
		if(edges == null) return new ArrayList<Edge>();
		else	return edges;
	}
	
	public AugmentedSpectrumGraph(Spectrum spec, IPE ipe, String para, int maxRank){
		read(para);
		updateTables(ipe.getAASet());
		Node.clearAccuracy();
		Edge.clearAccuracy();
		HashMap<Node, float[]> nodeIonProbs = getNodeIonProbs(spec, ipe, maxRank);
		nodes = getNodes(nodeIonProbs, weights, offset, spec.getCharge());
		sinkNode = nodes.get(nodes.size()-1);
		updateAdjList(spec.getCharge());
	}

	public void setSinkNode(Node sinkNode) {
		this.sinkNode = sinkNode;
	}

	public Node getSinkNode() {
		return sinkNode;
	}
	public static void setHiddenNodeTable(HashMap<Node, ArrayList<Integer>> hiddenNodeTable) {
		AugmentedSpectrumGraph.hiddenNodeTable = hiddenNodeTable;
	}
	public static HashMap<Node, ArrayList<Integer>> getHiddenNodeTable() {
		return hiddenNodeTable;
	}
	public static HashMap<Integer, AminoAcid> getAABinDiffs() {
		return aaBinDiffs;
	}
	public static HashMap<Node, Integer> getMinAANumTable() {
		return minAANumTable;
	}
	
	public static void main(String[] args){
		String specfilename ="/home/kwj/workspace/inputs/Training/shewLengthAll.mgf";
		String para = specfilename.substring(0, specfilename.lastIndexOf(".")) + ".par";
		
		Iterator<Spectrum> iterator;
		Tolerance tol = new Tolerance(0.5f, false);
		Tolerance pmtol = new Tolerance(0.5f, false);
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			IPE ipe = new IPE(para, tol, pmtol, aaSet).filter(new WindowFilter(8, 50));
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				AugmentedSpectrumGraph mag = new AugmentedSpectrumGraph(spec, ipe, para, 100);
				for(Node node : mag.getNodes()){
					System.out.println("Node " + node.getMass() + " acc " + node.getAccuracy());
					for(Edge e : mag.getEdges(node)){
						System.out.println("\t Edge " + e.getLeftNode().getMass() + " : " + e.getRightNode().getMass() + " acc " + e.getAccuracy());
					}
				}
				
			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
		
		
		//NewMassAccuracyGraph(Spectrum spec, IPE ipe, String para, int maxRank)
	}
	
	static public class Node implements Comparable<Node>{
		private float mass;
		private int binNumber;
		private float attachmentCost = Float.MAX_VALUE;
		
		static private HashMap<Node, Float >accuracy = new HashMap<Node, Float >();
	
		private Node(float mass){
			
			if(highPrecision){//TODO
				
			}else{
				binNumber =  Math.round(mass * Constants.INTEGER_MASS_SCALER);
				this.mass = binNumber/Constants.INTEGER_MASS_SCALER;
			}
		}
		
		private Node(int binNumber){
			this.binNumber = binNumber;
			if(highPrecision){//TODO
				
			}else{
				this.mass = binNumber/Constants.INTEGER_MASS_SCALER;
			}
		}
		
		static public Node getNode(float mass){ return new Node(mass); }
		static public Node getNode(int binNumber) { return new Node(binNumber); }
		static public Node getGapNode(Edge e){
			return Node.getNode(e.r.binNumber -e.l.binNumber);
		}
		static public void clearAccuracy() {accuracy.clear();}
		
		public float getMass(){
			return mass;
		}
		
		public int getBinNumber(){
			return binNumber;
		}
		
		public float getAccuracy(){
			return accuracy.get(this);
		}
		
		public void setAccuracy(float acc){
			accuracy.put(this, acc);
		}
		
		public int hashCode(){
			return getBinNumber();
		}
		
		@Override
		public boolean equals(Object obj) {
			if(obj == this) return true;
			if(obj instanceof Node){
				return ((Node)obj).getBinNumber() == this.getBinNumber();
			}else return false;
		}
		
		public String toString(){
			return new Float(this.getMass()).toString();
		}
		
		public int compareTo(Node o) {
			return ((Integer)this.getBinNumber()).compareTo(o.getBinNumber());
		}
		
		public float getAttachmentCost(){ // for prim's algorithm
			return attachmentCost;
		}
		
		public void setAttachmentCost(float attachmentCost){
			this.attachmentCost = attachmentCost;
		}
		
	}
	
	static public class Edge implements Comparable<Edge>{
		static private HashMap<Edge, Integer> roundedAccuracy = new HashMap<Edge, Integer>();
		static private HashMap<Edge, Float> accuracy = new HashMap<Edge, Float>();
		
		private Node l=null, r=null;
		public Edge(Node l, Node r){
			this.l = l; this.r = r;
		}
		private Edge(){
			// null edge
		}
		
		public int getRoundedLogAccuracy(){
			Integer acc = roundedAccuracy.get(this);
			if(acc == null){
				if(this.getMass() > maxEdgeMass) return getRoundedLogAccuracyFromAccuracy(r.getAccuracy());
				else return Integer.MIN_VALUE/2;
			}
			return acc;
		}
		
		public float getAccuracy(){
			Float acc = accuracy.get(this);
			if(acc == null){
				if(this.getMass() > maxEdgeMass) return r.getAccuracy();
				else return 0;
			}
			return acc;
		}
		
		public float getMass() { return r.getMass() - l.getMass(); }
		
		static public Edge getNullEdge(){
			return new Edge();
		}
		
		static public Edge getEdge(Node l, Node r){
			return new Edge(l,r);
		}
		
		static public void clearAccuracy() {roundedAccuracy.clear(); accuracy.clear();}
		
		private void setAccuracy(float acc){
			accuracy.put(this, acc);
			roundedAccuracy.put(this, getRoundedLogAccuracyFromAccuracy(acc));
		}
		public Node getLeftNode() {return l;}
		public Node getRightNode() {return r;}
		
		public int hashCode(){
			if(l==null || r==null) return 0;
			return l.hashCode() * r.hashCode();
		}
		
		@Override
		public boolean equals(Object obj) {
			if(obj == this) return true;
			if(obj instanceof Edge){
				return ((Edge)obj).l.equals(this.l) && ((Edge)obj).r.equals(this.r);
			}else return false;
		}
		
		public float weightForMST(){
			return l.getAccuracy() + r.getAccuracy() - this.getAccuracy() * l.getAccuracy();
		}
		
		@Override
		public int compareTo(Edge arg0) {
			return new Float(this.weightForMST()).compareTo(new Float(arg0.weightForMST()));
		}
		
	}
	
}
