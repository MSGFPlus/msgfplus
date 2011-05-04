package msgappednovo;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.PriorityQueue;

import msgappednovo.AugmentedSpectrumGraph.Edge;
import msgappednovo.AugmentedSpectrumGraph.Node;
import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.Composition;
import msutil.Peptide;

public class GappedPeptideForMSGappedNovo extends BitSet implements Comparable<GappedPeptideForMSGappedNovo>{ // excludes the source node

	private static final long serialVersionUID = 1L;
	private ArrayList<Node> allNodes;
	private float accuracy = 0;
	private String s;
	private static int alpha = 100;
	
	static public void setAlpha(int alpha){
		GappedPeptideForMSGappedNovo.alpha = alpha;
	}
	
	private void setHiddenNodes(){
		/*HashMap<NewNode, ArrayList<Integer>> hiddenNodeTable = AugmentedSpectrumGraph.getHiddenNodeTable();
		
		NewNode prevNode = NewNode.getNewNode(0f);
		
		for(NewNode node : getNodes()){
			NewNode gapNode = NewNode.getNewNode(node.getBinNumber()-prevNode.getBinNumber());
			if(hiddenNodeTable.containsKey(gapNode)){
				ArrayList<Integer> indices = hiddenNodeTable.get(gapNode);
				for(int index : indices){
					int i = Collections.binarySearch(allNodes, NewNode.getNewNode(prevNode.getBinNumber() + index));
					if(i>0) this.set(i);
				}
			}
			prevNode = node;
		}*/
	}
	
	public GappedPeptideForMSGappedNovo(ArrayList<Node> nodes, ArrayList<Node> allNodes){
		this.allNodes = allNodes;
	//	NewNode prev = NewNode.getNewNode(0f);
	//	int l = 1;
		//float acc = 1;
		//float minNodeAcc = 1;
		for(Node node : nodes){
			int i = Collections.binarySearch(allNodes, node);
			if(node.getBinNumber() > 0 && i>=0){
				this.set(i);
			//	acc *= Math.min(1, Edge.getEdge(prev, node).getAccuracy());// + prev.getAccuracy()/acc * 0.15 * l++);
			}
	
		//	prev = node;
		}
		//acc = (float) (1- Math.pow((1-acc), 0.8f*(1 + (this.cardinality()-4))));
		
		float acc = 0;
		
		for(int i=0;i<getNodes().size()-1;i++){
			Node node = getNodes().get(i);
			acc += node.getAccuracy();
		//	minNodeAcc = Math.min(minNodeAcc, node.getAccuracy());
		}
			
		for(Edge e : getMST2()){
			acc -= e.weightForMST();	
		}
		
		acc = Math.max(0, acc);
		acc = Math.min(1, acc);
		
		this.accuracy = acc;
		setHiddenNodes();
		
		/*
		if(this.length() == 6){
			float x= 0 ;
			for(Edge e : getMST2()){
				System.out.println("*" + e.getLeftNode().getBinNumber() + " " + e.getRightNode().getBinNumber() + " " + e.weightForMST());
				x +=   e.weightForMST();
			}
			System.out.println("********* " + x);
			x=0;
			for(Edge e : getEdges()){
				if(e.getLeftNode().equals(NewNode.getNewNode(0f))) continue;
				if(e.getRightNode().equals(getNodes().get(getNodes().size()-1))) continue;
				System.out.println(e.getLeftNode().getBinNumber() + " " + e.getRightNode().getBinNumber() + " " +  e.weightForMST());
				x += e.weightForMST();
			}
			System.out.println("*** " + x);
			System.exit(0);
			//System.out.println();
		}*/
	}
	
	private GappedPeptideForMSGappedNovo(GappedPeptideForMSGappedNovo one, GappedPeptideForMSGappedNovo two){
		super();
		super.or(one); super.or(two);
		
		this.allNodes = one.allNodes;
	
		float acc = 0.0f; 
		//float minNodeAcc = 1;
		Node prev = Node.getNode(0);
		
		for(int i=0;i<getNodes().size()-1;i++){
			Node node = getNodes().get(i);
			acc += node.getAccuracy(); 
			Edge e = Edge.getEdge(prev, node);
			if(e.getAccuracy() == 0){
				this.accuracy = 0.05f; // to be conservative
				return;
			}
			prev = node;
		//	minNodeAcc = Math.min(minNodeAcc, node.getAccuracy());
		}
		
		for(Edge e : getMST2()){
			acc -= e.weightForMST();	
		}

		acc = Math.max(0, acc);
		acc = Math.min(1, acc);

		this.accuracy = acc;
	}
	
	
	public GappedPeptideForMSGappedNovo getUnionGappedPeptide(GappedPeptideForMSGappedNovo other){
		return new GappedPeptideForMSGappedNovo(this, other);
	}
	
	public float getParentMass() { return (float) (this.allNodes.get(allNodes.size()-1).getMass() + Composition.H2O);}
	
	public ArrayList<Node> getNodes(){//excludes source
		ArrayList<Node> nodes = new ArrayList<Node>();
		for (int b = this.nextSetBit(0); b >= 0; b = this.nextSetBit(b+1)){
			nodes.add(allNodes.get(b));
		}
		return nodes;
	}
	
	public ArrayList<Edge> getEdges(){
		ArrayList<Edge> ret = new ArrayList<Edge>();
		Node prev = Node.getNode(0f);
		
		for(Node node : getNodes()){
			Edge e = Edge.getEdge(prev, node);
			ret.add(e);
			prev = node;
		}
		return ret;
	}
	
	public float getAccuracy() {return accuracy;}
	
	public boolean isCorrect(Peptide pep, Tolerance tol, Tolerance pmtol){
		if(Math.abs(this.getNodes().get(this.getNodes().size()-1).getMass() - pep.getMass()) > pmtol.getToleranceAsDa(pep.getMass()))
				return false;
		ArrayList<Node> nodes = this.getNodes();
		ArrayList<Float> prms = new ArrayList<Float>();
		
		int prmIndex = 0;
		boolean isCorrect = true;
		
		for(float prm : pep.getPRMMasses(true, 0)){
			prms.add(prm);
		}
		
		for(int i=0;i<nodes.size()-1;i++){
			float m = nodes.get(i).getMass();
			float delta = tol.getToleranceAsDa(Math.max(m, pep.getParentMass()-m)) * 2; 
			
			while(prmIndex < prms.size() && prms.get(prmIndex) < m - delta){
				prmIndex++;
			}
			
			if(prmIndex >= prms.size() || Math.abs(m - prms.get(prmIndex)) > delta){
				isCorrect = false; break;
			}
			
		}
		
		return isCorrect;
	}
	
	public String toString(){
		if(this.s!=null) return this.s;
		String s = "";
		Node prev = Node.getNode(0f);

		for (Node m : this.getNodes()) {
			Node gapNode = Node.getNode(m.getMass() - prev.getMass());
			AminoAcid aa = AugmentedSpectrumGraph.getAABinDiffs().get(gapNode.getBinNumber());
			 if(aa!= null){
				s += aa.getResidueStr();
			 }else
				s += "["+gapNode.getMass()+"]";
			 prev = m;
		 }
		 
		this.s =  s;
		return this.s;
	}
	
	public int length() {
		return this.cardinality();
	}
	
	private ArrayList<Edge> getMST(){//TODO optimize!!
	
		ArrayList<Edge> MST = new ArrayList<Edge>();
		HashSet<Node> treeNodes = new HashSet<Node>();
		ArrayList<Node> allNodes = getNodes();
		allNodes.remove(allNodes.size()-1); // exclude sink
	//	System.out.println(allNodes);
		PriorityQueue<Edge> queue = new PriorityQueue<Edge>();
		treeNodes.add(getNodes().get(0)); //beginning
		Node lastAdded = getNodes().get(0);
		
		while(MST.size() < length()-2){
			allNodes.remove(lastAdded);
			//System.out.println(" " + allNodes);
			for(Node node : allNodes){
				Edge e;
				if(node.getBinNumber() < lastAdded.getBinNumber()){
					e = Edge.getEdge(node, lastAdded);
				}else{
					e = Edge.getEdge(lastAdded, node);
				}
				queue.add(e);
			}
			
			Edge minEdge = null;
			while(true){
			minEdge = queue.poll();
				if(!treeNodes.contains(minEdge.getRightNode())){
					lastAdded = minEdge.getRightNode();
					break;
				}else if(!treeNodes.contains(minEdge.getLeftNode())){
					lastAdded = minEdge.getLeftNode();
					break;
				}else{
				//	System.out.println("NO" + lastAdded.getBinNumber());
				}
			}
			MST.add(minEdge);
			treeNodes.add(lastAdded);
		}
		
		
		return MST;
	}
	
	
	private class AttachmentComparator implements Comparator<Node>{
		@Override
		public int compare(Node o1, Node o2) {
			return new Float(o1.getAttachmentCost()).compareTo(new Float(o2.getAttachmentCost()));
		}
	}

	private ArrayList<Edge> getMST2(){//TODO optimize!!
		
		ArrayList<Node> nodes = getNodes();
		PriorityQueue<Node> queue = new PriorityQueue<Node>(nodes.size(), new AttachmentComparator());
		ArrayList<Edge> MST = new ArrayList<Edge>();
		ArrayList<Node> nodesNotInTree = new ArrayList<Node>();
		int[] nodeIndices = new int[nodes.get(nodes.size()-1).getBinNumber()+1];
		
		for(int i=0; i<nodes.size()-1;i++){// init
			Node node = nodes.get(i);
			if(i==0){
				node.setAttachmentCost(0f);
			}else{
				node.setAttachmentCost(Float.MAX_VALUE);
			}
			nodesNotInTree.add(node);
			queue.add(node);
		}
		
		while(!nodesNotInTree.isEmpty()){
			
			Node node = queue.poll();
			int connectingNodeIndex = nodeIndices[node.getBinNumber()];
			
			if(connectingNodeIndex != 0){
				Edge e;
				if(node.getBinNumber() < connectingNodeIndex){
					e = Edge.getEdge(node, Node.getNode(connectingNodeIndex));
				}else{
					e = Edge.getEdge(Node.getNode(connectingNodeIndex), node);
				}
				MST.add(e);
			}
			nodesNotInTree.remove(node);
			
			for(Node n : nodesNotInTree){
				Edge e;
				if(n.getBinNumber() < node.getBinNumber()){
					e = Edge.getEdge(n, node);
				}else{
					e = Edge.getEdge(node, n);
				}
				if(e.weightForMST() < n.getAttachmentCost()){// change key
					queue.remove(n);
					n.setAttachmentCost(e.weightForMST());
					queue.add(n);
					nodeIndices[n.getBinNumber()] = node.getBinNumber();
				}
			}
			
			
		}
		return MST;
	}
	
	@Override
	public int compareTo(GappedPeptideForMSGappedNovo arg0) {
	//	int c = new Float(this.length()).compareTo(new Float(arg0.length()));
		//if(c == 0)	
			return new Float(this.accuracy + alpha * this.length()).compareTo(new Float(arg0.accuracy + alpha * arg0.length()));
		//else return c;
		//	return new Float(this.accuracy).compareTo(new Float(arg0.accuracy));
			
	}
	
	
	
}
