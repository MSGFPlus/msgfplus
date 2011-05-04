package suffixtree.trees;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.TreeMap;


import sequences.FastaSequence;
import suffixtree.Constants;
import suffixtree.edges.Edge;
import suffixtree.edges.ByteEdge;
import suffixtree.nodes.InternalNode;
import suffixtree.nodes.Node;


public class GappedSuffixTree extends SuffixTree {

  
  /**
   * This is the class to hold the information for merged items.
   * @author jung
   *
   */
  public static class MergeInfo {
    private CompositeEdge mergedEdge;
    private ArrayList<TandemEdge> members;
    
    public MergeInfo(CompositeEdge mergedEdge, ArrayList<TandemEdge> members) {
      this.mergedEdge = mergedEdge;
      this.members = members;
    }
    
    public boolean isContained(HashSet<TandemEdge> others) {
      for (TandemEdge n : members) if(!others.contains(n)) return false;
      return true;
    }
  }
  
  
  public GappedSuffixTree(FastaSequence sequence) {
    super(sequence);
    System.out.println("Loaded sequence into memory");
    makeGapLinks();
  }
  
  public GappedSuffixTree(FastaSequence sequence, boolean noInsert) {
    super(sequence, noInsert);
  }
  
  /**
   * An CompositeEdge represents an edged with sorted label sequences and 
   * represents all the permutations of this label.
   * @author jung
   *
   */
  public class CompositeEdge extends Edge {
    
    // this edge can encode a label of 4 bytes
    private int label;

    public CompositeEdge(int label, InternalNode sink) {
      this.label = label;
      setSink(sink);
    } 
    
    public CompositeEdge(byte[] label, InternalNode sink) {
      if (label.length>4) 
        System.err.println("Attempting to create a CompositeEdge with a label greater than 4 bytes... truncating.");
      this.label = 0;
      Arrays.sort(label);
      for (int i=0; i<label.length; i++) {
        this.label = this.label<<8;
        this.label |= label[i];
      }
      setSink(sink);
    }
    
    @Override
    public int getLabelAt(int offset) {
      if (offset > 0)
        System.err.println("Attempting to retrieve a label that does not exist... returning first label.");
      return label;
    }

    @Override
    public int size() {
      return 1;
    }
    
    /**
     * Get the actual length of the label.
     * @return the number of letters (bytes) that this label represents.
     */
    public int length() {
      int mask = 0xFF;
      for (int i=0; i<4; i++) {
        if ((label&mask)==0) return i;
        mask = mask<<8;
      }
      return 4;
    }

    @Override
    public Edge split(int offset) {
      System.err.println("Composite edge does not support splitting");
      return null;
    }

    @Override
    public int getEnd() {
      System.err.println("Unsupported method: get end position of CompositeEdge");
      System.exit(-1);
      return 0;
    }
    
    @Override
    public int getStart() {
      System.err.println("Unsupported method: get start position of CompositeEdge");
      System.exit(-1);
      return 0;
    }
    
    @Override
    public String toString() {
      FastaSequence sequence = getSequence();
      String result = "CompositeEdge: [";
      int mask = 0xFF000000;
      for (int i=0; i<32; i+=8) {
        byte item = (byte)((this.label&mask)>>>(24-i));
        if (item != 0)
          result += sequence.toChar(item);
        mask = mask>>>8;
      }
      return result + "] : " + label;
    }

    @Override
    public int mass() {
      System.err.println("GappedSuffixTree.CompositeEdge does not support mass()");
      return 0;
    }
  }
  
  
  
  /**
   * An tandem edge is copy edge of another compressed edge with an implicit node
   * splitting the given node.
   * @author jung
   *
   */
  public class TandemEdge extends Edge {

    private int label;       // the label of the first half
    private int end;         // the index of the end of the second edge
    private int start;       // the index of the start of the second edge
    
    public TandemEdge(Node sink, byte[] label, int start, int end) {
      Arrays.sort(label);
      if (label.length>4) 
        System.err.println("Attempting to create a composite edge with a label greater than 4 bytes... truncating.");
      this.label = 0;
      for (int i=0; i<label.length; i++) {
        this.label = this.label<<8;
        this.label |= label[i];
      }
      this.start = start;
      this.end = end;
      setSink(sink);
    }
    
    /**
     * Get the actual length of the label.
     * @return the number of letters (bytes) that this label represents.
     */
    public int length() {
      int mask = 0xFF;
      for (int i=0; i<4; i++) {
        if ((label&mask)==0) return i;
        mask = mask<<8;
      }
      return 4;
    }
    
    @Override
    public int getLabelAt(int offset) {
      if (offset > 0)
        System.err.println("Attempting to retrieve a label that does not exist... returning first label.");
      return label;
    }

    @Override
    public int size() {
      return 1;
    }

    @Override
    public Edge split(int offset) {
      System.err.println("Unsupported operation: split tandem edge " + toString() + " at " + offset);
      return null;
    }

    @Override
    public int getEnd() {
      System.err.println("Unsupported operation: Getting end position of tandem edge");
      System.exit(-1);
      return start;
    }
    
    @Override
    public int getStart() {
      System.err.println("Unsupported operation: Getting start position of tandem edge");
      throw new RuntimeException();
      //System.exit(-1);
      //return 0;
    }
    
    @Override
    public String toString() {
      FastaSequence sequence = getSequence();
      String result = "TandemEdge: [";
      int mask = 0xFF000000;
      for (int i=0; i<32; i+=8) {
        byte item = (byte)((this.label&mask)>>>(24-i));
        if (item != 0)
          result += sequence.toChar(item);
        mask = mask>>>8;
      }
      result += "] -o- ";
      result += sequence.getSubsequence(start, Math.min(start+6, end));
      result += "["+start+","+end+") #" + hashCode();
      //result += ". Parent " + this.parent;
      //result += ". Goes to sink " + getSink();
      //result += ". Original ";
      //for (byte b : original) result += b + String.format("[%s]", sequence.toString(b)) + ", ";
      return result;
    }
    
    @Override
    public int hashCode() {
      return getSink().hashCode();
    }
    
    @Override
    public boolean equals(Object other) {
      TandemEdge o = (TandemEdge)other;
      return this.getSink()==o.getSink() && this.start==o.start && this.end==o.end;
    }
    
    @Override
    public int mass() {
      System.err.println("GappedSuffixTree.TandemEdge does not support mass()");
      return 0;
    }
  }
  
    
  /**
   * Helper method and entry point to collect all descent of the given node.
   * Note that gap edges (tandem edges) have size of at least 2.
   * @param n the node to collect the children from. n must only have simple
   *          outgoing edges (single byte labels).
   * @param children the data structure to store the results.
   */
  public void getChildren(Node n, TreeMap<TandemEdge,ArrayList<TandemEdge>> children) {
    
    // collect all possible children that are within the gap cut off
    for (int i=0; i<n.getDegree(); i++) {
      Edge e = n.getEdgeAt(i);

      // clone the prefix so we don't mix labels of sibling edges
      ArrayList<Byte> currentPrefix = new ArrayList<Byte>();
      
      // Tandem edges need have the first edge of at least 2
      currentPrefix.add((byte)e.getLabelAt(0));
      
      // add all the tandem edges representing gaps to the children structure
      for (int j=1; j<e.size(); j++) {
        
        currentPrefix.add((byte)e.getLabelAt(j));
        if (currentPrefix.size()>Constants.MAX_GAP) break;
        
        int start = e.getStart() + currentPrefix.size();
        
        // copy the prefix into a new array to create the tandem edge
        byte[] prefixArray = new byte[currentPrefix.size()];
        for (int k=0; k<prefixArray.length; k++) prefixArray[k] = currentPrefix.get(k);
        TandemEdge te = new TandemEdge(e.getSink(), prefixArray, start, e.getEnd());
        
        if (children.containsKey(te)) {
          children.get(te).add(te);
        }
        else {
          ArrayList<TandemEdge> edgeList = new ArrayList<TandemEdge>();
          edgeList.add(te);
          children.put(te, edgeList);
        }
      }
      
      if (currentPrefix.size()<Constants.MAX_GAP) {
        // need to recurse and keep collecting children
        getChildren(e.getSink(), children, currentPrefix);
      }
    }
  }

  
  /**
   * Helper method to get all descendant of a node within the gap cutoff.
   * @param n the node to operate
   * @param children the data structure to store the results
   * @param prefix the bytes seen so far
   */
  private void getChildren(Node n, TreeMap<TandemEdge,ArrayList<TandemEdge>> children, ArrayList<Byte> prefix) {
    
    // collect all possible children that are within the gap cut off
    for (int i=0; i<n.getDegree(); i++) {
      Edge e = n.getEdgeAt(i);
      
      // clone the prefix so we don't mix labels of sibling edges
      ArrayList<Byte> currentPrefix = new ArrayList<Byte>(prefix);
      
      // add all the tandem edges representing gaps to the children structure
      for (int j=0; j<e.size(); j++) {
        
        currentPrefix.add((byte)e.getLabelAt(j));
        if (currentPrefix.size()>Constants.MAX_GAP) break;
        
        int start = e.getStart() + j + 1;
        
        byte[] prefixArray = new byte[currentPrefix.size()];
        for (int k=0; k<prefixArray.length; k++) prefixArray[k] = currentPrefix.get(k);
        
        TandemEdge te = new TandemEdge(e.getSink(), prefixArray, start, e.getEnd());
        
        if (children.containsKey(te)) {
          children.get(te).add(te);
        }
        else {
          ArrayList<TandemEdge> edgeList = new ArrayList<TandemEdge>();
          edgeList.add(te);
          children.put(te, edgeList);
        }
      }
      
      if (currentPrefix.size()<Constants.MAX_GAP) {
        // need to recurse and keep collecting children
        getChildren(e.getSink(), children, currentPrefix);
      }
    }
  }
  
  

  /**
   * Helper class to be able to wrap the Node into a priority queue sorted by 
   * their depth 
   */
  private static class NodeContainer implements Comparable<NodeContainer> {
    private int depth;
    private Node node;
    //private String path; // the path to this node container
    public NodeContainer(Node node, int depth) {
      this.node = node; this.depth = depth;
    }
    @Override
    public int compareTo(NodeContainer other) {
      return this.depth-other.depth;
    }
  }
  
  
  /**
   * Helper method to find out which is largest prefix length among the edgeList
   * @param commonEdges the list of edge objects. All edges should be compressed
   *        edges here.
   * @return
   */
  private static int getLongestPrefix(ArrayList<Edge> commonEdges) {
    int maxCommonIndex = 1;  // the first item is guaranteed to be the same
    while (true) {
      boolean completed = true;
      
      // we only have to check here because shortest edge must come first
      if (maxCommonIndex>=commonEdges.get(0).size()) {
        completed = false;
        break;
      }
      
      byte keyByte = (byte)commonEdges.get(0).getLabelAt(maxCommonIndex);
      for (int i=1; i<commonEdges.size(); i++) {
        if (keyByte!=(byte)commonEdges.get(i).getLabelAt(maxCommonIndex)) {
          completed = false;
          break;
        }
      }
      
      if (!completed) break;
      
      maxCommonIndex++;
    }
    return maxCommonIndex;
  }
  
  
  /**
   * Given a set of sorted edges, add them into a single node. If they are 
   * conflicting edges, merge them, creating internal nodes.
   * @param edgeList the list of edges to merge.
   * @return the node that has all the edgeList edges as outgoing edges 
   */
  private InternalNode transformToTree(ArrayList<Edge> edgeList) {
    
    Collections.sort(edgeList);
    
    /*
    System.err.println("Edge list as follows:");
    for (Edge item : edgeList) {
      System.err.println(item);
    }
    */
    
    // the new node created after merging the edgeList
    InternalNode n = new InternalNode();
    Edge previous = null;
    ArrayList<Edge> commonEdges = new ArrayList<Edge>();
    for (Edge e : edgeList) {
      
      if (!(e instanceof CompressedEdge)) {
        System.err.println("Found: " + e + " is not a Compressed edge in the tree transformation");
        System.exit(-1);
      }
      
      if (previous!=null && e.getLabel()!=previous.getLabel()) {
        // process the edges in the commonEdge list
        if (commonEdges.size()==1) {
          // easy case, there is a unique outgoing edge for this label
          n.insert(commonEdges.get(0));
        }
        else {
          // find out the largest common index for this set of common edges
          int maxCommonIndex = getLongestPrefix(commonEdges);
          
          //System.err.println("Size of conflict " + commonEdges.size());
          
          // create the new set of edges as a TreeSet
          ArrayList<Edge> outEdges = new ArrayList<Edge>();
          ArrayList<Integer> positions = new ArrayList<Integer>();
          for (Edge currentEdge : commonEdges) {
            //System.err.println("Processing: " + currentEdge);
            //System.err.println("Edge size: " + currentEdge.size() + " CommonIndex: " + maxCommonIndex);
            if (currentEdge.size()==maxCommonIndex) {
              for (int pos : currentEdge.getSink().getPositions()) positions.add(pos);
              // insert all outgoing edges of the given internal node.
              for (int i=0; i<currentEdge.getSink().getDegree(); i++) {
                outEdges.add(currentEdge.getSink().getEdgeAt(i));  
              }
            }
            else {
              outEdges.add(currentEdge.split(maxCommonIndex));
            }
          }
          
          // create a new edge for the collapse of this common edges
          int start = commonEdges.get(0).getStart();
          if (outEdges.size() > 0) {
            InternalNode compressedNode = transformToTree(outEdges);
            compressedNode.addPositions(positions);
            n.insert(createCompressedEdge(compressedNode, start, start+maxCommonIndex));
          }
          else {
            InternalNode compressedNode = new InternalNode(positions);
            n.insert(createCompressedEdge(compressedNode, start, start+maxCommonIndex));
          }
        }
        commonEdges.clear();
      }
      
      commonEdges.add(e);
      previous = e;
    }
    
    // process the last set of edges
    if (commonEdges.size()==1) {
      n.insert(commonEdges.get(0));
    }
    else {
      // find out the largest common index for this set of common edges
      int maxCommonIndex = getLongestPrefix(commonEdges);
      
      // create the new set of edges
      ArrayList<Edge> outEdges = new ArrayList<Edge>();
      ArrayList<Integer> positions = new ArrayList<Integer>();
      for (Edge currentEdge : commonEdges) {
        if (currentEdge.size()==maxCommonIndex) {
          for (int pos : currentEdge.getSink().getPositions()) positions.add(pos);
          // insert all outgoing edges of the given internal node.
          for (int i=0; i<currentEdge.getSink().getDegree(); i++) {
            outEdges.add(currentEdge.getSink().getEdgeAt(i));  
          }
        }
        else {
          outEdges.add(currentEdge.split(maxCommonIndex));
        }
      }
      
      // create a new edge for the collapse of this common edges
      int start = commonEdges.get(0).getStart();
      if (outEdges.size()>0) {
        InternalNode compressedNode = transformToTree(outEdges);
        compressedNode.addPositions(positions);
        n.insert(createCompressedEdge(compressedNode, start, start+maxCommonIndex));
      }
      else {
        InternalNode compressedNode = new InternalNode(positions);
        n.insert(createCompressedEdge(compressedNode, start, start+maxCommonIndex));
      }
    }
    commonEdges.clear();
  
    return n;
    
  }
  
  
  public void search(ArrayList<ByteEdge> query, HashSet<Integer> results) {
  
    Node currentNode = getRoot();
    if (currentNode==null) return;
    ByteEdge e = query.get(0);
    int i = 0;
 
    while (true) {
      
      //System.err.println("Matching: " + i + " - " + e.getLabel());
      
      int matchIndex = currentNode.search(e);
      if (matchIndex < 0) {
        // no match
        /*
        for (int j=0; j<currentNode.getDegree(); j++) {
          System.err.println(currentNode.getEdgeAt(j));
        }
        System.err.println("Rejected because no match for " + e.getLabel());
        */
        return;
      }
      
      Edge matchingEdge = currentNode.getEdgeAt(matchIndex);
      
      // we matched the current edge, try to advance
      i++;
      currentNode = matchingEdge.getSink();
      if (i>=query.size()) break; // found matching node
      e = query.get(i);           // update next query edge
      
      if (matchingEdge instanceof CompositeEdge) {
        // simply move on to the next internal node and continue matching
        continue;
      }
      
      // the matching edge is either a tandem edge or compressed edge
      int queryEdgeSize = e.length();
      int start = 0, end = 0;
      if (matchingEdge instanceof TandemEdge) {
        TandemEdge te = (TandemEdge)matchingEdge;
        start = te.start;
        end = te.end;
      }
      else {
        CompressedEdge ce = (CompressedEdge)matchingEdge;
        start = ce.getStart() + 1;     // the first one was already matched
        end = ce.getEnd();
        //System.err.println(start + ":" + end);
      }
      
      boolean done = false;
      while (end-start >= queryEdgeSize) {
        // match all you can starting at start
        CompositeEdge matchingSlice = new CompositeEdge(getSequence().getBytes(start, start+queryEdgeSize), null);
        if (matchingSlice.compareTo(e)!=0) {
          //System.err.println(matchingSlice.size() + " vs " + e.size());
          //System.err.println(matchingSlice + " - " + e.getLabel());
          //System.err.println("Rejected because no match");
          return;
        }
        
        // advance the query edge
        start += queryEdgeSize;
        i++;
        if (i>=query.size()) {
          done = true;
          break;
        }
        e = query.get(i);
        queryEdgeSize = e.length();
      }
      
      if (done) break;
      
      if (start==end) {
        // move on to the next internal node
        continue;
      }
      
      //System.err.println(e.getLabel());
      // match and remove the leftovers
      e = e.remove(getSequence().getBytes(start, end));
      if (e==null) {
        // no match
        //System.err.println("Removing " + getSequence().toString(start, end));
        //System.err.println("Rejected because no match");
        return;
      }
    }
    
    //System.err.println("Matching node is " + currentNode);
    // We found the matching node, return starting positions
    //System.out.println("Retrieving results on node with degree " + currentNode.getDegree());
    currentNode.getAllPositions(results);
    return;
    
  }
  
  
  public static ArrayList<CompositeEdge> reduceEdges(ArrayList<TandemEdge> edges, 
                                                     HashMap<TandemEdge,ArrayList<MergeInfo>> mergedNodes,
                                                     ArrayList<TandemEdge> remainingEdges) {
    HashSet<TandemEdge> targets = new HashSet<TandemEdge>();
    targets.addAll(edges);
    
    ArrayList<CompositeEdge> results = new ArrayList<CompositeEdge>();
    while(true) {
      MergeInfo chosenCandidate = null;
      for (TandemEdge mergeMember : edges) {
        if (mergedNodes.containsKey(mergeMember)) {
          for (MergeInfo candidate : mergedNodes.get(mergeMember)) {
            if (candidate.isContained(targets)) {
              if (chosenCandidate==null) {
                chosenCandidate = candidate;
              }
              else {
                if (chosenCandidate.members.size() > candidate.members.size()) {
                  chosenCandidate = candidate;
                }
              }
            }
          }
        }
      }
      // no more seen merges
      if (chosenCandidate==null) break;
      
      // remove the members of the chosen one from the targets
      targets.removeAll(chosenCandidate.members);
      results.add(chosenCandidate.mergedEdge);
      
      // we are done removing
      if (targets.size()==0) break;
    }
    
    // add the remaining edges
    for (TandemEdge te : edges) if (targets.contains(te)) remainingEdges.add(te);
    
    return results;
  }
  
  
  /**
   * Helper method to merge tandem edges with a common prefix. 
   * @param targetEdges the collection of edges that have a common first edge
   *              as they were hashed into the same bin by the hash function.
   * @return the resulting merged edge as a composite with the label explicitly 
   *         represented.
   */
  private CompositeEdge mergeEdges(ArrayList<TandemEdge> targetEdges, HashMap<TandemEdge,ArrayList<MergeInfo>> mergedNodes) {
    
    //mergedNodes.clear();
    
    
    ArrayList<TandemEdge> remainingEdges = new ArrayList<TandemEdge>();
    ArrayList<CompositeEdge> reusedEdges = reduceEdges(targetEdges, mergedNodes, remainingEdges);
    
    // collect all the outgoing edges from the split tandem edge array
    ArrayList<Edge> nextEdges = new ArrayList<Edge>();
    
    ArrayList<Integer> positions = new ArrayList<Integer>();
    for (TandemEdge e : remainingEdges) {
      TandemEdge te = (TandemEdge)e;
      Node nextNode = te.getSink();
      if (te.start==te.end) {
        // empty continuation edge will add all edges of the next sink node
        for (int i=0; i<nextNode.getDegree(); i++) {
          nextEdges.add(nextNode.getEdgeAt(i));
        }
        // collect the ending positions
        for (int pos : nextNode.getPositions()) {
          positions.add(pos);
        }
      }
      else {
        // otherwise just add the tandem edge as a new compressed edge
        nextEdges.add(new CompressedEdge(nextNode, te.start, te.end));
      }  
    }
    
    // add the edges of the observed merges
    for (CompositeEdge ce : reusedEdges) {
      Node sink = ce.getSink();
      for (int i=0; i<sink.getDegree(); i++) {
        nextEdges.add(sink.getEdgeAt(i));
      }
      for (int pos : sink.getPositions()) {
        positions.add(pos);
      }
    }
    
    // nextEdges are actually all compressed edges with a collided parent edge
    int label = targetEdges.get(0).getLabel();
    
    CompositeEdge result = null;
    if (nextEdges.size()>0) {
      InternalNode collapsedNode = transformToTree(nextEdges);
      collapsedNode.addPositions(positions);
      result = new CompositeEdge(label, collapsedNode);
    }
    else {
      result = new CompositeEdge(label, new InternalNode(positions));
    }
    
    // record the successful merge
    ArrayList<TandemEdge> membersCopy = new ArrayList<TandemEdge>(targetEdges);
    MergeInfo mi = new MergeInfo(result, membersCopy);
    for (TandemEdge memberEdge : membersCopy) {
      if (!mergedNodes.containsKey(memberEdge)) {
        mergedNodes.put(memberEdge, new ArrayList<MergeInfo>());
      }
      mergedNodes.get(memberEdge).add(mi);
    }
    
    return result;
  }
   
 
  
  public static void countNodes(Node target, HashSet<Node> seen) {
    
    for (int i=0; i<target.getDegree(); i++) {
      seen.add(target.getEdgeAt(i).getSink());
      countNodes(target.getEdgeAt(i).getSink(), seen);
    }
    
  }
  
  
  
  /**
   * The general method to create the bridging edges that simulate a gap
   * from the current node to its descendant many generations away.
   */
  public void makeGapLinks() {
    
    if (getRoot()==null) return;
    
    PriorityQueue<NodeContainer> q = new PriorityQueue<NodeContainer>();
    HashSet<Node> seenNodes = new HashSet<Node>();
    HashMap<TandemEdge,ArrayList<MergeInfo>> mergedNodes = new HashMap<TandemEdge,ArrayList<MergeInfo>>();
    q.add(new NodeContainer(getRoot(), 0));
    
    TreeMap<TandemEdge,ArrayList<TandemEdge>> children = new TreeMap<TandemEdge,ArrayList<TandemEdge>>();
    int prevLevel = -1;
    while (!q.isEmpty()) {
      NodeContainer nc = q.poll();
      Node currentNode = nc.node;
      
      if (nc.depth!=prevLevel) {
        prevLevel = nc.depth;
        //mergedNodes.clear();
        //System.err.println("--- Clearing the merged edges");
      }
      
      if (seenNodes.contains(currentNode)) {
        continue;
      }
      seenNodes.add(currentNode);
      
      //System.err.println("At level " + nc.depth);
      
      // insert all the children of this node into the queue
      for (int i=0; i<currentNode.getDegree(); i++) {
        Node nextSink = currentNode.getEdgeAt(i).getSink();
        NodeContainer nextNodeContainer = new NodeContainer(nextSink, nc.depth+currentNode.getEdgeAt(i).size());
        //System.err.println("Inserting child " + currentNode.getEdgeAt(i));
        q.add(nextNodeContainer);
      }

      //if (currentNode.getDegree()<=1) continue;
      
     
      // get all reachable descendants nodes from the current node
      children.clear();
      getChildren(currentNode, children);

      // we will need to add each of this edges into the current node
      for (TandemEdge key : children.keySet()) {
        ArrayList<TandemEdge> commonEdges = children.get(key);
        if (commonEdges.size()==1) {
          // easy case, just add the only edge to the current node
          currentNode.insert(key);
        }
        else {
          /*
          System.err.println("Merging");
          for (Edge te : commonEdges) {
            System.err.println(te);
          }*/
          
          // collision, call the recursive merge method to fix the conflict
          CompositeEdge mergedEdge = mergeEdges(commonEdges, mergedNodes);
          NodeContainer nextNodeContainer = new NodeContainer(mergedEdge.getSink(), nc.depth+mergedEdge.length());
          
          q.add(nextNodeContainer);
          
          //System.err.println("Merged edge " + mergedEdge);
          //System.err.println(toString(2, mergedEdge.getSink()));
          
          currentNode.insert(mergedEdge);
        }
      }
    }
    
    
  }
  
  /*
  public boolean verifyEdges(Node root) {
    for (int i=0; i<root.getDegree(); i++) {
      if (!(root.getEdgeAt(i) instanceof CompressedEdge)) return false;
      if (!verifyEdges(root.getEdgeAt(i).getSink())) return false;
    }
    return true;
  }*/
  
  private String toString(int level, Node node) {
    StringBuffer result = new StringBuffer();
    String padding = "";
    for (int i=0; i<level; i++) {
      padding += "  ";
    }
    for (int i=0; i<node.getDegree(); i++) {
      result.append(padding + node.getEdgeAt(i).toString()+"\n");    
      result.append(toString(level+1, node.getEdgeAt(i).getSink()));
    }
    return result.toString();
  }
  
  @Override
  public String toString() {
    return toString(0, getRoot());
  }
  
}
