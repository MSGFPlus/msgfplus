package suffixtree.trees.deprecated;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.TreeMap;

import sequences.FastaSequence;
import sequences.ProteinFastaSequence;
import suffixtree.Constants;
import suffixtree.edges.Edge;
import suffixtree.matches.ExactMatchObject;
import suffixtree.nodes.InternalNode;
import suffixtree.nodes.Node;
import suffixtree.trees.IntegerSuffixTree;


/**
 * This is the gapped version of the suffix tree.
 * @author jung
 *
 */
public class IntegerGappedSuffixTree extends IntegerSuffixTree {

  /**
   * An tandem edge is copy edge of another compressed edge with an implicit node
   * splitting the given node.
   * @author jung
   *
   */
  class TandemEdge extends Edge {

    private int mass;        // the label of the first half
    private int end;         // the index of the end of the second edge
    private int start;       // the index of the start of the second edge
    
    /**
     * Constructor
     * @param sink the destination sink
     * @param mass the mass of the first edge
     * @param start the start position of the second edge
     * @param end the end position of the second edge
     */
    public TandemEdge(Node sink, int mass, int start, int end) {
      this.mass = mass;
      this.start = start;
      this.end = end;
      setSink(sink);
    }
    
    /**
     * Get the actual length of the label.
     * @return the number of letters (bytes) that this label represents.
     */
    public int length() {
      return 1;
    }
    
    @Override
    public int getLabelAt(int offset) {
      if (offset > 0)
        System.err.println("Attempting to retrieve a label that does not exist... returning first label.");
      return mass;
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
      String result = "TandemEdge: [" + this.mass + "] -o- ";
      result += sequence.getSubsequence(start, Math.min(start+6, end));
      result += "["+start+","+end+") ---> " + getSink().hashCode();
      return result;
    }
    
    @Override
    public int hashCode() {
      return getSink().hashCode();
    }
    
    @Override
    public boolean equals(Object other) {
      TandemEdge o = (TandemEdge)other;
      return this.getSink()==o.getSink() && 
             this.start==o.start && 
             this.end==o.end;
    }

    @Override
    public int mass() {
      System.err.println("Unsupported method: IntegerGappedSuffixTree.TandemEdge.mass()");
      return 0;
    }
  }  
  
  /**
   * An IntegerEdge represents the sum of multiple amino acid masses.
   * @author jung
   *
   */
  class IntegerEdge extends Edge {
    
    // this edge can encode a label of 4 bytes
    private int mass;

    public IntegerEdge(int mass, Node sink) {
      this.mass = mass;
      setSink(sink);
    } 
    
    @Override
    public int getLabelAt(int offset) {
      if (offset > 0)
        System.err.println("Attempting to retrieve a label that does not exist... returning first label.");
      return mass;
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
      return mass;
    }

    @Override
    public Edge split(int offset) {
      System.err.println(toString() + " does not support splitting at " + offset);
      System.exit(-1);
      return null;
    }

    @Override
    public int getEnd() {
      System.err.println("Unsupported method: get end position of IntegerEdge");
      System.exit(-1);
      return 0;
    }
    
    @Override
    public int getStart() {
      System.err.println("Unsupported method: get start position of IntegerEdge");
      System.exit(-1);
      return 0;
    }
    
    @Override
    public String toString() {
      return "IntegerEdge: " + this.mass + " Da";
    }

    @Override
    public int mass() {
      return mass;
    }
  }
  
  
  /**
   * This is the class to hold the information for merged items.
   * @author jung
   *
   */
  private static class MergeInfo {
    private IntegerEdge mergedEdge;
    private ArrayList<TandemEdge> members;
    
    public MergeInfo(IntegerEdge mergedEdge, ArrayList<TandemEdge> members) {
      this.mergedEdge = mergedEdge;
      this.members = members;
    }
    
    public boolean isContained(HashSet<TandemEdge> others) {
      for (TandemEdge te : members) if(!others.contains(te)) return false;
      return true;
    }
  }
  
  
  
  /**
   * Constructor taking the ProteinFastaSequence as the parameter. 
   * @param sequence the sequence to build the tree on.
   */
  public IntegerGappedSuffixTree(ProteinFastaSequence sequence) {
    super(sequence);
    System.out.println("Loaded sequence into memory");
    System.out.println(this.collectStats());
    makeGapLinks();
  }
  
  /**
   * Constructor that initializes the tree but does not insert any sequences.
   * @param sequence
   * @param noInsert
   */
  public IntegerGappedSuffixTree(ProteinFastaSequence sequence, boolean noInsert) {
    super(sequence, noInsert);
  }
  
  /**
   * Helper method and entry point to collect all descent of the given node.
   * Note that gap edges (tandem edges) have size of at least 2.
   * @param n the node to collect the children from. n must only have simple
   *          outgoing edges (single byte labels).
   * @param children the data structure to store the results.
   */
  private void getChildren(Node n, TreeMap<TandemEdge,HashSet<TandemEdge>> children) {
    
    // collect all possible children that are within the gap cut off
    for (int i=0; i<n.getDegree(); i++) {
      // the edge is a compressed edge
      Edge e = n.getEdgeAt(i);

      // cumulative mass so far
      int cumMass = 0;
      
      // Tandem edges need have the first edge of at least 2
      cumMass = e.getLabel();
      
      // add all the tandem edges representing gaps to the children structure
      for (int j=1; j<e.size(); j++) {
        
        cumMass += e.getLabelAt(j);
        if (cumMass>Constants.MAX_GAP_MASS) break;
        
        int start = e.getStart() + j + 1;
        
        // copy the prefix into a new array to create the tandem edge
        TandemEdge te = new TandemEdge(e.getSink(), cumMass, start, e.getEnd());
        
        if (children.containsKey(te)) {
          children.get(te).add(te);
        }
        else {
          HashSet<TandemEdge> edgeList = new HashSet<TandemEdge>();
          edgeList.add(te);
          children.put(te, edgeList);
        }
      }
      
      if (cumMass<Constants.MAX_GAP_MASS) {
        // need to recurse and keep collecting children
        getChildren(e.getSink(), children, cumMass);
      }
    }
  }

  /**
   * Helper class that collects the statistics of the tree.
   * @author jung
   *
   */
  private class TreeStats {
    private HashSet<Node> iNodes;
    private HashSet<Edge> tEdges;
    private HashSet<Edge> iEdges;
    private HashSet<Node> lNodes;
    private HashSet<Edge> cEdges;
    
    TreeStats() {
      iNodes = new HashSet<Node>();
      lNodes = new HashSet<Node>();
      iEdges = new HashSet<Edge>();
      cEdges = new HashSet<Edge>();
      tEdges = new HashSet<Edge>();
    }
    
    @Override
    public String toString() {
      StringBuffer sb = new StringBuffer();
      sb.append("Internal nodes " + iNodes.size() + "\n");
      sb.append("Leaf nodes " + lNodes.size() + "\n");
      sb.append("Integer edges " + iEdges.size() + "\n");
      sb.append("Tandem edges " + tEdges.size() + "\n");
      sb.append("Compressed edges " + cEdges.size() + "\n");
      return sb.toString();
    }
  }

  /**
   * Collect the statistics of this data structure. Node and edges.
   * @return the statistics.
   */
  public TreeStats collectStats() {
    TreeStats ts = new TreeStats();
    ts.iNodes.add(getRoot());
    //System.out.println("The root is " + getRoot());
    collectStats(getRoot(), new HashSet<Node>(), ts);
    return ts;
  }
  
  /**
   * Collect statistics about this graph.
   * @param n the node to operate on
   * @param seenNodes what we have seen already
   * @param stats the statistics object
   */
  private void collectStats(Node n, HashSet<Node> seenNodes, TreeStats stats) {
    if (n.getPositions().length == 0) {
      stats.iNodes.add(n);
    }
    else {
      /*
      if (n.getPositions().length > 2) {
      for (int start : n.getPositions()) {
        System.out.println(this.getSequence().toString(start, start+20));
      }
      System.out.println();
      }
      */
      stats.lNodes.add(n);
    }
    
    for (int index = 0; index < n.getDegree(); index++) {
      
      Node child = n.getEdgeAt(index).getSink();
      if (n.getEdgeAt(index) instanceof TandemEdge) {
        stats.tEdges.add(n.getEdgeAt(index)); 
      }
      else if (n.getEdgeAt(index) instanceof IntegerEdge) {
        stats.iEdges.add(n.getEdgeAt(index));
      }
      else {
        stats.cEdges.add(n.getEdgeAt(index));
      }
      
      if (!seenNodes.contains(child)) {
        seenNodes.add(child);
        collectStats(child, seenNodes, stats);
      }
    }
  }
  
  
  /**
   * Helper method to get all descendant of a node within the gap cutoff.
   * @param n the node to operate
   * @param children the data structure to store the results
   * @param prefix the bytes seen so far
   */
  private void getChildren(Node n, TreeMap<TandemEdge,HashSet<TandemEdge>> children, int cumMass) {
    
    // collect all possible children that are within the gap cut off
    for (int i=0; i<n.getDegree(); i++) {
      Edge e = n.getEdgeAt(i);
      
      int thisCumMass = cumMass;
      // add all the tandem edges representing gaps to the children structure
      for (int j=0; j<e.size(); j++) {
        
        thisCumMass += e.getLabelAt(j);
        if (thisCumMass>Constants.MAX_GAP_MASS) break;
        
        TandemEdge te = new TandemEdge(e.getSink(), thisCumMass, e.getStart()+j+1, e.getEnd());
        
        if (children.containsKey(te)) {
          children.get(te).add(te);
        }
        else {
          HashSet<TandemEdge> edgeList = new HashSet<TandemEdge>();
          edgeList.add(te);
          children.put(te, edgeList);
        }
      }
      
      if (thisCumMass<Constants.MAX_GAP_MASS) {
        // need to recurse and keep collecting children
        getChildren(e.getSink(), children, thisCumMass);
      }
    }
  }
  
  /**
   * Helper class to be able to wrap the Node into a priority queue sorted by 
   * their depth 
   */
  private static class NodeContainer implements Comparable<NodeContainer> {
    private int cumMass;
    private Node node;
    //private String path; // the path to this node container
    public NodeContainer(Node node, int cumMass) {
      this.node = node; this.cumMass = cumMass;
    }
    @Override
    public int compareTo(NodeContainer other) {
      return this.cumMass-other.cumMass;
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
  
  
  
  public void search(ArrayList<Integer> query, HashSet<ExactMatchObject> results) {
    
    Node currentNode = getRoot();
    if (currentNode==null) return;
    Integer queryEdge = query.get(0);
    int queryIndex = 0;
 
    while (true) {
      
      // we have reached the end of the tree with no match
      if (currentNode.getDegree()==0) return;
      
      //System.err.println("Matching: " + queryIndex + " - " + queryEdge.getLabel());
      
      int matchIndex = currentNode.search(queryEdge);
      int offset = 0;
      if (matchIndex < 0) {
        // no match, and the outgoing edge is not single
        if (currentNode.getDegree()>1) return;
        
        matchIndex = 0;
        offset = -1;
      }
      
      Edge matchingEdge = currentNode.getEdgeAt(matchIndex);

      //System.err.println("Matched " + matchingEdge);
      
      // we matched the current edge, try to advance
      currentNode = matchingEdge.getSink();
      if (offset==0) {
        // current edge of the query has been matched try to match the next edge.
        queryIndex++;
        if (queryIndex>=query.size()) {
          // matched completely the query, we got a match
          break;     // found matching node
        }
        queryEdge = query.get(queryIndex);       // update next query edge
      }
      
      if (matchingEdge instanceof IntegerEdge) {
        // simply move on to the next internal node and continue matching
        continue;
      }
      
      // the matching edge is either a tandem edge or compressed edge
      int start = 0, end = 0;
      if (matchingEdge instanceof TandemEdge) {
        TandemEdge te = (TandemEdge)matchingEdge;
        start = te.start;
        end = te.end;
      }
      else {
        CompressedEdge ce = (CompressedEdge)matchingEdge;
        start = ce.getStart() + 1 + offset;      // the first one was already matched
        end = ce.getEnd();
        //System.err.println(start + ":" + end);
      }

      if (start==end) {
        continue;
      }
      
      // need to match the unique edge with the query edges greedy
      //System.err.println("Matching " + e.getLabel() + " with " + getSequence().toString(start, end));
      
      int currentIndex = start;
      int cumMass = 0;
      boolean done = false;
      while (currentIndex < end) {
        cumMass += getSequence().getIntegerMass(currentIndex);
        currentIndex++;
        //System.err.println("Matching " + cumMass + " with " + e.getLabel());
        if (cumMass == queryEdge) {
          queryIndex++;
          if (queryIndex>=query.size()) {
            done = true;
            break;
          }
          cumMass = 0;
          queryEdge = query.get(queryIndex);
        }
        
        if (cumMass > queryEdge) {
          // no match
          return;
        }
      }
   
      
      if (done) break;
      
      if (start==end) {
        // move on to the next internal node
        continue;
      }
      
      //System.err.println(e.getLabel());
      // match and remove the leftovers
      int unmatchedMass = queryEdge - cumMass;
      if (unmatchedMass<=0) {
        // no match
        //System.err.println("Removing " + getSequence().toString(start, end));
        //System.err.println("Rejected because no match");
        return;
      }
      //System.err.println("Mass left: " + unmatchedMass);
      queryEdge = unmatchedMass;
    }
    
    //System.err.println("Matching node is " + currentNode);
    // We found the matching node, return starting positions
    //System.out.println("Retrieving results on node with degree " + currentNode.getDegree());
    HashSet<Integer> startPos = new HashSet<Integer>();
    currentNode.getAllPositions(startPos);
    int totalMass = 0;
    for (int e : query) totalMass += e;
    for (int start : startPos) {
      System.out.println("Fail attempt to create MatchObject at position " + start);  
      //results.add(new ExactMatchObject(getSequence(), start, query, totalMass, 0));
    }
    return;
    
  }
  
  
  
  public static ArrayList<IntegerEdge> reduceEdges(Collection<TandemEdge> edges, 
                                                   HashMap<TandemEdge,ArrayList<MergeInfo>> mergedNodes,
                                                   ArrayList<TandemEdge> remainingEdges) {
    HashSet<TandemEdge> targets = new HashSet<TandemEdge>();
    targets.addAll(edges);
    
    ArrayList<IntegerEdge> results = new ArrayList<IntegerEdge>();
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
                if (chosenCandidate.members.size() < candidate.members.size()) {
                  chosenCandidate = candidate;
                }
              }
            }
          }
        }
      }
      // no more seen merges
      if (chosenCandidate==null) break;
      
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
  private IntegerEdge mergeEdges(Collection<TandemEdge> targetEdges, 
                                 HashMap<TandemEdge,ArrayList<MergeInfo>> mergedNodes) {
    
    //mergedNodes.clear();
    
    // nextEdges are actually all compressed edges with a collided parent edge
    int label = targetEdges.iterator().next().getLabel();
    
    ArrayList<TandemEdge> remainingEdges = new ArrayList<TandemEdge>();
    ArrayList<IntegerEdge> reusedEdges = reduceEdges(targetEdges, mergedNodes, remainingEdges);
    
    if (remainingEdges.size()==0) {
      if (reusedEdges.size()==1) { 
        return new IntegerEdge(label, reusedEdges.get(0).getSink());
      }
    }
    
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
    for (IntegerEdge ie : reusedEdges) {
      Node sink = ie.getSink();
      for (int i=0; i<sink.getDegree(); i++) {
        nextEdges.add(sink.getEdgeAt(i));
      }
      for (int pos : sink.getPositions()) {
        positions.add(pos);
      }
    }
    
    IntegerEdge result = null;
    if (nextEdges.size()>0) {
      InternalNode collapsedNode = transformToTree(nextEdges);
      collapsedNode.addPositions(positions);
      result = new IntegerEdge(label, collapsedNode);
    }
    else {
      result = new IntegerEdge(label, new InternalNode(positions));
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
    
    /*
    HashSet<Integer> postPositions = new HashSet<Integer>();
    result.getSink().getAllPositions(postPositions);
    if (prePositions.size()!=postPositions.size()) {
      for (int l : prePositions) System.err.print(l+" ");
      System.err.println();
      for (int l : postPositions) System.err.print(l+" ");
      System.err.println();
      
      System.err.println("Merging... ");
      for (TandemEdge te1 : targetEdges) {
        System.err.println(te1);
        System.err.println(toString(2, te1.getSink()));
      }
      
      System.err.println("Merged edge " + result);
      System.err.println(toString(2, result.getSink()));
      
      System.exit(-1);
    }
    */
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

    long lastTime = System.currentTimeMillis();
    
    PriorityQueue<NodeContainer> q = new PriorityQueue<NodeContainer>();
    HashSet<Node> seenNodes = new HashSet<Node>();
    HashMap<TandemEdge,ArrayList<MergeInfo>> mergedNodes = new HashMap<TandemEdge,ArrayList<MergeInfo>>();
    q.add(new NodeContainer(getRoot(), 0));
    
    TreeMap<TandemEdge,HashSet<TandemEdge>> children = new TreeMap<TandemEdge,HashSet<TandemEdge>>();
    while (!q.isEmpty()) {
      NodeContainer nc = q.poll();
      Node currentNode = nc.node;
      
      /*
      if (!verifyEdges(currentNode)) {
        System.err.println("This node has been processed before!");
        System.exit(-1);
      }*/
      
      if (currentNode.getDegree()==1) {
        Edge e = currentNode.getEdgeAt(0);
        NodeContainer nextNodeContainer = new NodeContainer(e.getSink(), nc.cumMass+e.mass());
        q.add(nextNodeContainer);
        continue;
      }
      
      if (seenNodes.contains(currentNode)) { 
        continue;
      }
      seenNodes.add(currentNode);

      if (System.currentTimeMillis() - lastTime > 20000) {
        lastTime = System.currentTimeMillis();
        System.out.println("At mass " + nc.cumMass);
      }
      
      // insert all the children of this node into the queue
      for (int i=0; i<currentNode.getDegree(); i++) {
        Node nextSink = currentNode.getEdgeAt(i).getSink();
        NodeContainer nextNodeContainer = new NodeContainer(nextSink, nc.cumMass+currentNode.getEdgeAt(i).mass());
        q.add(nextNodeContainer);
      }
     
      // get all reachable descendants nodes from the current node
      children.clear();
      getChildren(currentNode, children);

      // we will need to add each of this edges into the current node
      for (TandemEdge key : children.keySet()) {
        HashSet<TandemEdge> commonEdges = children.get(key);
        int matchIndex = currentNode.search(key);
        if (commonEdges.size()==1) {
          // There is a no matching edge in the currentNode
          if (matchIndex < 0) {
            // easy case, just add the only edge to the current node
            currentNode.insert(key, -matchIndex-1);
          }
          else {
            // add the matching set of common edges
            Edge me = currentNode.getEdgeAt(matchIndex);
            int start = me.getStart()+1, end = me.getEnd();
            TandemEdge additionalEdge = new TandemEdge(me.getSink(), key.mass, start, end);
            commonEdges.add(additionalEdge);
            IntegerEdge mergedEdge = mergeEdges(commonEdges, mergedNodes);
            
            // replace the matching edge with the new merge
            currentNode.setEdgeAt(matchIndex, mergedEdge);

            // record the make gap link
            NodeContainer nextNodeContainer = new NodeContainer(mergedEdge.getSink(), nc.cumMass+key.mass);
            q.add(nextNodeContainer);
          }
          
        }
        else {
          //System.err.println("Merging");
          //for (Edge te : commonEdges) System.err.println(te);
          
          // There is a matching edge in the currentNode
          if (matchIndex < 0) {
            IntegerEdge mergedEdge = mergeEdges(commonEdges, mergedNodes);
            currentNode.insert(mergedEdge, -matchIndex-1);
            NodeContainer nextNodeContainer = new NodeContainer(mergedEdge.getSink(), nc.cumMass+key.mass);
            q.add(nextNodeContainer);
          }
          else {
            // add the matching edge to the common edges
            Edge me = currentNode.getEdgeAt(matchIndex);
            int start = me.getStart()+1, end = me.getEnd();
            TandemEdge additionalEdge = new TandemEdge(me.getSink(), key.mass, start, end);
            commonEdges.add(additionalEdge);
            IntegerEdge mergedEdge = mergeEdges(commonEdges, mergedNodes);
            
            // replace the matching edge with the new merge
            currentNode.setEdgeAt(matchIndex, mergedEdge);

            // record the new node to make more gap links
            NodeContainer nextNodeContainer = new NodeContainer(mergedEdge.getSink(), nc.cumMass+key.mass);
            q.add(nextNodeContainer);
          }
        }
        
        //Scanner scanner = new Scanner(System.in);
        //scanner.next();
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
  }
  */
  
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
