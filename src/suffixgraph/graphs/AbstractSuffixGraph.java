package suffixgraph.graphs;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import msutil.AminoAcid;
import msutil.Composition;
import msutil.GappedPeptide;
import msutil.Sequence;

import sequences.FastaSequence;
import suffixgraph.Constants;
import suffixgraph.misc.BitArray;
import suffixgraph.misc.LongPriorityQueue;
import suffixgraph.nodes.AbstractNode;
import suffixgraph.nodes.EdgeRuler;
import suffixgraph.nodes.MergedNodeInfo;
import suffixgraph.nodes.Node;



/**
 * The parent class of the composition graphs.
 * @author jung
 *
 */
public abstract class AbstractSuffixGraph<T extends Node> {

  private static long nextTime = 0;    // elapsed time
  
  
/***** MEMBER VARIABLES *****/
  ArrayList<AbstractNode> nodes;
  FastaSequence sequence;
  EdgeRuler er;

  
  /**
   * Read all the nodes from the file initializing the nodes member of this class.
   * This class has the same effect as the buildSuffixGraph class.
   * @param path filename the path of to the file written with the toFile method.
   */
  public abstract void fromFile(String path);

  
  /**
   * Search the this graph given a sequence of amino acids. This function
   * converts the list of amino acids into a list of compositions and queries
   * the graph. Therefore, amino acids with equal compositions will be treated
   * equally.
   * @param query the sequence of amino acids.
   * @return true or false for match or mismatch.
   */
  public ArrayList<String> search(Sequence<? extends AminoAcid> query) {
    ArrayList<Composition> compQuery = new ArrayList<Composition>();
    for (AminoAcid aa : query)    compQuery.add(aa.getComposition());
    return search(compQuery);
  }
  
  
  /**
   * Getter methods for the node at a index.
   * @param index the index of the node to retrieve.
   * @return the retrieved node at the index.
   */
  public AbstractNode getNodeAt(int index) {
    return nodes.get(index);  
  }
  
  /**
   * Search the this graph given a GappedPeptide object.
   * @param query the GappedPeptide object.
   * @return true or false for match or mismatch.
   */
  public ArrayList<String> search(GappedPeptide query) {
    return search(query.getCompositions(), query.getCount());  
  }
  
  
  /**
   * Search the this graph given a sequence of compositions.
   * @param query the sequence of amino acids.
   * @return true or false for match or mismatch.
   */
  public ArrayList<String> search(ArrayList<Composition> query) {
    return search(query, query.size());
  }
  
  
  /**
   * Search the this graph given a sequence of compositions.
   * @param query the sequence of amino acids.
   * @param depth the number of letters to retrieve when reading from the fasta sequence.
   * @return true or false for match or mismatch.
   */
  public abstract ArrayList<String> search(ArrayList<Composition> query, int depth);
  
  
  /**
   * The most general search method given a list of integers.
   * @param query the sequence of amino acids.
   * @param depth the number of letters to retrieve when reading from the fasta sequence.
   * @return true or false for match or mismatch.
   */
  public ArrayList<String> search(int[] query, int depth) {
    AbstractNode currentNode = nodes.get(0);
    for (int edge : query) {
      int matchIndex = currentNode.getEdge(edge);
      if(matchIndex < 0) {
        return null;
      }
      // navigate to the next node
      currentNode = nodes.get(currentNode.getNodeIdAt(matchIndex));
    }
    
    HashSet<Integer> matches = new HashSet<Integer>(); 
    getAllStartingPositions(currentNode, new HashSet<Integer>(), matches); 
    ArrayList<String> strMatches = new ArrayList<String>();
    for (int startPos : matches) {
      strMatches.add(sequence.getSubsequence(startPos, startPos+depth));
    }
    return strMatches;
  }
  
  
  /**
   * Return the number of nodes in this graph.
   * @return the number of nodes in this graph.
   */
  public int size() {
    return nodes.size();
  }
  
  
  /**
   * Get the total number of stating positions stored in the nodes. This is a 
   * good indicator of how redundant the graph is.
   * @return the total number of stating positions stored among all nodes.
   */
  public int getTotalStartPositions() {
    int positionCount = 0;
    for (AbstractNode n : nodes) {
      positionCount += n.getPositions().length;
    }
    return positionCount;
  }

  
  /**
   * Get the average number of outgoing edges per node.
   * @return the average number of outgoing edges per node.
   */
  public float getAverageOutDegree() {
    int edgeCount = 0;
    for (AbstractNode n : nodes) {
      edgeCount += n.getEdgeCount();
    }
    return edgeCount / (float)nodes.size();    
  }
  
  
  
  public String representAsGraphViz() {
    return representAsGraphViz(this.nodes);  
  }
  
  public String representAsGraphViz(ArrayList<? extends AbstractNode> nodes) {
    return null;
  }
  
  
  /**
   * Helper method to visit all the nodes.
   */
  private void retrieve(int nodeId, float mass, int[] nodeCounts, int[] degreeCounts, BitArray seenNodes) {
    AbstractNode n = nodes.get(nodeId);
    nodeCounts[(int)mass]++;
    degreeCounts[(int)mass] = n.getEdgeCount();
    
    for (int i=0; i<n.getEdgeCount(); i++) {
      if (!seenNodes.get(n.getNodeIdAt(i))) {
        seenNodes.set(n.getNodeIdAt(i));
        retrieve(n.getNodeIdAt(i), mass+er.toFloat(n.getKeyAt(i)), nodeCounts, degreeCounts, seenNodes);
      }
    }
  }
  
  
  /**
   * Goes over all the nodes and calculates the average degree per 
   * @return
   */
  public String summarizeNodeCounts() {
    int[] nodeCounts = new int[(int)Constants.MAX_DEPTH_MASS+200];
    int[] degreeCounts = new int[(int)Constants.MAX_DEPTH_MASS+200];
    BitArray seenNodes = new BitArray();
    final int lumpFactor = 20;
    
    retrieve(0, 0, nodeCounts, degreeCounts, seenNodes);
    // recursively visit all the nodes collecting statistics 
    StringBuffer sb = new StringBuffer();
    
    sb.append("Total nodes visited " + seenNodes.getSetItems() + "\n");
    
    for (int i=0; i<nodeCounts.length-lumpFactor; i++) {
      int totalNodes = 0;
      int totalDegrees = 0;
      for (int j=0; j<lumpFactor; j++, i++) {
        totalNodes += nodeCounts[i];
        totalDegrees += degreeCounts[i];
      }
      if (totalNodes>0) {
        sb.append((i-lumpFactor)+"\t"+totalNodes+"\t"+(totalDegrees/(float)totalNodes)+"\n");
      }
    }
    return sb.toString();
  }
  
  
  @Override
  public String toString() {
    long usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
    String retVal = "---- Memory used for the simple graph " + usedMem + "MB at " + nodes.size()/1000 + "K nodes.\n";
    retVal += "----- Average node out degree " + getAverageOutDegree() + '\n';
    retVal += "----- Start position / sequence length ratio " + getTotalStartPositions() / (float)sequence.getSize();
    return retVal;
  }
   
  
  
/***** HELPER METHODS *****/
  /**
   * Given the merged nodes information return a reduced set of merges that
   * are equivalent for the merge. This transformation needs to guaranteed the
   * correctness of the DFA
   * @param targets the ids of the node to merge
   * @param mergedNodes the information of the merged nodes
   * @return a set of ids that are equivalent when merged to the original targets
   */
  private static int[] reduceTargets(int[] targets, ArrayList<ArrayList<MergedNodeInfo>> mergedNodes) {
    
    // speed up the test membership routine
    HashSet<Integer> superset = new HashSet<Integer>();
    for(int target : targets) superset.add(target);
    
    // stores all candidate subsets, maps merged node index -> merging members
    HashMap<Integer,MergedNodeInfo> subsets = new HashMap<Integer,MergedNodeInfo>();
    
    // artificial trackers, to greedily solve the set cover problem
    int largestSetCount = -1, largestSetId = -1;
    
    // convert the array into a hash to optimize set membership query
    HashSet<Integer> targetsHash = new HashSet<Integer>();
    for (int target : targets) targetsHash.add(target);
    
    // iterate over the target nodes to collect possible merges
    for (int nodeId : targets) {
      // iterate over the list of merged items in which this node is member of
      ArrayList<MergedNodeInfo> mergedNodeInfoArray = mergedNodes.get(nodeId);
      if (mergedNodeInfoArray!=null) {
        for (MergedNodeInfo subset : mergedNodeInfoArray) {
          if (!subsets.containsKey(subset.getNodeId()) && subset.isContained(superset)) {
            subsets.put(subset.getNodeId(), subset);
            if (subset.getItems().length > largestSetCount) {
              largestSetCount = subset.getItems().length;
              largestSetId = subset.getNodeId();
            }
          }
        }
      }
    } // all nodes in the mergedSets are subsets of the targets 

    // return value
    ArrayList<Integer> refinedTargets = new ArrayList<Integer>();

    while (superset.size()>0 && largestSetCount>0) {  

      // subtract largest subset from the superset
      for (int subsetItem : subsets.remove(largestSetId).getItems()) {
        superset.remove(subsetItem);
      }
      // the subset is replaced by the merged node id
      refinedTargets.add(largestSetId);
      
      largestSetCount = -1; largestSetId = -1;
      HashMap<Integer,MergedNodeInfo> tempSubsets = new HashMap<Integer,MergedNodeInfo>();
      // make sure all the candidate subsets are still subsets of the reduced superset
      for (MergedNodeInfo subset : subsets.values()) {
        if (subset.isContained(superset)) {
          tempSubsets.put(subset.getNodeId(), subset);
          if (subset.getItems().length > largestSetCount) {
            largestSetCount = subset.getItems().length;
            largestSetId = subset.getNodeId();
          }
        }
      }
      subsets = tempSubsets;
    }
    
    // add the leftovers to the refined set and convert back to an array
    refinedTargets.addAll(superset);
    int[] result = new int[refinedTargets.size()];
    for (int i=0; i < result.length; i++) result[i] = refinedTargets.get(i);
    return result;
  }
  
  
  /**
   * Recursively retrieve the starting positions.
   */
  void getAllStartingPositions(AbstractNode start, HashSet<Integer> visited, HashSet<Integer> results) {
    if (start.getPositions() != null)  
      for (int position : start.getPositions()) 
        results.add(position);
    
    // recursively get all starting positions down the tree
    for (int i = 0; i < start.getEdgeCount(); i++) {
      int nodeId = start.getNodeIdAt(i);
      if (!visited.contains(nodeId)) {
        visited.add(nodeId);
        getAllStartingPositions(nodes.get(nodeId), visited, results);
      }
    }
  }
  
  
  /**
   * Write the sequence of nodes of this graph into the file.
   * @param path the path to the file to write.
   * @param nodes the composition nodes to write.
   */
  static void toFile(String path, ArrayList<Node> nodes) {
    try {
      DataOutputStream nodeFile = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)));
      for (Node n : nodes)        n.toFile(nodeFile);
      System.out.println("--- Number of nodes written to file: " + nodes.size());
      nodeFile.flush(); nodeFile.close();
    } 
    catch(IOException e) {
      System.err.println(e);
    }
  }
  
  
  /**
   * Method that builds the suffix trie with a given maximum depth mass.
   * @param sequence the protein fasta sequence.
   * @return the list of composition nodes built from the fasta sequence. 
   */
  ArrayList<Node> buildSuffixGraph(FastaSequence sequence) {
    return buildSuffixGraph(sequence, (int)sequence.getSize());
  }   
  
  
  /**
   * Method that builds the suffix trie with a given maximum depth mass.
   * @param sequence the protein fasta sequence.
   * @param length the cut off for the number of characters to use for the sequence.
   * @return the list of composition nodes built from the fasta sequence. 
   */
  ArrayList<Node> buildSuffixGraph(FastaSequence sequence, int length) {
    
    ArrayList<Node> nodes = new ArrayList<Node>();
    nodes.add(er.nodeFactory());  // add the root node
    
    // for every start position
    for(int start = 0; start<length; start++) {
      Node currentNode = nodes.get(0);
      Node previousNode = null;
      int prevDistance = 0;
      float cumMass = 0f;
      
      boolean advanced = false;
      // for every amino acid in the database try to extend until the MAX_DEPTH_MASS is reached.
      for(int index=start, depth=0; index<length; depth++, index++) {
        AminoAcid aa = AminoAcid.getStandardAminoAcid(sequence.getCharAt(index));
        if(aa == null) {
          // nothing to do, because most likely a separator.
          break;
        }
        else {          
          advanced = true;
          // print out the debugging info
          //System.out.println("Adding " + aa.toString() + "(" + ec.fromAa(aa) + ") at level " + depth);
          
          // try to add a new edge to the graph
          int thisDistance = er.fromAa(aa);
          int mIndex = currentNode.getEdge(thisDistance);
          if(mIndex<0) {
            // no node with this label found
            int newSinkId = -1;
            
            // there is a possibility that adding another edge, we can find a match
            for (int i=0; i<currentNode.getEdgeCount(); i++) {
              int firstDistance = currentNode.getKeyAt(i);
              
              if (er.compareEdges(thisDistance, firstDistance) < 0)    
                // no need to continue because thisDistance < firstDistance, no need to test second distance  
                break;
              
              Node secondNode = nodes.get(currentNode.getNodeIdAt(i)); 
              for (int j=0; j<secondNode.getEdgeCount(); j++) {
                int compResult = er.compareEdges(thisDistance, firstDistance+secondNode.getKeyAt(j));
                if (compResult == 0) {
                  // we found a matching node
                  newSinkId = secondNode.getNodeIdAt(j);
                }
                else if (compResult < 0) {
                  // no point in following this path
                  break;
                }
              }
            }
            
            // there is also a possibility that this node can be matched to a child of the parent of the current node
            if (newSinkId<0 && previousNode!=null) {
              int prevMatchindex = previousNode.getEdge(thisDistance + prevDistance);
              if (prevMatchindex >= 0) {
                newSinkId = previousNode.getNodeIdAt(prevMatchindex);
              }
            }
            
            if (newSinkId<0) {
              // add thisEdge to the currentNode
              Node newSink = er.nodeFactory();
              currentNode.addEdge(thisDistance, nodes.size(), -mIndex-1);
              previousNode = currentNode; currentNode = newSink;
              nodes.add(newSink);
            }
            else {
              currentNode.addEdge(thisDistance, newSinkId, -mIndex-1);
              previousNode = currentNode; currentNode = nodes.get(newSinkId);
            }
          }
          else {
            // there is already a node, just move on along that path
            previousNode = currentNode; currentNode = nodes.get(currentNode.getNodeIdAt(mIndex));
          }
          cumMass += aa.getMass();
          if(cumMass > Constants.MAX_DEPTH_MASS) {
            break;
          }
          prevDistance = thisDistance;
        }
      }
      // this is the last visited node and also a accepting state
      if (advanced)          currentNode.addPosition(start);
    }
    
    return nodes;
  }
  
  
  /**
   * Get all the nodes that are children of the current node, given the distance
   * restriction.
   * @param target the node to target.
   * @param nodes the node array.
   * @param sum the cumulative composition/mass of the root to the current root node.
   * @param maxMass the maximum distance (mass) allowed.
   * @param values a map composition -> unique set of nodes reachable within the
   *               specified distance. This is the return value.
   * @param er the edge ruler to allow this method to translate edges to distances.              
   */
  private static void getChildNodes(Node target, ArrayList<Node> nodes, int sum,
                                    float maxMass, HashMap<Integer,HashSet<Integer>> values, EdgeRuler er) {
    for (int i=0; i<target.getEdgeCount(); i++) {
      int cumulativedMass = sum + target.getKeyAt(i);
      if (er.toFloat(cumulativedMass) < maxMass) {
        int childNodeId = target.getNodeIdAt(i);
        if (!values.containsKey(cumulativedMass)) {
          HashSet<Integer> nodeIds = new HashSet<Integer>();
          nodeIds.add(childNodeId);
          values.put(cumulativedMass, nodeIds);
        }
        else {
          values.get(cumulativedMass).add(childNodeId);
        }
        getChildNodes(nodes.get(childNodeId), nodes, cumulativedMass, maxMass, values, er);
      }
    }
  }
  
  
  /**
   * This routine modifies the graph so that all the gaps are encoded while
   * conserving the DFA properties using a top down recursion.
   * The target must have unique outgoing edge labels.
   * @param nodes the array of nodes.
   * @param seenNodes the bit array marking which nodes have been processed.
   * @param er the object that translates edges to distances.
   */
  void makeGapLinks(ArrayList<Node> nodes, ArrayList<ArrayList<MergedNodeInfo>> mergedNodes, BitArray seenNodes) {
    
    // upper half represents the mass, lower half represents the node index.
    LongPriorityQueue q = new LongPriorityQueue();
    q.add(0L);
    
    HashMap<Integer,HashSet<Integer>> children = new HashMap<Integer,HashSet<Integer>>();
    while (!q.isEmpty()) {
      long current = q.poll();
      Node currentNode = nodes.get((int)current);
      
      // Gather all the child nodes
      children.clear();
      getChildNodes(currentNode, nodes, 0, Constants.MAX_GAP_MASS, children, er);
      
      for (int comp : children.keySet()) {
        int matchIndex = currentNode.getEdge(comp);
        HashSet<Integer> destNodes = children.get(comp);
        // remove existing direct links to the children node
        if (matchIndex>=0) destNodes.remove(currentNode.getNodeIdAt(matchIndex));
        
        if (destNodes.size()==1) {
          if (matchIndex<0) {
            // easy case, just add a link with no new nodes
            currentNode.addEdge(comp, destNodes.iterator().next(), -matchIndex -1);
            destNodes.clear();
          }
          else {
            // in place merging   
            int[] destNodeArray = new int[1];
            destNodeArray[0] = destNodes.iterator().next();
            mergeNodes(currentNode.getNodeIdAt(matchIndex), destNodeArray, nodes, mergedNodes);
          }
        }
        else if (destNodes.size()>1) {
          // Merge the matches
          if (matchIndex >= 0) {
            int[] destNodeArray = new int[destNodes.size()];
            int j = 0;
            for (int nodeIndex : destNodes) destNodeArray[j++] = nodeIndex;
            destNodes.clear();
            mergeNodes(currentNode.getNodeIdAt(matchIndex), destNodeArray, nodes, mergedNodes);
          }
          else {
            int[] destNodeArray = new int[destNodes.size()];
            int j = 0;
            for (int nodeIndex : destNodes) destNodeArray[j++] = nodeIndex;
            currentNode.addEdge(comp, mergeNodes(destNodeArray, nodes, mergedNodes), -matchIndex-1);
          }
        }
      }
      
      // current node is transformed, add its children to the queue
      for (int i=0; i < currentNode.getEdgeCount(); i++) {
        int nextNodeId = currentNode.getNodeIdAt(i);
        if (!seenNodes.get(nextNodeId)) {
          seenNodes.set(nextNodeId);
          int transition = (int)er.toFloat(currentNode.getKeyAt(i));
          long next = (((long)((current>>>32)+transition))<<32) | nextNodeId;
          q.add(next);
        }
      }  
      
      if (nextTime < System.currentTimeMillis()) {
        nextTime = System.currentTimeMillis()+5000;
        long usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
        System.out.print("----- Used " + usedMem + "MB at " + seenNodes.getSetItems()/1000 + "K nodes at mass " + (current>>>32) + "Da. ");
        System.out.printf("~%.2f%% done. Queue Size %dK.\n", 100.0*seenNodes.getSetItems()/nodes.size(), q.size()/1000);
      }
      
    // clear the stored merges for this node because it will not be queried anymore
    mergedNodes.set((int)current, null);
    }
  }
  
  
  /**
   * In site merge replacing the nId with the resulting merged node.
   * @param nId the id of the node to merge
   * @param targets the other nodes to merge
   * @param nodes the array of all nodes
   * @param mergedNodes keeps tracks of what nodes have been merged to avoid repeating merges
   */
  private void mergeNodes(int nId, int[] targets, ArrayList<Node> nodes, ArrayList<ArrayList<MergedNodeInfo>> mergedNodes) {
    
    //targets = reduceTargets(targets, mergedNodes);
    
    Node n = nodes.get(nId);
    
    // the list of nodes to merge
    Node[] mergingNodes = new Node[targets.length];
    int index = 0;
    for (int nodeIndex : targets) {
      mergingNodes[index++] = nodes.get(nodeIndex);
    }
    
    // holds the list of nodes with the same composition to merge
    HashSet<Integer> mergingTargetHash = new HashSet<Integer>();
    
    // the indices of the current pointer to the target edges. Initialized to 0.
    int[] currentIndeces = new int[targets.length];
    
    long prevEdge = 0;  // important to keep track of the last seen different composition
    long minEdge = 0L;  // registers the last edge added
    while (true) {
      // find the minimum item
      int minIndex = -1;
      for (int i=0; i<currentIndeces.length; i++) {
        if (currentIndeces[i]<mergingNodes[i].getEdgeCount()) {
          long currentEdge = mergingNodes[i].getEdgeAt(currentIndeces[i]);
          if (minIndex<0 || er.compareEdges((int)currentEdge, (int)minEdge)<0) {
            minIndex = i; minEdge = currentEdge;
          }
        }
      }

      // nothing more to do because all currentIndeces have reached their maximum
      if (minIndex<0)   break;
      
      // two cases, this is a new  composition, or repeating one
      if (mergingTargetHash.size()>0) {          // prevComp has been initialized  
        if (((int)minEdge)!=(int)prevEdge) {     // new composition observed, process all the mTargets
          
          // convert the current merging targets into an array of integers
          int nMatchIndex = n.getEdge((int)prevEdge);
          if (nMatchIndex>=0) mergingTargetHash.remove(n.getNodeIdAt(nMatchIndex));
          int[] mergingTargetArray = new int[mergingTargetHash.size()];
          int j = 0;
          for (int targetId : mergingTargetHash) mergingTargetArray[j++] = targetId;
          mergingTargetHash.clear();
          
          if (mergingTargetArray.length==1) {
            // easy case, we can just add the new edge to the node
            if (nMatchIndex<0) n.addEdge(prevEdge, -nMatchIndex-1);
            // recursively merge the destination conflicting nodes
            else mergeNodes(n.getNodeIdAt(nMatchIndex), mergingTargetArray, nodes, mergedNodes);
          }
          else if (mergingTargetArray.length>1){ // we will have to recurse because of conflicting edges
            if (nMatchIndex<0) n.addEdge((int)prevEdge, mergeNodes(mergingTargetArray, nodes, mergedNodes), -nMatchIndex-1);
            else mergeNodes(n.getNodeIdAt(nMatchIndex), mergingTargetArray, nodes, mergedNodes);
          }
          prevEdge = minEdge;
        }

      }
      else // record the index of the last node, nothing to compare it against
        prevEdge = minEdge;
      
      mergingTargetHash.add((int)(minEdge>>>32));    // record this nodeId  
      currentIndeces[minIndex]++;           // advance the pointer of the minimum element
    }
   
    // process the last batch of mTargets
    if (mergingTargetHash.size()>0) {
      // convert the current merging targets into an array of integers
      int nMatchIndex = n.getEdge((int)prevEdge);
      if (nMatchIndex>=0) mergingTargetHash.remove(n.getNodeIdAt(nMatchIndex));
      int[] mergingTargetArray = new int[mergingTargetHash.size()];
      int j = 0;
      for (int targetId : mergingTargetHash) mergingTargetArray[j++] = targetId;
      mergingTargetHash.clear();
      
      if (mergingTargetArray.length==1) {        
        // easy case, just add this edge to the finals
        if (nMatchIndex<0) n.addEdge(prevEdge, -nMatchIndex-1);
        // merging needs to be done
        else mergeNodes(n.getNodeIdAt(nMatchIndex), mergingTargetArray, nodes, mergedNodes);
      }
      else if (mergingTargetArray.length>1){     // we will have to recurse because of conflicting edges
        if (nMatchIndex<0) n.addEdge((int)prevEdge, mergeNodes(mergingTargetArray, nodes, mergedNodes), -nMatchIndex-1);
        else mergeNodes(n.getNodeIdAt(nMatchIndex), mergingTargetArray, nodes, mergedNodes);
      }
    }
    mergingTargetHash = null;
    
    // update the start positions to an array
    HashSet<Integer> startPos = new HashSet<Integer>();
    for (int position : n.getPositions())
      startPos.add(position);
    for (int target : targets) 
      for (int position : nodes.get(target).getPositions()) 
        startPos.add(position);
    int[] positionArray = new int[startPos.size()];
    Iterator<Integer> positionIterator = startPos.iterator();
    for (int i=0; i<positionArray.length; i++) 
      positionArray[i] = positionIterator.next();
    
    n.setPositions(positionArray);
    
    // create a new mergedInfo object and register it with the mergedNodes array
    /*
    MergedNodeInfo thisMerge = new MergedNodeInfo(nId, allTargets);
    for (int memberId : allTargets) {
      ArrayList<MergedNodeInfo> merges = mergedNodes.get(memberId);
      if (merges==null) {
        mergedNodes.set(memberId, new ArrayList<MergedNodeInfo>());
      }
      mergedNodes.get(memberId).add(thisMerge);
    }*/
  }
  
  
  /**
   * Merge all nodes specified by targets and return the id of the resulting 
   * node.
   * @param nId the id of the node to merge
   * @param targets the other nodes to merge
   * @param nodes the array of all nodes
   * @param er the object for comparison of edges.
   * @return the id of the newly merged node.
   */
  private int mergeNodes(int[] targets, ArrayList<Node> nodes, ArrayList<ArrayList<MergedNodeInfo>> mergedNodes) {
    
    // reduce the set of merged nodes by reusing previous merges
    targets = reduceTargets(targets, mergedNodes);
    
    // create a new node
    int mergedNodeId = nodes.size();
    nodes.add(null);
    mergedNodes.add(null);
    
    // the list of nodes to merge
    Node[] mergingNodes = new Node[targets.length];
    int index = 0;
    for (int nodeIndex : targets) {
      mergingNodes[index++] = nodes.get(nodeIndex);
    }

    // holds the list of nodes with the same composition to merge
    HashSet<Integer> mergingTargetHash = new HashSet<Integer>();
    
    // final list of edges coming out of the merge
    ArrayList<Long> outEdges = new ArrayList<Long>();
    
    // the indices of the current pointer to the target edges. Initialized to 0.
    int[] currentIndeces = new int[mergingNodes.length];
    
    long prevEdge = 0;  // important to keep track of the last seen different composition
    long minEdge = 0L;  // registers the last edge added
    while (true) {
      // find the minimum item
      int minIndex = -1;
      for (int i=0; i<currentIndeces.length; i++) {
        if (currentIndeces[i]<mergingNodes[i].getEdgeCount()) {
          long currentEdge = mergingNodes[i].getEdgeAt(currentIndeces[i]);
          if (minIndex<0 || er.compareEdges((int)currentEdge, (int)minEdge)<0) {
            minIndex = i; minEdge = currentEdge;
          }
        }
      }

      // nothing more to do because all currentIndeces have reached their maximum
      if (minIndex<0)   break;
      
      // two cases, this is a new  composition, or repeating one
      if (mergingTargetHash.size()>0) {          // prevComp has been initialized  
        if (((int)minEdge) != (int)prevEdge) {   // new composition observed, process all the mTargets
          // the current merging targets into an array of integers
          int[] mergingTargetArray = new int[mergingTargetHash.size()];
          int j = 0;
          for (int targetId : mergingTargetHash) mergingTargetArray[j++] = targetId;
          mergingTargetHash.clear();
          
          // easy case, just add this edge to the finals
          if (mergingTargetArray.length==1) outEdges.add(prevEdge);
          // we will have to recurse because of conflicting edges
          else outEdges.add((((long)mergeNodes(mergingTargetArray, nodes, mergedNodes))<<32)|((int)prevEdge));
          prevEdge = minEdge;
        }
      }
      // record the index of the last node, nothing to compare it against
      else prevEdge = minEdge;
      
      mergingTargetHash.add((int)(minEdge>>>32));// record this nodeId  
      currentIndeces[minIndex]++;                // advance the pointer of the minimum element
    }
   
    // process the last batch of mTargets
    if (mergingTargetHash.size()>0) {
      int[] mergingTargetArray = new int[mergingTargetHash.size()];
      int j = 0;
      for (int targetId : mergingTargetHash) mergingTargetArray[j++] = targetId;
      mergingTargetHash.clear();
      
      // easy case, just add this edge to the finals
      if (mergingTargetArray.length==1) outEdges.add(prevEdge);
      // we will have to recurse because of conflicting edges
      else outEdges.add((((long)mergeNodes(mergingTargetArray, nodes, mergedNodes))<<32)|((int)prevEdge));
    }
    
    // convert the sorted edges from an ArrayList to an array
    long[] sortedEdges = new long[outEdges.size()];
    for (int i=0; i<outEdges.size(); i++) sortedEdges[i] = outEdges.get(i);
    
    // convert the start positions to an array
    HashSet<Integer> startPos = new HashSet<Integer>();
    for (int target : targets) 
      for (int position : nodes.get(target).getPositions()) 
        startPos.add(position);
    int[] positionArray = new int[startPos.size()];
    Iterator<Integer> positionIterator = startPos.iterator();
    for (int i=0; i<positionArray.length; i++) 
      positionArray[i] = positionIterator.next();
      
    Node retNode = er.nodeFactory(sortedEdges, outEdges.size(), positionArray);
    nodes.set(mergedNodeId, retNode);
    
    // create a new mergedInfo object and register it with the mergedNodes array
    MergedNodeInfo thisMerge = new MergedNodeInfo(mergedNodeId, targets);
    for (int memberId : targets) {
      ArrayList<MergedNodeInfo> merges = mergedNodes.get(memberId);
      if (merges==null) {
        mergedNodes.set(memberId, new ArrayList<MergedNodeInfo>());
      }
      mergedNodes.get(memberId).add(thisMerge);
    }
    
    return mergedNodeId;
  }
}
