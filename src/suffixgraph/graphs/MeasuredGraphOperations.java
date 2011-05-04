package suffixgraph.graphs;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import msutil.AminoAcid;

import sequences.FastaSequence;
import suffixgraph.Constants;
import suffixgraph.misc.BitArray;
import suffixgraph.misc.GraphvizGraph;
import suffixgraph.misc.LongPriorityQueue;
import suffixgraph.nodes.EdgeRuler;
import suffixgraph.nodes.Node;


/**
 * Provides a re-implementation of the graph operations, but with statistics 
 * printed to a file.
 * @author jung
 *
 */
public class MeasuredGraphOperations {
  
  private static long nextTime = 0;
  
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
   * Method that builds the suffix trie with a given maximum depth mass.
   * @param sequence the protein fasta sequence.
   * @param length the cut off for the number of characters to use for the sequence.
   * @return the list of composition nodes built from the fasta sequence. 
   */
  static ArrayList<Node> buildSuffixGraph(FastaSequence sequence,  EdgeRuler er, int length) {
    
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
            
            int newSinkId = -1;
            
            // there is a possibility that there might be a match two nodes down
            for (int i=0; i<currentNode.getEdgeCount(); i++) {
              int firstDistance = currentNode.getKeyAt(i);
              if (er.compareEdges(thisDistance, firstDistance) < 0)    break;
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
   * Method that builds the suffix trie with a given maximum depth mass.
   * @param sequence the protein fasta sequence.
   * @return the list of composition nodes built from the fasta sequence. 
   */
  static ArrayList<Node> buildSuffixGraph(FastaSequence sequence, EdgeRuler ec, int length, GraphvizGraph g) {
    
    ArrayList<Node> nodes = new ArrayList<Node>();
    nodes.add(ec.nodeFactory());  // add the root node
    g.setNormalNode(0);
    
    // for every start position
    for(int start = 0; start<length; start++) {
      int currentNodeIndex = 0;
      int prevDistance = 0;
      float cumMass = 0f;
      
      boolean advanced = false;
      
      Node currentNode = nodes.get(currentNodeIndex);
      Node previousNode = null;
      
      // for every amino acid in the database try to extend until the MAX_DEPTH_MASS is reached.
      for(int index=start, depth=0; index<length; depth++, index++) {
        
        AminoAcid aa = AminoAcid.getStandardAminoAcid(sequence.getCharAt(index));
        if(aa == null) {
          // nothing to do, because most likely a separator.
          break;
        }
        else {          
          advanced = true;
     
          // try to add a new edge to the graph
          int thisDistance = ec.fromAa(aa);
          int mIndex = currentNode.getEdge(thisDistance);
          if(mIndex<0) {
            
            int newSinkId = -1;
            
            // there is a possibility that there might be a match two nodes down
            for (int i=0; i<currentNode.getEdgeCount(); i++) {
              int firstDistance = currentNode.getKeyAt(i);
              if (ec.compareEdges(thisDistance, firstDistance) < 0)    break;
              Node secondNode = nodes.get(currentNode.getNodeIdAt(i)); 
              for (int j=0; j<secondNode.getEdgeCount(); j++) {
                int compResult = ec.compareEdges(thisDistance, firstDistance+secondNode.getKeyAt(j));
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
              // add thisEdge to the currentNode, new node, new edge
              Node newSink = ec.nodeFactory();
              currentNode.addEdge(thisDistance, nodes.size(), -mIndex-1);
              previousNode = currentNode; currentNode = newSink;
              
              g.setNormalNode(nodes.size());
              g.setNormalEdge(currentNodeIndex, nodes.size());
              g.setEdgeLabel(currentNodeIndex, nodes.size(), aa.toString());
              currentNodeIndex = nodes.size();

              nodes.add(newSink);
            }
            else {
              // new edge old node
              currentNode.addEdge(thisDistance, newSinkId, -mIndex-1);
              previousNode = currentNode; currentNode = nodes.get(newSinkId);
              
              g.setNormalEdge(currentNodeIndex, newSinkId);
              g.setEdgeLabel(currentNodeIndex, newSinkId, aa.toString());
              currentNodeIndex = newSinkId;
            }
          }
          else {
            // increment the edge count here
            g.setNormalEdge(currentNodeIndex, currentNode.getNodeIdAt(mIndex));
            g.setEdgeLabel(currentNodeIndex, currentNode.getNodeIdAt(mIndex), aa.toString());
            currentNodeIndex = currentNode.getNodeIdAt(mIndex);
            
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
      
      // print the intermediate
      try {
        PrintWriter pw = new PrintWriter("/tmp/graph"+start+".dot");
        pw.println(g.toString());
        pw.close();
      } catch (FileNotFoundException e) {
        e.printStackTrace();
      }
    }
    
    return nodes;
  }
  
  
  /**
   * This routine modifies the graph so that all the gaps are encoded while
   * conserving the DFA properties using a top down recursion.
   * The target must have unique outgoing edge labels.
   * @param nodes the array of nodes.
   * @param seenNodes the bit array marking which nodes have been processed.
   * @param er the object that translates edges to distances.
   */
  static void makeGapLinks(ArrayList<Node> nodes, BitArray seenNodes, EdgeRuler er, PrintWriter statFile) {
    final float MASS_MULT = 10000.0f; 
 
    // upper half represents the mass, lower half represents the node index.
    LongPriorityQueue q = new LongPriorityQueue();
    q.add(0L);
    
    // the number on nodes in this bin
    int[] nodeCount = new int[(int)Constants.MAX_DEPTH_MASS+200];
    // the cumulative degree of all nodes here
    int[] cumDegree = new int[(int)Constants.MAX_DEPTH_MASS+200];
    // the number of nodes that are merged
    int[] mergedNodes = new int[(int)Constants.MAX_DEPTH_MASS+200];
    // the number of nodes that are merged
    int[] newNodes = new int[(int)Constants.MAX_DEPTH_MASS+200];
    // collisions at this mass
    int[] collisions = new int[(int)Constants.MAX_DEPTH_MASS+200];
    
    HashMap<Integer,HashSet<Integer>> children = new HashMap<Integer,HashSet<Integer>>();
    while (!q.isEmpty()) {
      long current = q.poll();
      Node currentNode = nodes.get((int)current);
      
      float thisMass = (current>>>32)/MASS_MULT;
      nodeCount[(int)thisMass]++;
      
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
            // just edge increase, no node increase
          }
          else {
            // in place merging   
            int[] destNodeArray = new int[1];
            destNodeArray[0] = destNodes.iterator().next();
            destNodes.clear();
            mergeNodes(currentNode.getNodeIdAt(matchIndex), destNodeArray, nodes, er, (int)thisMass, newNodes, mergedNodes);
            collisions[(int)thisMass] += 2;
          }
        }
        else if (destNodes.size()>1) {
          // Merge the matches
          if (matchIndex >= 0) {
            int[] destNodeArray = new int[destNodes.size()];
            int j = 0;
            for (int nodeIndex : destNodes) destNodeArray[j++] = nodeIndex;
            destNodes.clear();
            mergeNodes(currentNode.getNodeIdAt(matchIndex), destNodeArray, nodes, er, (int)thisMass, newNodes, mergedNodes);
            collisions[(int)thisMass] += 1+destNodeArray.length;
          }
          else {
            int[] destNodeArray = new int[destNodes.size()];
            int j = 0;
            for (int nodeIndex : destNodes) destNodeArray[j++] = nodeIndex;
            currentNode.addEdge(comp, mergeNodes(destNodeArray, nodes, er, (int)thisMass, newNodes, mergedNodes), -matchIndex-1);
            collisions[(int)thisMass] += 1+destNodeArray.length;
          }
        }
      }
      
      // current node is transformed, add its children to the queue
      for (int i=0; i < currentNode.getEdgeCount(); i++) {
        int nextNodeId = currentNode.getNodeIdAt(i);
        if (!seenNodes.get(nextNodeId)) {
          seenNodes.set(nextNodeId);
          int transition = (int)(er.toFloat(currentNode.getKeyAt(i))*MASS_MULT);
          long next = (((long)((current>>>32)+transition))<<32) | nextNodeId;
          q.add(next);
        }
      }  
      
      cumDegree[(int)thisMass] += currentNode.getEdgeCount();
      
      if (nextTime < System.currentTimeMillis()) {
        nextTime = System.currentTimeMillis()+5000;
        long usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
        System.out.print("----- Used " + usedMem + "MB at " + seenNodes.getSetItems()/1000 + "K nodes at mass " + ((float)(current>>>32))/MASS_MULT + ". ");
        System.out.printf("~%.2f%% done. Queue Size %dK.\n", 100.0*seenNodes.getSetItems()/nodes.size(), q.size()/1000);
      }
    }
    
    int lumpFactor = 20;
    for (int i=0; i+lumpFactor<nodeCount.length; i++) {
      int totalNodes = 0;
      int totalDegrees = 0;
      int totalCollisions = 0;
      for (int j=0; j<lumpFactor; j++, i++) {
        totalNodes += nodeCount[i];
        totalDegrees += cumDegree[i];
        totalCollisions += collisions[i];
      }
      
      if (totalNodes>0) {
        statFile.printf("%d\t%d\t%f\t%f\n", i-lumpFactor, totalNodes, totalDegrees/(float)totalNodes, totalCollisions/(float)totalNodes);
      }
    }
  }
  
  /**
   * In site merge replacing the nId with the resulting merged node.
   * @param nId the id of the node to merge
   * @param targets the other nodes to merge
   * @param nodes the array of all nodes
   * @param er the object that allows edge comparisons.
   */
  private static void mergeNodes(int nId, int[] targets, ArrayList<Node> nodes, EdgeRuler er, int mass, int[] newNodes, int[] mergedNodes) {
        
    Node n = nodes.get(nId);
    mergedNodes[mass]++;
    if (nId == 1) {
      System.out.println(targets.length);
    }
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
            else mergeNodes(n.getNodeIdAt(nMatchIndex), mergingTargetArray, nodes, er, mass, newNodes, mergedNodes);
          }
          else if (mergingTargetArray.length>1){ // we will have to recurse because of conflicting edges
            if (nMatchIndex<0) n.addEdge((int)prevEdge, mergeNodes(mergingTargetArray, nodes, er, mass, newNodes, mergedNodes), -nMatchIndex-1);
            else mergeNodes(n.getNodeIdAt(nMatchIndex), mergingTargetArray, nodes, er, mass, newNodes, mergedNodes);
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
        else mergeNodes(n.getNodeIdAt(nMatchIndex), mergingTargetArray, nodes, er, mass, newNodes, mergedNodes);
      }
      else if (mergingTargetArray.length>1){     // we will have to recurse because of conflicting edges
        if (nMatchIndex<0) n.addEdge((int)prevEdge, mergeNodes(mergingTargetArray, nodes, er, mass, newNodes, mergedNodes), -nMatchIndex-1);
        else mergeNodes(n.getNodeIdAt(nMatchIndex), mergingTargetArray, nodes, er, mass, newNodes, mergedNodes);
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
  private static int mergeNodes(int[] targets, ArrayList<Node> nodes, EdgeRuler er, int mass, int[] newNodes, int[] mergedNodes) {
    
    // create a new node
    int mergedNodeId = nodes.size();
    nodes.add(null);
    
    newNodes[mass]++;
    
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
          else outEdges.add((((long)mergeNodes(mergingTargetArray, nodes, er, mass, newNodes, mergedNodes))<<32)|((int)prevEdge));
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
      else outEdges.add((((long)mergeNodes(mergingTargetArray, nodes, er, mass, newNodes, mergedNodes))<<32)|((int)prevEdge));
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
    
    return mergedNodeId;
  } 
}
