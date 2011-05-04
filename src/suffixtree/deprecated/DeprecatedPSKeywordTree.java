package suffixtree.deprecated;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;

import msutil.AminoAcid;


import sequences.MassSequence;
import suffixtree.edges.DirectedMassEdge;
import suffixtree.edges.Edge;
import suffixtree.matches.MatchObject;
import suffixtree.matches.deprecated.PrefixSuffixMatchObject;
import suffixtree.misc.ProgressMeter;
import suffixtree.nodes.ComplexInternalNode;
import suffixtree.nodes.Node;


/**
 * Prefix Suffix Keyword tree for mutation matching.
 * @author jung
 *
 */
public class DeprecatedPSKeywordTree {
  
  
  private ComplexInternalNode fRoot;        // the root of the forward tree
  private ComplexInternalNode rRoot;        // the root of the reverse tree
  private DirectedMassEdge[] queries;       // the list of queries
  private ComplexInternalNode[] traceBacks; // the bridge from one tree to the other
  private MassSequence db;                  // the database
  private int minPartialMatchCount;         // minimum number of edges before start matching
  private int minPartialMatchMass;          // minimum mass before start matching
  
  private int partialMatches = 0;           // variable to keep track of partial matches
  private ProgressMeter collectPM;
  private int collected;
  
  /**
   * Constructor taking a list of list of integer masses as the series of queries
   * to build the keyword tree for.
   * @param queries the array of array of masses to initialize the PSKeyword tree.
   * @param database the protein database object
   * @param minPartialMatchMass only gap peptides with at least this mass will be evaluated
   */
  public DeprecatedPSKeywordTree(ArrayList<ArrayList<Integer>> queries, MassSequence database, int minPartialMatchMass) {
    
    // initialize
    this.fRoot = new ComplexInternalNode(); // the forward tree does not need to have ComplexInternalNodes
    this.rRoot = new ComplexInternalNode();
    this.queries = new DirectedMassEdge[queries.size()];
    this.traceBacks = new ComplexInternalNode[queries.size()];  
    this.db = database;
    this.minPartialMatchMass = minPartialMatchMass;
    this.minPartialMatchCount = 3;
    
    DirectedMassEdge[] rQueries = new DirectedMassEdge[queries.size()];
    
    int queryIndex = 0;
    String msg = String.format("Building PS-KeywordTrees of %d queries", queries.size());
    ProgressMeter pm = new ProgressMeter(msg, queries.size(), System.out);
    for (ArrayList<Integer> iArray : queries) {
      
      // create the reverse version of the query 
      ArrayList<Integer> reversed = new ArrayList<Integer>(iArray);
      Collections.reverse(reversed);
      
      ComplexInternalNode reverseLeafNode = new ComplexInternalNode(queryIndex);
      //reverseLeafNode.setParentNode(this.rRoot);
      DirectedMassEdge rEdge = new DirectedMassEdge(reversed, reverseLeafNode);
      rQueries[queryIndex] = rEdge.duplicate();
      this.rRoot.insert(rEdge);
      
      //System.out.println("Inserting " + edge);
      ComplexInternalNode forwardLeafNode = new ComplexInternalNode(queryIndex);
      DirectedMassEdge edge = new DirectedMassEdge(iArray, forwardLeafNode);
      this.queries[queryIndex] = edge.duplicate();
      this.fRoot.insert(edge);
      
      // display progress
      pm.update(++queryIndex);
    }
    System.out.println();

    // initialize the tracebacks once the full tree is constructed, reconnect the parents, for the reverse edges
    queryIndex = 0;
    for (DirectedMassEdge edge : rQueries) {
      ComplexInternalNode n = this.rRoot;
      ComplexInternalNode prev = n;
      for (int edgeIndex=0; edgeIndex<edge.size();) {
        Edge matchedEdge = n.getEdgeAt(n.search(edge.getLabelAt(edgeIndex))); 
        prev = n;
        n = (ComplexInternalNode)matchedEdge.getSink();
        edgeIndex += matchedEdge.size();
        n.setParentNode(prev);
      }
      this.traceBacks[queryIndex++] = n;
    }
    
    // populate both trees
    this.populateForward();
    this.populateReverse();
  }
  
  
  
  /**
   * Collect prefix-suffix matches on this tree and database. There minimum mass
   * value for a prefix and suffix anchor is specified in Constants class.
   * @return the list of matched objects
   */
  public ArrayList<MatchObject> collectPrefixSuffixMatches() {
    this.collectPM = new ProgressMeter("Collecting", this.queries.length, System.out);
    this.collected = 0;
    HashMap<PrefixSuffixMatchObject,PrefixSuffixMatchObject> matches = new HashMap<PrefixSuffixMatchObject,PrefixSuffixMatchObject>();
    //collectPrefixSuffixMatches(this.fRoot, new Stack<ComplexInternalNode>(), matches);
    collectPrefixSuffixMatches(this.fRoot, new TreeMap<Integer,Integer>(), 0, Integer.MAX_VALUE, new TreeMap<Integer,ArrayList<Coor>>(), matches);
    return new ArrayList<MatchObject>(matches.values());
  }


  
  //the class to store the coordinates
  private class Coor {
    private long start, end;
    private Coor(long start, long end) { this.start=start; this.end=end; }
  }
  
  /**
   * Helper method to collect prefix-suffix matches.
   * @param node the node to recurse on, in the prefix tree
   * @param path the stack of the path taken so far
   * @param prefixes the coordinates collected so far
   * @param matches store the matches in this object
   */
  private void collectPrefixSuffixMatches(ComplexInternalNode node, 
                                          TreeMap<Integer,Integer> masses,
                                          int massToRoot,
                                          int massCutOff,
                                          TreeMap<Integer,ArrayList<Coor>> prefixes,
                                          HashMap<PrefixSuffixMatchObject,PrefixSuffixMatchObject> matches) {
 
    // these are the masses added by adding the current node, we will have to 
    // remove this at the end of this function
    ArrayList<Integer> addedMasses = new ArrayList<Integer>();
    
    // coordinates added by this node
    ArrayList<Coor> coors = new ArrayList<Coor>();
    for (int i=0; i<node.getPrefixMatchCount(); i++) {
      long start = node.getPrefixStartAtIndex(i), end = start+node.getPrefixExtendAtIndex(i);
      int dbCumMass = 0;
      for (long currentStart=start; currentStart<end; currentStart++) {
        dbCumMass += this.db.getIntegerMass(currentStart);
        
        // we need a minimum mass
        if (dbCumMass < massCutOff) continue;
        
        addedMasses.add(dbCumMass);
        if (!prefixes.containsKey(dbCumMass)) prefixes.put(dbCumMass, new ArrayList<Coor>());
        prefixes.get(dbCumMass).add(new Coor(start, currentStart+1));
      }
      
      // try to extend one more edge, but do not reach the next edge
      if (masses.containsKey(dbCumMass)) {
        int nextEdge = masses.get(dbCumMass);
        int cumMass = 0;
        for (long currentStart=end; currentStart<this.db.getSize(); currentStart++) {
          if (!this.db.hasMass(currentStart)) break;
          
          cumMass += this.db.getIntegerMass(currentStart);
          if (cumMass >= nextEdge) break;
          
          int shiftedMass = dbCumMass+cumMass;
          
          addedMasses.add(shiftedMass);
          if (!prefixes.containsKey(shiftedMass)) prefixes.put(shiftedMass, new ArrayList<Coor>());
          prefixes.get(shiftedMass).add(new Coor(start, currentStart+1));
        }
      }
      else {
        // this is the last match, we extend it later
        coors.add(new Coor(start, end));
      }
      
      //System.out.println("Added prefix " + this.db.getSubsequence(start, end) + " " + dbCumMass);
    }
    
    // when there is a query or are many queries, we construct the match objects
    if (node.getPositions().length>0) {
      
      this.collectPM.update(++this.collected);
      
      // Test all the queries
      for (int queryIndex : node.getPositions()) {
        //System.out.println("Retrieve reverse path from query index " + queryIndex);
        
        DirectedMassEdge query = this.queries[queryIndex];
        int addedMatches = 0;
        int newMatches = 0;
        int parentMass = query.getTotalMass();
        ComplexInternalNode currentNode = this.traceBacks[queryIndex];
        while (true) {
          // collect all the suffixes
          for (int i=0; i<currentNode.getPrefixMatchCount(); i++) {
            long start = currentNode.getPrefixStartAtIndex(i), end = start+currentNode.getPrefixExtendAtIndex(i);
            
            // we collect all the partial matches
            int cumMass = 0;
            int suffixMassEdgeIndex = query.size()-1;
            int matchedEdges = 0;
            int suffixMassEdgeToMatch = query.getLabelAt(suffixMassEdgeIndex);
            boolean matchOneMore = false;
            for (long currentEnd=end-1; currentEnd>=0; currentEnd--) {
              
              if (!this.db.hasMass(currentEnd)) break;
              
              // add the mass to the cumulative masses
              cumMass += this.db.getIntegerMass(currentEnd);
              
              // we have reached (or exceeded) the parent mass
              if (cumMass >= parentMass) break;
              
              // we have matched the last edge to match
              if (matchOneMore && suffixMassEdgeToMatch<=cumMass) break;
              
              // we have matched the current edge, update to the next one
              if (suffixMassEdgeToMatch==cumMass) {
                matchedEdges++;
                suffixMassEdgeToMatch += query.getLabelAt(--suffixMassEdgeIndex);
              }
              
              // signal that we only have to accumulate one more edge (for bounding)
              if (currentEnd==start) matchOneMore = true; 
              
              // have not reached the minimum requirements
              if (cumMass < this.minPartialMatchMass) continue;
              
              // have not reached the minimum matched edges
              if (matchedEdges < this.minPartialMatchCount) continue;
              
              // finally add check whether there is a prefix match
              int targetMass = parentMass - cumMass;
              if (prefixes.containsKey(targetMass)) {
                for (Coor c : prefixes.get(targetMass)) {
                  PrefixSuffixMatchObject m = new PrefixSuffixMatchObject(c.start, 
                                                                          c.end, 
                                                                          currentEnd, 
                                                                          end,
                                                                          this.db,
                                                                          query.getLabels(),
                                                                          queryIndex);
                  long relStart = m.getStart();// m.getRelativeStart();
                  long relEnd = m.getEnd(); //m.getRelativeEnd();
                  if (relEnd - relStart < 800 && relEnd - relStart > 30) {
                    addedMatches++;
                    // equality is defined by the start/end coordinates
                    if (!matches.containsKey(m)) {
                      matches.put(m, m);
                      newMatches++;
                    }
                    else {
                      matches.get(m).addMiddle(c.end, currentEnd);
                    }
                    
                    // check the correctness
                    // turn on the checks
                    if (!m.getPeptide().isCorrect(query.getLabels())) {
                      for (AminoAcid aa : m.getPeptide()) {
                        System.out.println(aa.getNominalMass());
                      }
                      System.out.println(this.db.getCharAt(c.end) + " " + this.db.hasMass(c.end));
                      System.out.println(this.db.getCharAt(currentEnd) + " " + this.db.hasMass(currentEnd));
                      System.out.printf("Peptide %s is not matched to query %s\n", m.getMatchAsString(), query.getLabels().toString());
                      System.exit(-9);
                    }
                  }            
                }
              }
              
            }
          }
          
          // we have reached the root
          if (currentNode.getParentNode()==null) break;
          
          currentNode = currentNode.getParentNode();
        }
      }
    }
    
    // recursively call the function
    for (int i=0; i<node.getDegree(); i++) {
      
      // first thing we need to do is to extend the last node's matches until right
      // before we reach the mass of the edge we are going to recurse on
      Edge e = node.getEdgeAt(i);
      ComplexInternalNode nextNode = (ComplexInternalNode)e.getSink();
      int massToMatch = node.getEdgeAt(i).getLabel(); 
      //System.out.println("Size of the edge " + node.getEdgeAt(i).size());
      ArrayList<Integer> extendMasses = new ArrayList<Integer>();     // will remove later
      for (Coor coor : coors) {
        int cumMass = 0;
        for (long currentStart=coor.end; currentStart<this.db.getSize(); currentStart++) {
          
          // no need to add anything else
          if (!this.db.hasMass(currentStart)) break;
          
          cumMass += this.db.getIntegerMass(currentStart);
          if (cumMass >= massToMatch) break;
          
          // add the extension
          int shiftedMass = cumMass + massToRoot;
          extendMasses.add(shiftedMass);
          if (!prefixes.containsKey(shiftedMass)) prefixes.put(shiftedMass, new ArrayList<Coor>());
          prefixes.get(shiftedMass).add(new Coor(coor.start, currentStart+1));
        }
      }
      
    
      
      // prepare the masses array
      int cumTotalMass = massToRoot;
      int nextMassCutOff = massCutOff;
      for (int edgeIndex=0; edgeIndex<e.size(); edgeIndex++) {
        masses.put(cumTotalMass, e.getLabelAt(edgeIndex));
        cumTotalMass += e.getLabelAt(edgeIndex);
        if (masses.size()==this.minPartialMatchCount) {
          nextMassCutOff=cumTotalMass;
        }
      }

      
      // recurse
      collectPrefixSuffixMatches(nextNode, masses, cumTotalMass, nextMassCutOff, prefixes, matches);

      
      // remove the last mass
      masses.tailMap(massToRoot, true).clear(); // remove everything added
      
      // now remove the added masses for the next recursion
      for (int massToRemove : extendMasses) {
        if (prefixes.get(massToRemove).size()==1) {
          // remove the entry
          prefixes.remove(massToRemove);
        }
        else {
          // remove the last item of the array
          ArrayList<Coor> values = prefixes.get(massToRemove);
          values.remove(values.size()-1);
        }
      }
    }
    
    // remove the added masses to the prefix
    for (int massToRemove : addedMasses) {
      if (prefixes.get(massToRemove).size()==1) {
        // remove the entry
        prefixes.remove(massToRemove);
      }
      else {
        // remove the last item of the array
        ArrayList<Coor> values = prefixes.get(massToRemove);
        values.remove(values.size()-1);
      }
    } // prefix data structure is restored to the original state
    
  }
  
  
  /**
   * Reverse populating the tree
   */
  private void populateReverse() {
    
    this.partialMatches = 0;
    int rootMaxEdge = this.rRoot.getMaximumEdge().getLabel();
    
    String msg = String.format("Initializing suffix partial matches with %d db", this.db.getSize());
    ProgressMeter pm = new ProgressMeter(msg, this.db.getSize(), System.out);
    long compCount = 0, searchCount = 0;
    for (long start=this.db.getSize()-1; start>=0; start--) {
      
      if (!this.db.hasMass(start)) {
        continue;  // skip this position
      }
      
      int cumMass = 0;
      int divider = 0;
      for (long i=start; i>=0 && this.db.hasMass(i); i--) {
        cumMass += this.db.getIntegerMass(i);
            
        // we have reached an impossible path
        if (cumMass > rootMaxEdge) break;

        int matchIndex = this.rRoot.search(cumMass, divider);
        compCount++;
        if (matchIndex>=0) {
          divider = matchIndex + 1;
          //System.out.println("Matched mass " + cumMass);
          compCount += this.storeMatchesReverse(start, i-1, this.rRoot.getEdgeAt(matchIndex), 1, cumMass, 1);
        }
      }
      
      searchCount++;
      pm.update(this.db.getSize()-start-1);
    }
    System.out.printf("\nStored %d partial matches coordinates in the tree\n", this.partialMatches);
  }
  
 
  
  /**
   * Register all partial matches into the keyword tree
   */
  private void populateForward() {
  
    this.partialMatches = 0;
    int rootMaxEdge = this.fRoot.getMaximumEdge().getLabel();
    
    String msg = String.format("Initializing prefix partial matches with %d db", this.db.getSize());
    ProgressMeter pm = new ProgressMeter(msg, this.db.getSize(), System.out);
    long compCount = 0, searchCount = 0;
    for (long start=0; start<this.db.getSize(); start++) {
      
      if (!this.db.hasMass(start)) {
        continue; // skip this position
      }
     
      int cumMass = 0;
      int divider = 0;
      // try to match the first edge here
      for (long i=start; i<this.db.getSize() && this.db.hasMass(i); i++) {
        cumMass += this.db.getIntegerMass(i);
            
        // we have reached an impossible path
        if (cumMass > rootMaxEdge) break;

        int matchIndex = this.fRoot.search(cumMass, divider);
        compCount++;
        if (matchIndex>=0) {
          divider = matchIndex+1;
          compCount += this.storeMatchesForward(start, i+1, this.fRoot.getEdgeAt(matchIndex), 1, cumMass, 1);
        }
      }
      
      searchCount++;
      pm.update(start);
    }
    System.out.printf("\nStored %d partial matches coordinates in the tree\n", this.partialMatches);
    //System.out.println("\nAverage number of comparisons per position " + compCount/searchCount);
  }

  

  /**
   * Helper method to store all partial matches when searching the db backwards
   * @param dbStart the start position of the db where matching started
   * @param dbIndex the index of the db to match
   * @param e the edge to match
   * @param eIndex the index of the current edge where matching should occur. 
   *               This item could be equal to e.size(), which means, we have 
   *               matched the edge completely, and we should move on to the
   *               node and evaluate the edges.
   * @param matchedMass the matched mass so far
   * @param matchedEdges the number of edges matched so far
   * @return the number of comparisons done for this query.
   */
  private int storeMatchesReverse(long dbStart, long dbIndex, Edge e, int eIndex, int matchedMass, int matchedEdges) {
    int compCount = 0;
    
    ComplexInternalNode lastNode = null;    // the node to store the prefix matches
    boolean hasMatch = false;               // flag to see whether there was match
    
    Node sink = e.getSink();
    
    // We have reached the leaf of this branch OR we ran out of database
    if ((eIndex==e.size() && sink.getDegree()==0) || dbIndex < 0 || !this.db.hasMass(dbIndex)) {
      // base case 
      lastNode = (ComplexInternalNode)sink;
    }
    else {
      
      int cumMass = 0;
      if (eIndex==e.size()) {
        // we have consumed a complete edge from the compressed tree, and we 
        // arrived to a node, so we need to find a matching edge from the node
        
        // optimization stuff
        int lower = sink.getMinimumEdge().getLabel();
        int upper = sink.getMaximumEdge().getLabel();
      
        // branch many possible gaps
        int divider = 0;
        for (long i=dbIndex; i>=0 && this.db.hasMass(i); i--) {
          cumMass += this.db.getIntegerMass(i);
  
          // optimization
          if (cumMass < lower) continue;
          if (cumMass > upper) break;
          
          int matchIndex = sink.search(cumMass, divider);
          compCount++;
          if (matchIndex>=0) {
            divider = matchIndex + 1;
            // recurse
            compCount += storeMatchesReverse(dbStart, i-1, sink.getEdgeAt(matchIndex), 1, matchedMass+cumMass, matchedEdges+1);
            hasMatch = true;
          }
        }
        lastNode = (ComplexInternalNode)sink;
      }
      else {
        
        // we can do a greedy match because there is only one edge
        int currentMass = e.getLabelAt(eIndex);
        for (long i=dbIndex; i>=0 && this.db.hasMass(i); i--) {
          cumMass += this.db.getIntegerMass(i);
          
          if (cumMass < currentMass) continue;
          // we have reached an impossible path, by over shooting
          if (cumMass > currentMass) break;
            
          // There is a match
          compCount += storeMatchesReverse(dbStart, i-1, e, eIndex+1, matchedMass+cumMass, matchedEdges+1);
          hasMatch = true;
        }
        lastNode = ((ComplexInternalNode)e.getSink());
      }
    }
    
    // add only those with a the minimum number of matched edges
    if (!hasMatch && matchedMass >= this.minPartialMatchMass && matchedEdges >= this.minPartialMatchCount) {
      //System.out.println("Added suffix " + this.db.getSubsequence(dbIndex+1, dbStart));
      lastNode.addPartialMatch(dbIndex+1, (int)(dbStart-dbIndex));
      //System.out.printf("\n%d\t%d\n", dbIndex+1, dbStart);
      this.partialMatches++;
    }
    return compCount;
  }
  
  
  /**
   * Helper method to store all partial matches searching the db forward
   * @param dbStart the start position of the db where matching started
   * @param dbIndex the index of the db to match
   * @param e the edge to match
   * @param eIndex the index of the current edge where matching should occur. 
   *               This item could be equal to e.size(), which means, we have 
   *               matched the edge completely, and we should move on to the
   *               node and evaluate the edges.
   * @param matchedMass the matched mass so far
   * @param matchedEdges the number of edges matched so far
   * @return the number of comparisons done for this query.
   */
  private int storeMatchesForward(long dbStart, long dbIndex, Edge e, int eIndex, int matchedMass, int matchedEdges) {
    int compCount = 0;
    
    ComplexInternalNode lastNode = null;    // the node to store the prefix matches
    boolean hasMatch = false;               // flag to see whether there was match
    
    Node sink = e.getSink();
    
    // We have reached the leaf of this branch OR we ran out of database
    if ((eIndex==e.size() && sink.getDegree()==0) || dbIndex >= this.db.getSize() || !this.db.hasMass(dbIndex)) {
      // base case 
      lastNode = (ComplexInternalNode)sink;
    }
    else {
      
      int cumMass = 0;
      if (eIndex==e.size()) {
        // we have consumed a complete edge from the compressed tree, and we 
        // arrived to a node, so we need to find a matching edge from the node
        
        // optimization stuff
        int lower = sink.getMinimumEdge().getLabel();
        int upper = sink.getMaximumEdge().getLabel();
      
        // branch many possible gaps
        int divider = 0;
        for (long i=dbIndex; i<db.getSize() && this.db.hasMass(i); i++) {
          cumMass += this.db.getIntegerMass(i);
  
          // optimization
          if (cumMass < lower) continue;
          if (cumMass > upper) break;
          
          int matchIndex = sink.search(cumMass, divider);
          compCount++;
          if (matchIndex>=0) {
            divider = matchIndex + 1;
            // recurse
            compCount += storeMatchesForward(dbStart, i+1, sink.getEdgeAt(matchIndex), 1, matchedMass+cumMass, matchedEdges+1);
            hasMatch = true;
          }
        }
        lastNode = (ComplexInternalNode)sink;
      }
      else {
        
        // we can do a greedy match because there is only one edge
        int currentMass = e.getLabelAt(eIndex);
        for (long i=dbIndex; i<this.db.getSize() && this.db.hasMass(i); i++) {
          cumMass += this.db.getIntegerMass(i);
          
          if (cumMass < currentMass) continue;
          // we have reached an impossible path, by over shooting
          if (cumMass > currentMass) break;
            
          // There is a match
          compCount += storeMatchesForward(dbStart, i+1, e, eIndex+1, matchedMass+cumMass, matchedEdges+1);
          hasMatch = true;
        }
        lastNode = ((ComplexInternalNode)e.getSink());
      }
    }
    
    // add only those with a the minimum number of matched edges
    if (!hasMatch && matchedMass >= this.minPartialMatchMass && matchedEdges >= this.minPartialMatchCount) {
      //System.out.println("Added prefix " + this.db.toString(this.currentStart, queryIndex));
      lastNode.addPartialMatch(dbStart, (int)(dbIndex-dbStart));
      this.partialMatches++;
    }
    return compCount;
  }
  
  
  @Override
  public String toString() {
    StringBuffer sb = new StringBuffer();
    for (int i=0; i<fRoot.getDegree(); i++) {
      sb.append(fRoot.getEdgeAt(i)+"\n");
    }
    return sb.toString();
  }
  

}
