package suffixtree.trees;

import java.util.ArrayList;
import java.util.HashSet;

import sequences.MassSequence;
import suffixtree.edges.Edge;
import suffixtree.edges.MassEdge;
import suffixtree.matches.ExactMatchObject;
import suffixtree.misc.ProgressMeter;
import suffixtree.nodes.FinalInternalNode;
import suffixtree.nodes.InternalNode;
import suffixtree.nodes.Node;


/**
 * Data structure for matching gapped peptides against a database. This approach
 * bundles the queries together and the database is fed one time, once the 
 * KeywordTree is built.
 * @author jung
 *
 */
public class KeywordTree {
  
  
  private Node root;
  private MassEdge[] queries;
  
  /**
   * Helper that finalizes part of the query tree to speed-up look ups
   * @param n the starting node to run this process
   * @param level the level of the node
   * @return the converted node
   */
  private static Node finalizeNode(Node n, int level) {
    
    for (int i=0; i<n.getDegree(); i++) {
      Edge e = n.getEdgeAt(i);
      e.setSink(finalizeNode(e.getSink(), level+1));
    }
    if (n.getDegree() > 50 || level < 4) {
      return new FinalInternalNode(n);
    }
    return n;
  }
  
  
  /**
   * Get the root of this tree
   * @return the root of this tree
   */
  public Node getRoot() { return this.root; }
  
  
  /**
   * Get the query with the given index
   * @param index the index of the query to retrieve
   * @return the query object
   */
  public ArrayList<Integer> getQueryAt(int index) {
    return queries[index].getLabels();
  }
  
  
  /**
   * Constructor taking a list of list of integer masses.
   * @param queries the array of array of masses to initialize the Keyword tree.
   */
  public KeywordTree(ArrayList<ArrayList<Integer>> queries) {
    this(queries, true);
  }

  
  /**
   * Constructor taking a list of list of integer masses.
   * @param queries the array of array of masses to initialize the Keyword tree.
   * @param optimize if set, the tree will hash the first edges for speed up the search times
   */
  public KeywordTree(ArrayList<ArrayList<Integer>> queries, boolean optimize) {
    // initialize
    this.root = new InternalNode();
    this.queries = new MassEdge[queries.size()];
    
    System.out.printf("Building keyword tree of %d queries:", this.queries.length);
    int count = 0, step = 0;
    for (ArrayList<Integer> iArray : queries) {
      //System.out.println("Inserting " + edge);
      MassEdge edge = new MassEdge(iArray, new InternalNode(count));
      this.queries[count] = edge.duplicate();
      this.insert(edge);
      
      // display progress
      if (count++ >= (step*queries.size())/20.0) {
        step++;
        System.out.printf(" %d%%", step*5);
      }
    }
    System.out.println();
    
    // optimize nodes
    if (optimize) this.root = finalizeNode(this.root, 0);
  }
  
  
  /**
   * Helper class that collects the statistics of the tree.
   * @author jung
   *
   */
  public class TreeStats {
    private int iNodes;
    private int edges;
    private int lNodes;
    
    @Override
    public String toString() {
      StringBuffer sb = new StringBuffer();
      sb.append("Internal nodes " + iNodes + "\n");
      sb.append("Leaf nodes " + lNodes + "\n");
      sb.append("Edges " + edges + "\n");
      return sb.toString();
    }
  }
  
  
  /**
   * Collect the statistics of this data structure. Node and edges.
   * @return the statistics.
   */
  public TreeStats collectStats() {
    TreeStats ts = new TreeStats();
    ts.iNodes++;
    collectStats(this.root, new HashSet<Node>(), ts);
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
      if (n.getDegree()>1) 
        stats.iNodes++;
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
      stats.lNodes++;
    }
  
    for (int index = 0; index < n.getDegree(); index++) {
      
      Node child = n.getEdgeAt(index).getSink();
      stats.edges++;
      
      if (!seenNodes.contains(child)) {
        seenNodes.add(child);
        collectStats(child, seenNodes, stats);
      }
    }
  }

 
  
  /**
   * Register all partial matches into the keyword tree
   */
  /*
  public void mutationMatch(MassSequence db, ArrayList<MutMatchObject> results) {
  
    int rootMaxEdge = root.getMaximumEdge().getLabel();
    //int rootMinEdge = this.root.getMinimumEdge().getLabel();
    
    System.out.print("Searching the tree with 1 mutation: ");
    long compCount = 0, searchCount = 0, nextStop = 0;
    for (long start=0; start<db.getSize(); start++) {
      
      if (!db.hasMass(start)) {
        // shift to the next protein or segment
        continue;
      }
      
      // To avoid translating the db into masses again and again, translate
      // everything into an array and use a shift index to indicate where to start
          
      
      int cumMass = 0;
      // try to match the first edge here
      for (long subStart=start; db.hasMass(subStart); subStart++) {
        //System.out.println("Looking up " + this.db.getSubsequence(start, subStart+1));
        cumMass += db.getIntegerMass(subStart);
          
        // we have reached an impossible path
        if (cumMass > rootMaxEdge) break;

        int matchIndex = root.search(cumMass);
        compCount++;
        if (matchIndex>=0) {
          //System.out.println("Found a match with mass " + cumMass);
          compCount += this.mutationMatch(db, start, subStart+1, root.getEdgeAt(matchIndex), 1, 1, null, results);
        }
          
      }
      
      // for mutation search we have evaluate all possible outgoing edges of this node
      for (int edgeIndex=0; edgeIndex<root.getDegree(); edgeIndex++) {
        Edge targetEdge = root.getEdgeAt(edgeIndex);
        if (targetEdge==null) continue;
        ArrayList<Mutation> muts = new ArrayList<Mutation>();
        Matching.matchDbWithMutation(db, start, targetEdge.getLabel(), muts);
        for (Mutation mutation : muts) {
          
          // this is special case and some sanity checks need to be done
          if (mutation.isInsertion()) {
            // make sure that that insertion is not equal to the previous amino acid
            if (db.hasMass(start-1)) {
              if (db.getIntegerMass(start-1)==Constants.AA.getAminoAcid(mutation.getMutation()).getNominalMass()) continue;
              //if (db.getIntegerMass(start-1)==AminoAcid.getStandardAminoAcid(mutation.getMutation()).getNominalMass()) continue;
            }
          }
          
          compCount += mutationMatch(db, start, mutation.getNextStart(), targetEdge, 1, 1, mutation, results);
        }
      }
      
      searchCount++;
      
      int percDone = (int)(start*100/db.getSize());
      if (percDone>=nextStop) {
        nextStop += 5;
        System.out.printf(" %d%%", percDone);
      }
     
    }
    System.out.println("\nAverage number of comparisons per position " + compCount / searchCount);
  }
  */

  
  
  /**
   * This is the general matching method allowing for 1 mutation.
   * @param db the amino acid sequence 
   * @param dbStart the start position of the database for the current match
   * @param dbCurrent the position to match in the database at this iteration
   * @param edge the edge in the tree to match
   * @param edgeSubIndex the edge subIndex to match because an edge object
   *                     might represent multiple mass edges 
   * @param matchedEdgeCount the number of matched edges so far
   * @param m the mutation object if it has been collected
   * @param results the results
   * @return the number of calls to the search function for the given nodes
   */
  /*
  private int mutationMatch(MassSequence db,
                            long dbStart, 
                            long dbCurrent, 
                            Edge edge, 
                            int edgeSubIndex, 
                            int matchedEdgeCount, 
                            Mutation m, 
                            ArrayList<MutMatchObject> results) {
    
    
    // store the matches
    if (edgeSubIndex==edge.size() && edge.getSink().getPositionsCount()>0 && m!=null) {
      // we have found matches that have a mutation
      for (int queryIndex : edge.getSink().getPositions()) {
        MutMatchObject smmo = new MutMatchObject(dbStart, 
                                                                       dbCurrent,
                                                                       m,
                                                                       db, 
                                                                       this.getQueryAt(queryIndex), 
                                                                       queryIndex);
        results.add(smmo);
      }
    }
    
    // base case
    if (edgeSubIndex==edge.size() && edge.getSink().getDegree()==0) return 0;
    
    int compCount = 0;  // counts the number of times we have to call the search function
    int cumMass = 0;
    if (edgeSubIndex==edge.size()) {
      
      // we have consumed a complete edge from the compressed tree, and we 
      // arrived to a node, so we need to find a matching edge from the node
      Node sink = edge.getSink();
      
      // optimization stuff
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();
    
      // branch many possible gaps
      for (long i=dbCurrent; db.hasMass(i); i++) {
        cumMass += db.getIntegerMass(i);

        // optimization
        if (cumMass < lower) continue;
        if (cumMass > upper) break;
        
        int matchIndex = sink.search(cumMass);
        compCount++;
        if (matchIndex>=0) {
          // recurse
          //System.out.printf("Found a match with mass %d at index %d\n", cumMass, matchedEdgeCount);
          compCount += mutationMatch(db, dbStart, i+1, sink.getEdgeAt(matchIndex), 1, matchedEdgeCount+1, m, results);
        }
      }
      
      // for mutation search we have evaluate all possible outgoing edges of this node
      if (m==null) {
        for (int edgeIndex=0; edgeIndex<sink.getDegree(); edgeIndex++) {
          Edge targetEdge = sink.getEdgeAt(edgeIndex);
          if (targetEdge==null) continue;
          ArrayList<Mutation> muts = new ArrayList<Mutation>();
          Matching.matchDbWithMutation(db, dbCurrent, targetEdge.getLabel(), muts);
          for (Mutation mutation : muts) {
            compCount += mutationMatch(db, dbStart, mutation.getNextStart(), targetEdge, 1, matchedEdgeCount+1, mutation, results);
          }
        }
      }
      
    }
    else {
      
      // this case allows for the direct matching without calling the search
      // method of the node because of the unique edge
      int currentMass = edge.getLabelAt(edgeSubIndex);
      for (long i=dbCurrent; db.hasMass(i); i++) {
        cumMass += db.getIntegerMass(i);
        
        if (cumMass < currentMass) continue;
        // we have reached an impossible path, by over shooting
        if (cumMass > currentMass) break;
        
        // There is a match
        //System.out.printf("Found a match with mass %d at index %d\n", cumMass, matchedEdgeCount);
        compCount += mutationMatch(db, dbStart, i+1, edge, edgeSubIndex+1, matchedEdgeCount+1, m, results);
      }
      
      // the mutation search entails searching the only edge against the db
      if (m==null) {
        ArrayList<Mutation> muts = new ArrayList<Mutation>();
        Matching.matchDbWithMutation(db, dbCurrent, currentMass, muts);
        for (Mutation mutation : muts) {
          //System.out.println("***** Found mutation " + mutation);
          compCount += mutationMatch(db, dbStart, mutation.getNextStart(), edge, edgeSubIndex+1, matchedEdgeCount+1, mutation, results);
        }
      }
      
    }
    
    return compCount;
  }
  */
  
  
  /**
   * Feed in the database and store the matches in an array list.
   * @param db the amino acid fasta database
   * @param results the set of match objects representing the matching
   */
  public void match(MassSequence db, ArrayList<ExactMatchObject> results) {
    
    // no sequences in this tree
    if (this.root.getDegree()==0) return;
    
    int rootMaxEdge = this.getRoot().getMaximumEdge().getLabel();
    int rootMinEdge = this.getRoot().getMinimumEdge().getLabel();
    
    String msg = String.format("Searching db of size %d:", db.getSize());
    ProgressMeter pm = new ProgressMeter(msg, db.getSize(), System.out);
    
    // From the performance standpoint, it is much faster to have a database
    // buffered and read chunks of the db at a time, and proceed to shift
    // this integer array
    ArrayList<Integer> dbb = new ArrayList<Integer>();     // database buffer
    long edgeCount = 0, position = 0;
    while(position < db.getSize()) {
    
      dbb.clear();
      long absStart = position;
      while (db.hasMass(position)) {
        dbb.add(db.getIntegerMass(position));
        position++;
      }

      int dbbSize = dbb.size(); 
      if (dbbSize > 0) {
        
        Integer[] masses = dbb.toArray(new Integer[0]);
        
        // test the different start positions
        for (int start=0; start<dbbSize; start++) {
          int cumMass = 0;
          // try to match the first edge here
          for (int subStart=start; subStart<dbbSize; subStart++) {
              
            cumMass += masses[subStart];
                
            // we have reached an impossible path
            if (cumMass < rootMinEdge) continue;
            if (cumMass > rootMaxEdge) break;
    
            int matchIndex = this.root.search(cumMass);
            if (matchIndex>=0) {
              edgeCount += this.search(db, masses, start, absStart+start, subStart+1, this.root.getEdgeAt(matchIndex), 1, results);
            }
          }
        }
        
        pm.update(position);
      }
      
      position++;     // skip the index with no mass
    }
    
    System.out.println();
    System.out.println("\nAverage number of comparisons per position " + edgeCount/db.getSize());
  }
  
  
  
  /**
   * Feed in the database and retrieve the matches
   * @param database the amino acid fasta database
   * @param outfile the path of the file to write the results
   */
  /*
  public void exactMatch(ProteinFastaSequence database, String outfile) {
    
    // no sequences in this tree
    if (this.root.getDegree()==0) return;
    
    try {
      BufferedWriter os = new BufferedWriter(new FileWriter(outfile));
       
      this.currentDatabase = database;
      long compCount = 0, searchCount = 0, nextStop = 0;
      ArrayList<Integer> masses = new ArrayList<Integer>();
      ArrayList<ExactMatchObject> results = new  ArrayList<ExactMatchObject>();
      for (long start=0; start<database.getSize();) {
        
        if (database.isTerminator(start)) {
          // shift to the next protein
          masses.clear();
          start++;
          continue;
        }
        
        if (masses.size()==0) {      
          
          // build the protein that will be used for many queries
          masses.add(database.getIntegerMass(start));
          long end;
          for (end=start+1; end<database.getSize(); end++) {
            if (database.isTerminator(end)) break;
            masses.add(database.getIntegerMass(end));
          }  
          
          // try for every possible start of the protein
          for (long subStart=start; subStart<end; subStart++) {
            int shift = (int)(subStart-start);
            this.currentStart = subStart;
            this.currentShift = shift;
            
            int cumMass = 0;
            for (int i=shift; i<masses.size(); i++) {
              cumMass += masses.get(i);
              
              // we have reached an impossible path
              //if (cumMass > Constants.MAX_GAP_MASS) break;
              if (cumMass > this.root.getMaximumEdge().getLabel()) break;
              
              int matchIndex = this.root.search(cumMass);
              compCount++;
              if (matchIndex>=0) {
                compCount += this.search(masses, i+1, this.root.getEdgeAt(matchIndex), 1, results);
              }
            }
           
            searchCount++;
          }
          start = end;
        }
        
        if (start>nextStop) {
          nextStop += 100000;
          System.out.println("--- Processed: " + start + " letters. Matches: " + results.size());
        }
        
        // write out the results to file
        for (ExactMatchObject mo : results) {
          os.write(mo.shortSummary()+"\n");
        }
        results.clear();
      }
      System.out.println("Average number of comparisons per position " + compCount / searchCount);
      os.close();

    }
    catch (IOException e) {
      e.printStackTrace();
      System.exit(-9);
    }
  }
  */
  
  
  /**
   * Helper method for the insertion of a mass edge.
   * @param edge the edge to insert to this object.
   */
  private void insert(MassEdge edge) {
    this.root.insert(edge);
  }
    
  
  /**
   * The helper search method
   * @param db the database
   * @param dbb the buffer used to search
   * @param offset the number of positions shifted from dbb
   * @param dbStart the absolute start of the matching in the db 
   * @param dbbCurrent the current shift of dbb for the next match
   * @param edge the compressed edge object to match
   * @param edgeSubIndex the index of the sub edge in the compressed edge
   * @param results the results array
   * @return the number of navigated edges
   */
  private int search(MassSequence db,
                     Integer[] dbb,
                     int offset,
                     long dbStart,
                     int dbbCurrent,
                     Edge edge,
                     int edgeSubIndex,
                     ArrayList<ExactMatchObject> results) {
    
    // commonly used variables
    Node sink = edge.getSink();
    int degree = sink.getDegree();
    
    // there is a query in this node and it matched completely this edge
    if (sink.getPositionsCount()>0 && edgeSubIndex==edge.size()) {
      for (int match : sink.getPositions()) {
        ExactMatchObject mo = new ExactMatchObject(db,
                                                   dbStart, 
                                                   dbStart+dbbCurrent-offset, 
                                                   this.queries[match].getLabels(), 
                                                   match);
        results.add(mo);
      }
    }
    
    
    // base case, leaf node and nothing else to match
    if (degree==0 && edgeSubIndex==edge.size()) return 0;  
    
    
    int edgeCount = 1;
    int cumMass = 0;
    if (edgeSubIndex==edge.size()) {
      // we have consumed a complete edge from the compressed tree, and we 
      // arrived to a node, so we need to find a matching edge from the node
      // This is the slow matching situation because we have many branches to match
      
      // optimization stuff
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();

      // branch many possible paths
      int lowerIndex = 0;
      for (int i=dbbCurrent; i<dbb.length; i++) {
        cumMass += dbb[i];
        
        // optimizations
        if (cumMass < lower) continue;
        if (cumMass > upper) break;
        
        int matchIndex = sink.search(cumMass, lowerIndex, degree);
        if (matchIndex>=0) {
          // recurse
          edgeCount += search(db, dbb, offset, dbStart, i+1, sink.getEdgeAt(matchIndex), 1, results);
          lowerIndex = matchIndex+1;   // resume the match further up
        }
      }
    }
    else {
      // this case allows for the direct matching without calling the search
      // method of the node because of the unique edge. In other words, we can
      // do a greedy match
      
      int currentMass = edge.getLabelAt(edgeSubIndex);
      for (int i=dbbCurrent; i<dbb.length; i++) {
        cumMass += dbb[i];
        
        if (cumMass < currentMass)     continue;
        
        // we have reached an impossible path
        if (cumMass > currentMass)     break;
        
        // There is a match
        edgeCount += search(db, dbb, offset, dbStart, i+1, edge, edgeSubIndex+1, results);
      }
    }
    
    return edgeCount;
  }
  
  
  
  
  /**
   * Helper method to match an edge against an edge.
   * @param masses the entire query
   * @param queryIndex the next index to try to match
   * @param currentMatchingEdge the matching edge
   * @param currentEdgeIndex the next index of the matching edge. This item could be 
   *                  equal to e.size(), which means, we have matched the edge
   *                  completely, and query.getLabelAt(query) should be matched
   *                  to an outgoing edge of e.getSink()
   * @param matches the resulting matches
   * @return the number of comparisons done for this query.
   */
  /*
  private int search(ArrayList<Integer> masses, int queryIndex, Edge currentMatchingEdge, int currentEdgeIndex, ArrayList<ExactMatchObject> matches) {
    int compCount = 0;
    
    // the sequence has matched this path 
    if (currentMatchingEdge.getSink().getPositions().length>0 && currentEdgeIndex==currentMatchingEdge.size()) {
      for (int match : currentMatchingEdge.getSink().getPositions()) {
        ExactMatchObject mo = new ExactMatchObject(this.currentDatabase, 
                                                   this.currentStart, 
                                                   this.currentStart+queryIndex-this.currentShift, 
                                                   this.queries[match].getLabels(), 
                                                   match);
        matches.add(mo);
      }
    }
    
    // leave node and nothing else to match
    if (currentMatchingEdge.getSink().getDegree()==0 && currentEdgeIndex==currentMatchingEdge.size()) {
      // base case
      return 0;  
    }
    
    // base case, run out of query, nothing is matched
    if (queryIndex >= masses.size()) {
      return 0;
    }
    
    int cumMass = 0;
    if (currentEdgeIndex==currentMatchingEdge.size()) {
      // special case try to match all branches in the next node
      Node sink = currentMatchingEdge.getSink();
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();
      for (int i=queryIndex; i<masses.size(); i++) {
        cumMass += masses.get(i);
        if (cumMass < lower) continue;
        if (cumMass > upper) return compCount;
        
        int matchIndex = sink.search(cumMass);
        compCount++;
        if (matchIndex>=0) {
          // recurse
          //System.out.println(i + " Found " + cumMass);
          compCount += search(masses, i+1, sink.getEdgeAt(matchIndex), 1, matches);
        }
      }
    }
    else {
      int currentMass = currentMatchingEdge.getLabelAt(currentEdgeIndex);
      for (int i=queryIndex; i<masses.size(); i++) {
        cumMass += masses.get(i);
        if (cumMass < currentMass) continue;
        
        if (cumMass > currentMass) {
          // we have reached an impossible path
          return compCount;
        }
        
        //System.out.println(queryIndex + " Querying " + cumMass);
        //System.out.println(i + " Found by extension " + cumMass);
        compCount += search(masses, i+1, currentMatchingEdge, currentEdgeIndex+1, matches);
      }
    }
    return compCount;
  }
  */
  
  
  @Override
  public String toString() {
    StringBuffer sb = new StringBuffer();
    for (int i=0; i<root.getDegree(); i++) {
      sb.append(root.getEdgeAt(i)+"\n");
    }
    return sb.toString();
  }
  
}
