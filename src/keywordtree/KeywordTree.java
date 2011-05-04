package keywordtree;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import msgap.results.GappedPeptideResults;

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
 * KeywordTree is built. This version does not require extra memory for storing
 * matches
 * @author jung
 *
 */
public class KeywordTree {
  
  private static final int MATCH_FILE_COUNT = 100; 
  public static final String MATCH_FILE_PREFIX = "matchFile";
  
  
  // member variables
  private Node root;
  private MassEdge[] queries;
  private GappedPeptideResults gpr;
  private PrintWriter[] matchFiles;
  private HashMap<Integer,Integer> specId2matchFileIndex;
  private int matchedQueries;
  
  
  
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
   * Constructor taking the gapped peptide results object
   * @param gpr contains the queries and other meta information.
   */
  public KeywordTree(GappedPeptideResults gpr) {
    // initialize
    this.root = new InternalNode();
    this.queries = new MassEdge[gpr.getSequences().size()];
    
    ProgressMeter pm = new ProgressMeter("Building keyword tree of "+this.queries.length, this.queries.length, System.out);
    int count = 0;
    for (ArrayList<Integer> iArray : gpr.getSequences()) {
      //System.out.println("Inserting " + edge);
      MassEdge edge = new MassEdge(iArray, new InternalNode(count));
      this.queries[count] = edge.duplicate();
      this.insert(edge);
      
      pm.update(++count);
    }
    System.out.println();
    
    // optimize nodes
    this.root = finalizeNode(this.root, 0);
    this.gpr = gpr;
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
   * Helper method to delete file recursively from a directory
   * @param path the path to delete
   * @return true if the deletion was successful
   */
  static private boolean deleteDirectory(File path) {
    if(path.exists()) {
      File[] files = path.listFiles();
      for(int i=0; i<files.length; i++) {
        if(files[i].isDirectory()) {
          deleteDirectory(files[i]);
        }
        else {
          files[i].delete();
        }
      }
    }
    return(path.delete());
  }
  
  
  /**
   * Close each of the open files
   */
  public void closeMatchFiles() {
    for (PrintWriter pw : this.matchFiles) {
      if (pw!=null) pw.close();
    }
  }
  
  
  /**
   * Open the match files
   * @param matchDirectory
   *
   */
  private void openMatchFiles(String matchDirectory, boolean append) {
    
    try {
      File path = new File(matchDirectory); 
      if (path.exists()) {
        // delete the file or directory
        if (!append) deleteDirectory(path);
      }
      // create the directory
      path.mkdirs();
   
      // create the list of PrintWriters
      this.matchFiles = new PrintWriter[MATCH_FILE_COUNT];
      
      // initialize the specId to file array index
      this.specId2matchFileIndex = new HashMap<Integer,Integer>();
      int fileIndex = 0, specCount = 0;
      float nextStop = gpr.getSpecIds().size()/(float)MATCH_FILE_COUNT;
      for (int specId : gpr.getSpecIds()) {
        this.specId2matchFileIndex.put(specId, fileIndex);

        if (this.matchFiles[fileIndex]==null) {
          // create the print writer object
          this.matchFiles[fileIndex] = new PrintWriter(new FileWriter(new File(matchDirectory, String.format(MATCH_FILE_PREFIX+"%04d.txt", fileIndex)), append));
        }
        
        specCount++;
        if (specCount > nextStop) {
          fileIndex++;
          nextStop += gpr.getSpecIds().size()/(float)MATCH_FILE_COUNT;
        }
      }
    }
    catch(IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
  }
  
  
  
  /**
   * Feed in the database and store the matches in an array list.
   * @param db the amino acid fasta database
   * @param matchDir the directory where to store the matches
   * @param append determines the mode of match writing. Append means do not clear the match directory
   */
  public int match(MassSequence db, String matchDir, boolean append) {
    
    // no sequences in this tree
    if (this.root.getDegree()==0) return 0;
    
    long time = System.currentTimeMillis();
    
    this.matchedQueries = 0;
    
    openMatchFiles(matchDir, append);
    
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
              edgeCount += this.search(db, masses, start, absStart+start, subStart+1, this.root.getEdgeAt(matchIndex), 1);
            }
          }
        }
        
        pm.update(position);
      }
      
      position++;     // skip the index with no mass
    }
    
    System.out.println();
    System.out.println("\nAverage number of comparisons per position " + edgeCount/db.getSize());
    time = System.currentTimeMillis() - time;
    System.out.println("Time elapse for matching " + time/1000 + " seconds");
    
    return this.matchedQueries;
  }
  
  
  
  
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
   * @return the number of navigated edges
   */
  private int search(MassSequence db,
                     Integer[] dbb,
                     int offset,
                     long dbStart,
                     int dbbCurrent,
                     Edge edge,
                     int edgeSubIndex) {
    
    // commonly used variables
    Node sink = edge.getSink();
    int degree = sink.getDegree();
    
    // there is a query in this node and it matched completely this edge
    if (sink.getPositionsCount()>0 && edgeSubIndex==edge.size()) {
      for (int queryIndex : sink.getPositions()) {
        ExactMatchObject mo = new ExactMatchObject(db,
                                                   dbStart, 
                                                   dbStart+dbbCurrent-offset, 
                                                   this.gpr.getQueryAt(queryIndex), 
                                                   queryIndex);
        int matchFileIndex = this.specId2matchFileIndex.get(this.gpr.getSpecId(queryIndex));
        this.matchFiles[matchFileIndex].println(mo.toString());
        this.matchedQueries++;
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
          edgeCount += search(db, dbb, offset, dbStart, i+1, sink.getEdgeAt(matchIndex), 1);
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
        edgeCount += search(db, dbb, offset, dbStart, i+1, edge, edgeSubIndex+1);
      }
    }
    
    return edgeCount;
  }
  
  
  
  @Override
  public String toString() {
    StringBuffer sb = new StringBuffer();
    for (int i=0; i<root.getDegree(); i++) {
      sb.append(root.getEdgeAt(i)+"\n");
    }
    return sb.toString();
  }
  
}
