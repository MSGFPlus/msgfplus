package suffixgraph.graphs;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import msutil.Composition;
import msutil.GappedPeptide;

import sequences.FastaSequence;
import suffixgraph.Constants;
import suffixgraph.misc.BitArray;
import suffixgraph.nodes.FinalIntegerNode;
import suffixgraph.nodes.IntegerNode;
import suffixgraph.nodes.MergedNodeInfo;
import suffixgraph.nodes.Node;
import suffixgraph.nodes.AbstractNode;


/**
 * This is the general class that allows querying gap peptides in time 
 * proportional to the length of the query (excluding the time taken to
 * build the graph). A GappedSuffixGraph is composed of nodes and vertices and
 * has the structure of a direct acyclic graph. This object is satisfies
 * the requirements of a deterministic finite automaton (DFA), in the sense
 * that for each node, all out-going edges have unique labels and no epsilon
 * transitions are allowed.
 * @author jung
 */
public class IntegerGappedSuffixGraph extends AbstractSuffixGraph<IntegerNode> {
 
  
  
/***** DEFINITIONS *****/
  /**
   * Default constructor taking a fasta sequence object and the binary file.
   * @param sequence the protein fasta sequence to construct this trie from.
   * @param path the file path and name of the file. If it doesn't exist create it.
   */
  public IntegerGappedSuffixGraph(FastaSequence sequence, String path) {
    this.sequence = sequence;
    this.er = new IntegerNode().getEdgeRuler();
    
    if (!(new File(path)).exists() || true) {
      ArrayList<Node> nodes = buildSuffixGraph(sequence);
      long usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
      System.out.println("---- Memory used for the simple integer graph " + usedMem + "MB at " + nodes.size()/1000 + "K nodes.");
  
      makeGapLinks(nodes, new ArrayList<ArrayList<MergedNodeInfo>>(nodes.size()), new BitArray());
      
      usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
      System.out.println("---- Memory used for the gapped integer graph " + usedMem + "MB at " + nodes.size()/1000 + "K nodes.");
      toFile(path, nodes);
      nodes.clear();
    }
    fromFile(path);
    System.out.println(toString());
  }
  
  
  /**
   * This constructor does not record the graph file into permanent storage and
   * collect stats of the graph.
   * @param sequence the protein fasta sequence to construct this trie from.
   * @param length the restricted size of this graph.
   * @param printStats determines whether to print extra statistics.
   */
  public IntegerGappedSuffixGraph(FastaSequence sequence, int length, boolean printStats) {
    this.sequence = sequence;
    this.er = new IntegerNode().getEdgeRuler();
    
    String path = Constants.TEMP_GRAPH;
    ArrayList<Node> nodes = null;
    if (printStats) nodes = MeasuredGraphOperations.buildSuffixGraph(sequence, new IntegerNode().getEdgeRuler(), length);
    else nodes = buildSuffixGraph(sequence, length);
    long usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
    System.out.println("---- Memory used for the simple integer graph " + usedMem + "MB at " + nodes.size()/1000 + "K nodes.");
    try {
      String statFile = System.getProperty("user.home")+"/Desktop/intStats.txt";
      PrintWriter pw = new PrintWriter(statFile);
      if (!printStats)  MeasuredGraphOperations.makeGapLinks(nodes, new BitArray(), this.er, pw);
      else              makeGapLinks(nodes, new ArrayList<ArrayList<MergedNodeInfo>>(nodes.size()), new BitArray());
      pw.close();
    } 
    catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    toFile(path, nodes);
    nodes.clear();
    fromFile(path);
  }
  
  
  /**
   * Search the this graph given a GappedPeptide object.
   * @param query the GappedPeptide object.
   * @return true or false for match or mismatch.
   */
  public ArrayList<String> search(GappedPeptide query) {
    return search(query.getCompositions(), query.size());  
  }
  
  
  /**
   * Search the this graph given a sequence of compositions.
   * @param query the sequence of amino acids.
   * @return true or false for match or mismatch.
   */
  public ArrayList<String> search(ArrayList<Composition> query) {
    return search(query, query.size());
  }
  
  
  @Override
  public  ArrayList<String> search(ArrayList<Composition> query, int depth) {
    int[] queryArray = new int[query.size()];
    for (int i=0; i<queryArray.length; i++)
      queryArray[i]=query.get(i).getNominalMass();
    return search(queryArray, depth);
  }
  
  
  @Override
  public void fromFile(String path) {
    this.nodes = new ArrayList<AbstractNode>();
    try {
      DataInputStream nodeFile = new DataInputStream(new BufferedInputStream(new FileInputStream(path)));
      try {        // read all the nodes
        while(true) {
          this.nodes.add(FinalIntegerNode.finalIntegerNodeFactory(nodeFile));
        }
      }
      catch(EOFException _) {
        // do nothing
        nodeFile.close();
        //System.out.println("Read " + this.nodes.size() + " nodes.");
      }
    }
    catch(IOException e) {
      System.err.println(e);
    }
  }
   
}
