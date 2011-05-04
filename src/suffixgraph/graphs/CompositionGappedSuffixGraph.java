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

import sequences.FastaSequence;

import suffixgraph.Constants;
import suffixgraph.misc.BitArray;
import suffixgraph.nodes.AbstractNode;
import suffixgraph.nodes.FinalNode;
import suffixgraph.nodes.MergedNodeInfo;
import suffixgraph.nodes.Node;


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
public class CompositionGappedSuffixGraph extends AbstractSuffixGraph<Node> {
 
  
/***** DEFINITIONS *****/
  /**
   * Default constructor taking a fasta sequence object and the binary file.
   * @param sequence the protein fasta sequence to construct this trie from.
   * @param path the file path and name of the file. If it doesn't exist create it.
   */
  public CompositionGappedSuffixGraph(FastaSequence sequence, String path) {
    this.sequence = sequence;
    this.er = new Node().getEdgeRuler();
    
    if (!(new File(path)).exists()) {
      ArrayList<Node> nodes = buildSuffixGraph(sequence);
      ArrayList<ArrayList<MergedNodeInfo>> mergedNodes = new ArrayList<ArrayList<MergedNodeInfo>>(nodes.size());
      for (int i=0; i<nodes.size(); i++) mergedNodes.add(null);
      makeGapLinks(nodes, mergedNodes, new BitArray());
      toFile(path, nodes);
      nodes.clear();
    }
    fromFile(path);
  }
   
  
  /**
   * This constructor does not record the graph file into permanent storage and
   * collect stats of the graph.
   * @param sequence the protein fasta sequence to construct this trie from.
   * @param length build the graph up to this length.
   * @param printStats determines whether to print extra statistics.
   */
  public CompositionGappedSuffixGraph(FastaSequence sequence, int length, boolean printStats) {
    this.sequence = sequence;
    this.er = new Node().getEdgeRuler();
    
    String path = Constants.TEMP_GRAPH;
    ArrayList<Node> nodes = null;
    
    if (printStats) {
      nodes = MeasuredGraphOperations.buildSuffixGraph(sequence, new Node().getEdgeRuler(), length);
    }
    else {
      nodes = buildSuffixGraph(sequence, length);
    }
    long usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
    System.out.println("---- Memory used for the simple composition graph " + usedMem + "MB at " + nodes.size()/1000 + "K nodes.");
    try {
      String statFile = System.getProperty("user.home")+"/Desktop/compStats.txt";
      PrintWriter pw = new PrintWriter(statFile);
      if (printStats) MeasuredGraphOperations.makeGapLinks(nodes, new BitArray(), this.er, pw);
      else {
        makeGapLinks(nodes, new ArrayList<ArrayList<MergedNodeInfo>>(nodes.size()), new BitArray());
      }
      pw.close();
    } 
    catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    toFile(path, nodes);
    nodes.clear();
    fromFile(path);
  }
  
  
  @Override
  public void fromFile(String path) {
    this.nodes = new ArrayList<AbstractNode>();
    try {
      DataInputStream nodeFile = new DataInputStream(new BufferedInputStream(new FileInputStream(path)));
      try {        // read all the nodes
        while(true) {
          this.nodes.add(FinalNode.finalNodeFactory(nodeFile));
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
  
  @Override
  public ArrayList<String> search(ArrayList<Composition> query, int depth) {
    int[] queryArray = new int[query.size()];
    for (int i=0; i<queryArray.length; i++)
      queryArray[i]=query.get(i).getNumber();
    return search(queryArray, depth);
  }
}
