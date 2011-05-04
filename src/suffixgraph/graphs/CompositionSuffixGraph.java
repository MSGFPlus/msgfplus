package suffixgraph.graphs;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;

import msutil.Composition;

import sequences.FastaSequence;
import suffixgraph.Constants;
import suffixgraph.misc.GraphvizGraph;
import suffixgraph.nodes.FinalNode;
import suffixgraph.nodes.Node;
import suffixgraph.nodes.AbstractNode;


/**
 * This class defines a space efficient implementation a SuffixTrie. The trade 
 * off is that queries are limited to a certain mass length defined by a 
 * constant.
 * @author jung
 *
 */
public class CompositionSuffixGraph extends AbstractSuffixGraph<Node> {
  
  
/***** DEFINITIONS *****/  
  /**
   * Constructor taking a graph file path. The user is responsible for providing
   * the same fasta sequence used to construct the graph.
   * @param sequence the fasta sequence object used to construct the graph.
   * @param path the file path and name of the file. It is assumed that the
   *             file name will have .nodes and .edges extension. Do not provide these.
   */
  public CompositionSuffixGraph(FastaSequence sequence, String path) {
    this.sequence = sequence;
    this.er = new Node().getEdgeRuler();
    
    if (!(new File(path)).exists()) {
      ArrayList<Node> nodes = buildSuffixGraph(sequence);
      toFile(path, nodes);
    }
    fromFile(path);
  }
    
  
  /**
   * Constructor taking a limited length for the sequence. This constructor does
   * not write a permanent file into the disk. 
   * @param sequence the fasta sequence object used to construct the graph.
   * @param path restrict the graph to the first length letters
   * @param printStats flag to indicate whether to print statistics
   */
  public CompositionSuffixGraph(FastaSequence sequence, int length, boolean printStats, GraphvizGraph g) {
    this.sequence = sequence;
    this.er = new Node().getEdgeRuler();
    
    String path = Constants.TEMP_GRAPH;
    ArrayList<Node> nodes = null;
    if (printStats) {
      if (g==null) nodes = MeasuredGraphOperations.buildSuffixGraph(sequence, new Node().getEdgeRuler(), length);
      else nodes = MeasuredGraphOperations.buildSuffixGraph(sequence, new Node().getEdgeRuler(), length, g);
    }
    else {
      nodes = buildSuffixGraph(sequence, length);
    }
    toFile(path, nodes);
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
