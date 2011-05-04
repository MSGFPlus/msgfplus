package suffixgraph.test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

import sequences.FastaSequence;
import suffixgraph.Constants;
import suffixgraph.graphs.CompositionGappedSuffixGraph;
import suffixgraph.graphs.CompositionSuffixGraph;
import suffixgraph.graphs.IntegerGappedSuffixGraph;
import suffixgraph.graphs.IntegerSuffixGraph;
import suffixgraph.misc.GraphvizGraph;

public class StatsCollection {

  
  /**
   * Compares the statistics of a composition and integer gapped graph.
   */
  public static void collectCompStats() {
    String userHome = System.getProperty("user.home");
    String fastaFile = null;
    fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    FastaSequence sequence = new FastaSequence(fastaFile);
    final int length = 200000;
    
    CompositionGappedSuffixGraph gsg = new CompositionGappedSuffixGraph(sequence, length, true);
    System.out.println(gsg);
    gsg = null;
    
    IntegerGappedSuffixGraph igsg = new IntegerGappedSuffixGraph(sequence, length, true);
    System.out.println(igsg);
    
    CompositionSuffixGraph csg = new CompositionSuffixGraph(sequence, length, true, null);
    System.out.println(csg.summarizeNodeCounts());
    
  }
 
  
  /**
   * Collect important information as the memory consumption, node count, degree
   * average.
   */
  public static void compareGraphStats() {
    String userHome = System.getProperty("user.home");
    String fastaFile = null;
    fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    FastaSequence sequence = new FastaSequence(fastaFile);
    final int length = 200000;
    
    /*
    CompositionGappedSuffixGraph gsg = new CompositionGappedSuffixGraph(sequence, length, true);
    String outFile = userHome+"/Desktop/compGappedGraph" + length + ".txt";
    try {
      PrintWriter pw = new PrintWriter(outFile);
      pw.print(gsg.summarizeNodeCounts());
      pw.close();
    } 
    catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    */
    
    CompositionSuffixGraph sg = new CompositionSuffixGraph(sequence, length, true, null);
    String outFile = userHome+"/Desktop/compGraph" + length + ".txt";
    try {
      PrintWriter pw = new PrintWriter(outFile);
      pw.print(sg.summarizeNodeCounts());
      pw.close();
      System.out.println(sg);
    } 
    catch (FileNotFoundException e) {
      e.printStackTrace();
    }
  }
  
  
  /**
   * Get the data of graph sizes as a function of the sequence size.
   */
  public static void compareGraphSizes() {
    String userHome = System.getProperty("user.home");
    String fastaFile = null;
    fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    FastaSequence sequence = new FastaSequence(fastaFile);
    String outFile = userHome+"/Desktop/graphSizes.txt";
    
    PrintWriter pw = null;
    try {
      pw = new PrintWriter(outFile);
    } 
    catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    
    for (int length = 5; length < 101; length += 5) {
      CompositionGappedSuffixGraph gsg = new CompositionGappedSuffixGraph(sequence, length*10000, false);
      File f = new File(Constants.TEMP_GRAPH);
      pw.println(length*10000+"\t"+(f.length()/1048576.0)+"\t"+gsg.getAverageOutDegree()+"\tcompGap");
      gsg = null;
      
      CompositionSuffixGraph sg = new CompositionSuffixGraph(sequence, length*10000, false, null);
      f = new File(Constants.TEMP_GRAPH);
      pw.println(length*10000+"\t"+(f.length()/1048576.0)+"\t"+sg.getAverageOutDegree()+"\tcomp");
      sg = null;
      
      /*
      IntegerGappedSuffixGraph igsg = new IntegerGappedSuffixGraph(sequence, length*10000, false);
      f = new File(Constants.TEMP_GRAPH);
      pw.println(length*10000+"\t"+(f.length()/1048576.0)+"\t"+igsg.getAverageOutDegree()+"\tintGap");
      igsg = null;*/
    }
    pw.close();
  }
  
  public static void makeGraphVizCompositions() {
    String userHome = System.getProperty("user.home");
    String fastaFile = null;
    fastaFile = userHome+"/Data/Databases/small.fasta";
    String outFile = userHome+"/Desktop/compGraph.txt";
    
    FastaSequence sequence = new FastaSequence(fastaFile);
    
    int length = 100;
    GraphvizGraph g = new GraphvizGraph();
    CompositionSuffixGraph csg = new CompositionSuffixGraph(sequence, length, true, g);
    CompositionGappedSuffixGraph cgsg = new CompositionGappedSuffixGraph(sequence, length, false);
    try {
      PrintWriter pw = new PrintWriter(outFile);
      pw.print(g.toString(cgsg, csg));
      pw.close();
    } 
    catch (FileNotFoundException e) {
      e.printStackTrace();
    }
  }
  
  public static void makeGraphVizIntegers() {
    String userHome = System.getProperty("user.home");
    String fastaFile = null;
    fastaFile = userHome+"/Data/Databases/small.fasta";
    String outFile = userHome+"/Desktop/integerGraph.txt";
    
    FastaSequence sequence = new FastaSequence(fastaFile);
    
    int length = 100;
    GraphvizGraph g = new GraphvizGraph();
    IntegerSuffixGraph csg = new IntegerSuffixGraph(sequence, length, true, g);
    IntegerGappedSuffixGraph cgsg = new IntegerGappedSuffixGraph(sequence, length, false);
    try {
      PrintWriter pw = new PrintWriter(outFile);
      pw.print(g.toString(cgsg, csg));
      pw.close();
    } 
    catch (FileNotFoundException e) {
      e.printStackTrace();
    }
  }
  
  public static void generateExample() {
    String userHome = System.getProperty("user.home");
    String fastaFile = null;
    
    fastaFile = userHome+"/Data/Databases/example.fasta";
    FastaSequence sequence = new FastaSequence(fastaFile);
    GraphvizGraph g1 = new GraphvizGraph();
    CompositionSuffixGraph example = new CompositionSuffixGraph(sequence, (int)sequence.getSize(), true, g1);
    System.out.println(g1.toString(example));
    
  }
  
  public static void main(String[] args) {
    //collectCompStats();
    compareGraphStats();
    //compareGraphSizes();
    //makeGraphVizCompositions();
    //makeGraphVizIntegers();
    //generateExample();
  }
}
