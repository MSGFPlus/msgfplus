package suffixtree.test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

import sequences.FastaSequence;
import suffixtree.Constants;
import suffixtree.edges.ByteEdge;
import suffixtree.trees.GappedSuffixTree;

public class GappedSuffixTreeTest {
 
  
  public static String gappedQueryString(FastaSequence sequence, ArrayList<ByteEdge> query) {
    StringBuffer sb = new StringBuffer();
    for(ByteEdge qe : query) {
      sb.append("[" + sequence.toString(qe.toByteArray()) + "," + qe.getLabel() + "] ");
    }
    return sb.toString();
  }
  
  
  public static long query(GappedSuffixTree st, ArrayList<ByteEdge> query) {
    long time = System.currentTimeMillis();
    //System.out.println("Querying: " + gappedQueryString(st.getSequence(), query));
    HashSet<Integer> matches = new HashSet<Integer>();
    st.search(query, matches);
    time = System.currentTimeMillis() - time;
    
    if (matches.size()==0) {
      System.out.print("Querying: " + gappedQueryString(st.getSequence(), query));
      System.out.println(" - Not found!");
      return -1;
    }
    else {
      /*
      if (matches.size() > 10 || time > 1000) {
        System.out.println(" in " + time + " milisecs");
        System.out.print(" - [" + matches.size() + "] positions: ");
        for (int i : matches) {
          System.out.println("     " + i + " : " + st.getSequence().toString(i, i+MIN_QUERY_LENGTH));
        }
        System.out.println();
      }*/
    }
    
    return time;
  }
  
  
  public static void queryAll(FastaSequence sequence, GappedSuffixTree st) {
    Random r = new Random();
    
    int queryCount = 0;
    long cumTime = 0;
    // query all the items
    for (int start=0; start < sequence.getSize(); start++) {
      if (sequence.isTerminator(start)) continue;
      
      // for each start position, find the end
      int end = start;
      for (int i=start+1; i < sequence.getSize(); i++) {
        if (sequence.isTerminator(i)) {
          end = i;
          break;
        }
      }
      
      end = Math.min(end, start+Constants.MAX_QUERY_CHAR);
      if (end-start >= Constants.MIN_QUERY_CHAR) {
      
        // intact query
        ArrayList<ByteEdge> query = new ArrayList<ByteEdge>();
        for (int i=start; i<Math.min(end, start+Constants.MIN_QUERY_CHAR); i++) {
          query.add(new ByteEdge(sequence.getByteAt(i)));
        }
        queryCount++;
        cumTime += query(st, query);
        
        //gapped query
        ArrayList<ByteEdge> gappedQuery = new ArrayList<ByteEdge>();
        for (int i=start; i<Math.min(end, start+Constants.MIN_QUERY_CHAR);) {
          int gapSize = r.nextInt(Constants.MAX_GAP)+1;
          if (i+gapSize >= end) break;
          gappedQuery.add(new ByteEdge(sequence.getBytes(i, i+gapSize)));
          i += gapSize;
        }        
        
        if (gappedQuery.size()==0) continue;
        
        queryCount++;
        cumTime += query(st, gappedQuery);
      }
    }
    
    System.out.printf("-- %d queries in %.2f seconds\n", queryCount, cumTime/1000.0);
    System.out.printf("-- Average %.2f ms per 1000 query", 1000.0*cumTime/queryCount);
  }
  
  
  public static void main(String[] args) {
    String fastaFile;
    String userHome = System.getProperty("user.home");
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    //fastaFile = userHome+"/Data/Databases/small.fasta";
    //fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/repeat.fasta";
    //fastaFile = userHome+"/Data/Databases/half.fasta";
    //fastaFile = userHome+"/Data/Databases/million.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    fastaFile = userHome+"/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants.fasta";
    
    //String graphFile = fastaFile.replaceAll(".fasta$", "")+".st";
   
    long time = System.currentTimeMillis();
    FastaSequence sequence = new FastaSequence(fastaFile);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("--- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    GappedSuffixTree st = new GappedSuffixTree(sequence);
    System.out.println("-- Loading SuffixTree file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    //System.out.print(st);

    queryAll(sequence, st);
    
    /*
    //int[] starts = { 41104, 42405, 47865, 47866, 47867, 49662, 63308, 63309, 63310, 86643, 97851, 113315, 113316, 113317, 113318, 113319, 113333, 117435, 117436, 119528, 119839, 126582, 128163, 154715, 154716, 185363, 186491, 186492, 208631, 208632, 222353, 222819, 223376, 237833, 237834, 237835, 237836, 238651, 238652, 238653, 238654, 239469, 239470, 239471, 239472, 240273, 240402, 240403, 244193, 245358, 246114, 246175, 253778, 254597, 263224, 263225, 263360, 263979, 263980, 281267, 281467, 281468, 281848, 282246, 282446, 282447, 286627, 286636, 286658, 286659, 288836, 288866, 288867, 288868, 288869, 288870, 288871, 288872, 288873, 288874, 288875, 288876, 288877, 289698, 292845, 292846, 292851, 293082, 293892, 293893, 293898, 294669, 294670, 294675, 294906, 295646, 295647, 295652, 295883, 304603, 325705, 326159, 334473, 334474, 334475, 334476, 335413, 335414, 335419, 335420, 335421, 335422, 335423, 335424, 335439, 335440, 335441, 335442, 335443, 335450, 335451, 335452, 335453, 335472, 335473, 335474, 335500, 335505, 335506, 335507, 335508, 335602, 335603, 335608, 335609, 335610, 335611, 335612, 335613, 335628, 335629, 335630, 335631, 335632, 335639, 335640, 335641, 335642, 335661, 335662, 335663, 335689, 335694, 335695, 335696, 335697, 337295, 337296, 350621, 350668, 357607, 360723, 365918, 365965, 369015, 369815, 369816, 369817, 370715, 388137, 389524, 389812, 396175, 396176, 403956, 404483, 409263, 409320, 420383, 431607, 431614, 431619, 431620, 431627, 431628, 431629, 431630, 431637, 431638, 431645, 431646, 431688, 431718, 431729, 431730, 431731, 431732, 431745, 431746, 431747, 431748, 431749, 431750, 431751, 431752, 431753, 431754, 431755, 431756, 431757, 431758, 431759, 431760, 431761, 431762, 431763, 431764, 431774, 431775, 431781, 431789, 431790, 431791, 438917, 439227, 439228, 439229, 439230, 439231, 439232, 439233, 439234, 439235, 439236, 439237, 439238, 439239, 439240, 439241, 439242, 439243, 439244, 439245, 439246, 439247, 439248, 439254, 439255, 439256, 439257, 439258, 439259, 439260, 439261, 439262, 439263, 439264, 439277, 439278, 439279, 439280, 439281, 439282, 439283, 439284, 439285, 439286, 439287, 439288, 439289, 439302, 439307, 439308, 439309, 439310, 439311, 439312, 439352, 439353, 439358, 439363, 439364, 439379, 439380, 439381, 439382, 439383, 439384, 439399, 439400, 439401, 439402};
    int[] starts = {439312, 439352, 439353, 439358, 439363, 439364, 439379, 439380, 439381, 439382, 439383, 439384, 439399, 439400, 439401, 439402};
    GappedSuffixTree st = new GappedSuffixTree(sequence, true);
    for (int start : starts) {
      st.insert(start);
      //System.err.println(sequence.toString(start, start+300));
    }*/
    
    //System.err.print(st);
    
    /*
    HashSet<Node> seen = new HashSet<Node>();
    GappedSuffixTree.countNodes(st.getRoot(), seen);
    //System.err.println("Number of internal Nodes " + seen.size());
    //System.out.println("Making gap links");
    st.makeGapLinks();
    
    // [DR,3090] [VI,1287] [VDQ,330768] [DE,3085] [FVE,197901] [L,4] [LPE,264461] [K,2] [DE,3085] [EER,855314]
    int[] qEdgesInts = {3090, 1287, 330768, 3085, 197901, 4, 264461, 2, 3085, 855314};
    ArrayList<QueryEdge> query = new ArrayList<QueryEdge>();
    for (int i : qEdgesInts) {
      query.add(new QueryEdge(i));
    }
    */
    
    /*
    String queryString = "PILG";
    ArrayList<QueryEdge> query = new ArrayList<QueryEdge>();
    for (int i = 0; i<queryString.length(); i++) {
      query.add(new QueryEdge(sequence.toByte(queryString.charAt(i))));
    }
    */
    
    /*
    HashSet<Integer> matches = st.search(query);
    if (matches!=null) {
      System.err.println("Matches [" + matches.size() + "] at: ");
      for (int i : matches) {
        System.err.print(i + " ");
      }
      System.err.println();
    }
    else {
      System.err.println("Not found");
    }
    */
  }
}
