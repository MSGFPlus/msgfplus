package suffixtree;

import msutil.AminoAcidSet;

public class Constants {
  
  // the maximum number of characters to lump together in the data structure
  public final static int MAX_GAP = 4;
  
  // the minimum number of edges to lump together
  //public final static int MIN_PARTIAL_MATCH = 2;
  
  // maximum size of the gap
  public final static int MAX_GAP_MASS = 500;
  
  // the maximum queriable mass
  public final static int MAX_QUERY_MASS = 3000;

  // the minimum queriable mass
  public final static int MIN_QUERY_MASS = 750;
  
  // minimum queriable character count
  public static final int MIN_QUERY_CHAR = 20;
  
  // maximum queriable character count
  public static final int MAX_QUERY_CHAR = 50;
  
  // required minimum mass to start storing partial matches
  //public static final int MIN_PARTIAL_MASS_MATCH = 500;
  
  //required minimum mass to start storing partial matches
  //public static final int MIN_PARTIAL_EDGE_MATCH = 3;
  
  // the smallest amino acid mass possible
  public static final int MIN_AMINO_ACID_MASS = 57;
  
  // the largest amino acid mass
  public static final int MAX_AMINO_ACID_MASS = 186;
  
  // the constant that represents an empty amino acid
  public static final char EMPTY_AA = '*';
  
  // when in prefix-suffix mode do not look for partial matches separated by more than this in the db
  public static final int MAX_GAP_MUTATED_SEARCH = 10;
  
  public static final int MAX_QUERY_BUNDLING_COUNT = 5000000;
  
  // the maximum probability to accept
  public static final float PROB_CUTOFF = 1E-10f;
  
  //the maximum probability to accept
  public static final float MUT_PROB_CUTOFF = 1E-10f;
  
  // the minimum modification mass, inclusive
  public static final int MIN_MOD = -50;
  
  // the maximum modification mass, inclusive
  public static final int MAX_MOD = 200;
  
  // the amino acid set to use
  //public static AminoAcidSet AA = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
  public static AminoAcidSet AA = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCysWithTerm();

}