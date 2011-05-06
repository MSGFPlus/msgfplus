package suffixtree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import msutil.AminoAcid;
import msutil.AminoAcidSet;

import sequences.MassSequence;
import sequences.ProteinFastaSequence;


/**
 * This class provides the static methods to match the database against masses.
 * @author jung
 *
 */
public class Matching {
  
  private static AminoAcidSet alpha = Constants.AA;
  private static final char[] EMPTY_CHAR_ARRAY = new char[0];
  
  // Variable that stores possible mutations from an amino acid to another, given a mass
  private static HashMap<Character,HashMap<Integer,char[]>> mutTable;
  static {
    HashMap<Character,HashMap<Integer,ArrayList<Character>>> mutTableList = new HashMap<Character,HashMap<Integer,ArrayList<Character>>>();
    AminoAcidSet aas = alpha;
    
    Iterator<AminoAcid> it1 = aas.iterator();
    while (it1.hasNext()) {
      AminoAcid aa = it1.next();
      mutTableList.put(aa.getResidue(), new HashMap<Integer,ArrayList<Character>>());
      Iterator<AminoAcid> it2 = aas.iterator();
      while (it2.hasNext()) {
        AminoAcid mutant = it2.next();
        if (mutant.getResidueStr()!=aa.getResidueStr()) {
          int massDiff = aa.getNominalMass() - mutant.getNominalMass();
          if (!mutTableList.get(aa.getResidueStr()).containsKey(massDiff)) {
            mutTableList.get(aa.getResidueStr()).put(massDiff, new ArrayList<Character>());
          }
          mutTableList.get(aa.getResidueStr()).get(massDiff).add(mutant.getResidue()); 
        }
      }
      // add the deletion possibility 
      if (!mutTableList.get(aa.getResidueStr()).containsKey(aa.getNominalMass())) {
        mutTableList.get(aa.getResidueStr()).put(aa.getNominalMass(), new ArrayList<Character>());
      } 
      mutTableList.get(aa.getResidueStr()).get(aa.getNominalMass()).add(Constants.EMPTY_AA); 
    }
    
    // convert the ArrayList into arrays
    mutTable = new HashMap<Character,HashMap<Integer,char[]>>();
    for (char aa : mutTableList.keySet()) {
      HashMap<Integer,char[]> newDeltaTable = new HashMap<Integer,char[]>();
      HashMap<Integer,ArrayList<Character>> deltaTable = mutTableList.get(aa);
      for (int delta : deltaTable.keySet()) {
        char[] mutations = new char[deltaTable.get(delta).size()];
        for (int i=0; i<mutations.length; i++) {
          mutations[i] = deltaTable.get(delta).get(i);
        }
        newDeltaTable.put(delta, mutations);
      }
      mutTable.put(aa, newDeltaTable);
    }
    
  }
  
  
  
  /**
   * Tries to match the db with the given mass.
   * @param db the sequence of amino acids
   * @param start the start position to start the matching (include this position)
   * @param mass the mass to match
   * @return the position of the next place to start matching next, so we can
   *         chain this function or -1 for no match. The returning position is
   *         <b>exclusive</b> i.e. the next index to start matching.
   */
  public static long matchDb(ProteinFastaSequence db, long start, int mass) {
    int cumMass = 0;
    for (long position = start; position < db.getSize() && cumMass<=mass; position++) {
      if (!db.hasMass(position)) break;
      cumMass += db.getIntegerMass(position);
      if (cumMass==mass) return position+1;
    }
    return -1;
  }
  
  
  
  /**
   * Tries to match the db with the given mass in reverse.
   * @param db the sequence of amino acids
   * @param start the start position to start the matching (include this position)
   * @param mass the mass to match
   * @return the position of the next place to start matching next, so we can
   *         chain this function or -2 for no match. The returning position
   *         is <b>exclusive</b> i.e. the next index to start matching.
   */
  public static long matchDbR(ProteinFastaSequence db, long start, int mass) {
    int cumMass = 0;
    for (long position = start; position >= 0 && cumMass<=mass; position--) {
      if (!db.hasMass(position)) break;
      cumMass += db.getIntegerMass(position);
      if (cumMass==mass) return position-1;
    }
    return -2;
  }
  
  
  
  /**
   * This is the helper method that matching the given mass by allowing a single
   * mutation. If the mass of the targetMass is equal to the mass of any amino
   * acids, the insertion mutation will not be added.
   * @param idbb the array of integer masses
   * @param cdbb the array of characters
   * @param dbbIndex the start coordinate in the array (inclusive).
   * @param targetMass the mass to match
   * @param muts the resulting mutations
   */
  public static void matchDbWithMutation(Integer[] idbb, 
                                         Character[] cdbb, 
                                         long dbStart,
                                         int dbbIndex,// item NOT matched
                                         int targetMass, 
                                         ArrayList<Mutation> muts) {
    
    int currentCumMass = 0;
    
    // SPECIAL CASE: match target mass with an amino acids insertion to the DB
    if (targetMass <= Constants.MAX_AMINO_ACID_MASS) {
      for (AminoAcid aa : alpha.getAminoAcids(targetMass)) {
        muts.add(new Mutation(dbStart+dbbIndex, dbbIndex, Constants.EMPTY_AA, aa.getResidue()));
      }
    }
    
    // This is for optimization purposes. 
    // As we go through amino acids in the database, we keep track of the maximum seen mass. 
    // This will allow us to break once the cumulative database mass exceeds the target mass
    // by this maximum mass
    int maxMass = Integer.MIN_VALUE; // this is the upper bound for the database deletion
    
    int delta;
    
    // Try to extend until endPosition (inclusive)
    for (int endPosition=dbbIndex; endPosition<idbb.length; endPosition++) {
      
      int currentMass = idbb[endPosition];
      if (currentMass>maxMass) maxMass = currentMass;
      currentCumMass += currentMass;
      
      delta = currentCumMass - targetMass; 
      
      // no deletion can reconcile the masses difference
      if (delta > maxMass) break;
      
      // there is no insertion that reconcile this delta
      if (delta < -alpha.getMaxNominalMass()) continue;
      
      // do not allow mutations that do not change the mass
      if (delta == 0) continue;
      
      //System.out.println("Database mass " + currentCumMass + " delta " + delta);
      
      // we can try to mutate all the characters up to the endPosition
      for (int mutationPosition=dbbIndex; mutationPosition<=endPosition; mutationPosition++) {
          
        char original = cdbb[mutationPosition];
        // there are candidate mutations
        for (char mutation : mutateTo(original, delta)) {
            
          // SPECIAL CASE: Deleting the first amino in the db is the same as
          // deleting the last amino acid in the previous iteration, so we
          // do not add a duplicate to the matches. We don't add a deletion
          // on the last amino acid
          if (!(mutationPosition==endPosition && mutation==Constants.EMPTY_AA)) {
            muts.add(new Mutation(dbStart+mutationPosition, endPosition+1, original, mutation));
          }
          // technically if mutationPosition==start && mutation==Constants.EMPTY_AA
          // we have an out of edge mutation, but for the forward matching we don't
          // care about this
        }
      }
      
      // evaluate the possibility of an insertion
      for (AminoAcid aa : alpha.getAminoAcids(-delta)) {
        muts.add(new Mutation(dbStart+dbbIndex, endPosition+1, Constants.EMPTY_AA, aa.getResidue())); 
      }
    }
  }
  
  
  
  /**
   * This is the helper method that matching the given mass by allowing a single
   * mutation. If the mass of the targetMass is equal to the mass of any amino
   * acids, the insertion mutation will not be added.
   * @param idbb the array of integer masses
   * @param cdbb the array of characters
   * @param dbbIndex the start coordinate in the array (inclusive).
   * @param targetMass the mass to match
   * @param muts the resulting mutations
   */
  public static void matchDbWithMutationR(Integer[] idbb, 
                                          Character[] cdbb, 
                                          long dbStart,
                                          int dbbIndex,// item NOT matched
                                          int targetMass, 
                                          ArrayList<Mutation> muts) {
    
    int currentCumMass = 0;
    
    // SPECIAL CASE: match target mass with an amino acids insertion to the DB
    if (targetMass <= Constants.MAX_AMINO_ACID_MASS) {
      for (AminoAcid aa : alpha.getAminoAcids(targetMass)) {
        muts.add(new Mutation(dbStart-dbbIndex+1, dbbIndex, Constants.EMPTY_AA, aa.getResidue()));
      }
    }
    
    // This is for optimization purposes. 
    // As we go through amino acids in the database, we keep track of the maximum seen mass. 
    // This will allow us to break once the cumulative database mass exceeds the target mass
    // by this maximum mass
    int maxMass = Integer.MIN_VALUE; // this is the upper bound for the database deletion
    
    int delta;
    
    // Try to extend until endPosition (inclusive)
    for (int endPosition=dbbIndex; endPosition<idbb.length; endPosition++) {
      
      int currentMass = idbb[endPosition];
      if (currentMass>maxMass) maxMass = currentMass;
      currentCumMass += currentMass;
      
      delta = currentCumMass - targetMass; 
      
      // no deletion can reconcile the masses difference
      if (delta > maxMass) break;
      
      // there is no insertion that reconcile this delta
      if (delta < -alpha.getMaxNominalMass()) continue;
      
      // do not allow mutations that do not change the mass
      if (delta == 0) continue;
      
      //System.out.println("Database mass " + currentCumMass + " delta " + delta);
      
      // we can try to mutate all the characters up to the endPosition
      for (int mutationPosition=dbbIndex; mutationPosition<=endPosition; mutationPosition++) {
          
        char original = cdbb[mutationPosition];
        // there are candidate mutations
        for (char mutation : mutateTo(original, delta)) {
            
          // SPECIAL CASE: Deleting the first amino in the db is the same as
          // deleting the last amino acid in the previous iteration, so we
          // do not add a duplicate to the matches. We don't add a deletion
          // on the last amino acid
          if (!(mutationPosition==endPosition && mutation==Constants.EMPTY_AA)) {
            muts.add(new Mutation(dbStart-mutationPosition, endPosition+1, original, mutation));
          }
          // technically if mutationPosition==start && mutation==Constants.EMPTY_AA
          // we have an out of edge mutation, but for the forward matching we don't
          // care about this
        }
      }
      
      // evaluate the possibility of an insertion
      for (AminoAcid aa : alpha.getAminoAcids(-delta)) {
        muts.add(new Mutation(dbStart-dbbIndex+1, endPosition+1, Constants.EMPTY_AA, aa.getResidue())); 
      }
    }
  }
 
  
  
  /**
   * This function allows to query amino acid transformation that yield a mass
   * delta. For example, G -> A will yield a -14 delta. So inputting G, -14 will
   * return A as the sole element of the ArrayList.
   * @param from the mass to mutate. '*' is accepted as the empty character.
   * @param mass the mass to add to the from mass
   * @return the list of masses resulting after such transformation. '*' 
   *         represents the empty character. null is returned if no mutation is possible.
   */
  public static char[] mutateTo(char from, int mass) {
    if (mutTable.get(from).containsKey(mass)) return mutTable.get(from).get(mass);
    return EMPTY_CHAR_ARRAY;
  }
  
  
  
  /**
   * This is the helper method that matching the given mass by allowing a single
   * modification of any mass.
   * @param db the sequence of amino acids
   * @param dbStart the start coordinate in the database.
   * @param dbb the db buffer
   * @param dbbIndex the index of the item to match
   * @param targetMass the mass to match
   * @param edgeIndex the index of the edge with the target mass in the query
   * @param mods the resulting modifications
   */
  public static void matchDbWithModification(MassSequence db, 
                                             long dbStart, 
                                             Integer[] dbb,
                                             int dbbIndex,
                                             int targetMass, 
                                             ArrayList<Modification> mods) {
    
    int currentCumMass = 0;
    int delta;
    
    // Try to extend until endPosition (inclusive)
    for (int endPosition=dbbIndex; endPosition<dbb.length; endPosition++) {

      int currentMass = dbb[endPosition];
      currentCumMass += currentMass;
      
      delta = currentCumMass - targetMass; 
      
      // there is no mutation that reconcile this delta
      if (delta < -Constants.MAX_MOD) continue;
      
      // no need to consider anything else
      if (delta > -Constants.MIN_MOD) break;
      
      // do not allow mutations that do not change the mass
      if (delta == 0) continue;
      
      mods.add(new Modification(endPosition+1, -delta, dbStart+dbbIndex, dbStart+endPosition+1));
    }
  }
 
  
  
  /**
   * This is the helper method that matching the given mass by allowing a single
   * mutation. If the mass of the targetMass is equal to the mass of any amino
   * acids, the insertion mutation will not be added. The methods reverses the
   * coordinates and matches in reverse
   * @param db the sequence of amino acids
   * @param dbStart the start coordinate in the database (inclusive).
   * @param dbb the db buffer
   * @param dbbIndex the index of the item to match
   * @param targetMass the mass to match
   * @param mods the resulting modifications
   */
  public static void matchDbWithModificationR(MassSequence db, 
                                              long dbStart,
                                              Integer[] dbb,
                                              int dbbIndex, 
                                              int targetMass, 
                                              ArrayList<Modification> mods) {

    int currentCumMass = 0;
    int delta;
    
    // Try to extend until endPosition (inclusive)
    for (int endPosition=dbbIndex; endPosition<dbb.length; endPosition++) {
      
      int currentMass = dbb[endPosition];
      currentCumMass += currentMass;
      
      delta = currentCumMass - targetMass; 
      
      // there is no mutation that reconcile this delta
      if (delta < -Constants.MAX_MOD) continue;
      
      // no deletion can reconcile the masses difference
      if (delta > -Constants.MIN_MOD) break;
      
      // do not allow mutations that do not change the mass
      if (delta == 0) continue;
     
      mods.add(new Modification(endPosition+1, -delta, dbStart-endPosition, dbStart-dbbIndex+1));
    }
  }
  
}
