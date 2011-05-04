package suffixtree.trees;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Stack;
import java.util.TreeMap;

import msutil.AminoAcid;

import sequences.ProteinFastaSequence;
import suffixtree.Constants;
import suffixtree.Mutation;
import suffixtree.edges.DirectedMassEdge;
import suffixtree.edges.Edge;
import suffixtree.matches.MatchObject;
import suffixtree.matches.MutMatchObject;
import suffixtree.nodes.ComplexInternalNode;
import suffixtree.nodes.Node;


/**
 * Data structure for matching gapped peptides against a database implementing
 * the mutated search
 * @author jung
 *
 */
public class ComplexKeywordTree {
  
  
  private ComplexInternalNode root;         // the root of the forward tree
  private ComplexInternalNode leaf;         // the root of the reverse tree
  private DirectedMassEdge[] queries;       // the list of queries
  private ComplexInternalNode[] traceBacks; // the bridge from one tree to the other
  private ProteinFastaSequence db;          // the original sequence
  private int minMatchLength;               // for partial matches this is the minimum length to match (inclusive)
  
  private long currentStart;
  private int currentShift;
  //private int count = 0;
  
  

  /**
   * Constructor taking a list of list of integer masses.
   * @param queries the array of array of masses to initialize the Keyword tree.
   * @param database the protein database object
   */
  public ComplexKeywordTree(ArrayList<ArrayList<Integer>> queries, ProteinFastaSequence database) {
    // initialize
    this.root = new ComplexInternalNode();
    this.leaf = new ComplexInternalNode();
    this.queries = new DirectedMassEdge[queries.size()];
    this.traceBacks = new ComplexInternalNode[queries.size()]; 
    this.db = database;
    
    DirectedMassEdge[] rQueries = new DirectedMassEdge[queries.size()];
    
    int count = 0, step = 0, minQueryLength = Integer.MAX_VALUE;
    System.out.print("Building keyword trees: ");
    for (ArrayList<Integer> iArray : queries) {

      // compare minimum length entry
      if (iArray.size()<minQueryLength)     minQueryLength = iArray.size();
      
      // create the reverse version of the query 
      ArrayList<Integer> reversed = new ArrayList<Integer>(iArray);
      Collections.reverse(reversed);
      
      ComplexInternalNode rootNode = new ComplexInternalNode(count);
      rootNode.setParentNode(this.leaf);
      DirectedMassEdge rEdge = new DirectedMassEdge(reversed, rootNode);
      rQueries[count] = rEdge.duplicate();
      this.leaf.insert(rEdge);
      this.traceBacks[count] = rootNode;
      
      //System.out.println("Inserting " + edge);
      ComplexInternalNode leafNode = new ComplexInternalNode(count);
      leafNode.setParentNode(this.root);
      DirectedMassEdge edge = new DirectedMassEdge(iArray, leafNode);
      this.queries[count] = edge.duplicate();
      this.root.insert(edge);
      
      // display progress
      if (count++ >= (step*queries.size())/20.0) {
        step++;
        System.out.printf(" %d%%", step*5);
      }
    }
    System.out.println();

    // re-initialize the tracebacks once the full tree is constructed, reconnect the parents
    int queryIndex = 0;
    for (DirectedMassEdge edge : rQueries) {
      
      ComplexInternalNode n = this.leaf;
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
    
    
    this.minMatchLength = minQueryLength / 2;
    System.out.println("Minimum partial match size " + this.minMatchLength);
    
    // populate both trees
    this.populate();
    this.reversePopulate();
  }
  
  
  /**
   * Collect prefix-suffix matches on this tree and database. There minimum mass
   * value for a prefix and suffix anchor is specified in Constants class.
   * @return the list of matched objects
   */
  public ArrayList<MatchObject> collectPrefixSuffixMatches() {
    ArrayList<MatchObject> matches = new ArrayList<MatchObject>();
    collectPrefixSuffixMatches(this.root, new Stack<ComplexInternalNode>(), matches);
    return matches;
  }
  
  
  /**
   * Collect prefix-suffix matches on this tree and database with a mass bridge.
   * @return the list of mass objects.
   */
  public ArrayList<MutMatchObject> collectMatchesWithOneMutation() {
    //this.count = 0;
    ArrayList<MutMatchObject> matches = new ArrayList<MutMatchObject>();
    collectMatchesWithOneMutation(this.root, new TreeMap<Integer,ArrayList<Long>>(), matches);
    return matches;
  }
 
  
  /**
   * This is the helper method that matching the given mass by allowing a single
   * mutation. If the mass of the targetMass is equal to the mass of any amino
   * acids, the insertion mutation will not be added.
   * @param start the start coordinate in the database.
   * @param targetMass the mass to match
   * @param muts the resulting mutations
   */
  private void matchDbWithMutation(long start, int targetMass, ArrayList<Mutation> muts) {
    
    int currentCumMass = 0;
    
    // SPECIAL CASE: match target mass with an amino acids insertion to the DB
    if (targetMass <= Constants.MAX_AMINO_ACID_MASS) {
      for (AminoAcid aa : AminoAcid.getAminoAcids(targetMass)) {
        muts.add(new Mutation(start, start, Constants.EMPTY_AA, aa.getResidue()));
      }
    }
    
    // This is for optimization purposes. 
    // As we go through amino acids in the database, we keep track of the maximum seen mass. 
    // This will allow us to break once the cumulative database mass exceeds the target mass
    // by this maximum mass
    int maxMass = Integer.MIN_VALUE; // this is the upper bound for the database deletion
    
    int delta;
    
    // Try to extend until endPosition (inclusive)
    for (long endPosition=start; endPosition<this.db.getSize(); endPosition++) {
      
      // we are done if this is true
      if (!this.db.hasMass(endPosition)) break;
      
      int currentMass = this.db.getIntegerMass(endPosition);
      if (currentMass>maxMass) maxMass = currentMass;
      currentCumMass += currentMass;
      
      delta = currentCumMass - targetMass; 
      
      // no deletion can reconcile the masses difference
      if (delta > maxMass) break;
      
      // there is no mutation that reconcile this delta
      if (delta < -Constants.MAX_AMINO_ACID_MASS) continue;
      
      // do not allow mutations that do not change the mass
      if (delta == 0) continue;
      
      //System.out.println("Database mass " + currentCumMass + " delta " + delta);
      
      // we can try to mutate all the characters up to the deltaPosition
      for (long mutationPosition=start; mutationPosition<=endPosition; mutationPosition++) {
          
        char original = this.db.getCharAt(mutationPosition);
        ArrayList<Character> mutations = mutateTo(original, delta);
        if (mutations!=null) {
          // there are candidate mutations
          for (char mutation : mutations) {
            
            // SPECIAL CASE: Deleting the first amino in the db is the same as
            // deleting the last amino acid in the previous iteration, so we
            // do not add a duplicate to the matches
            if (!(mutationPosition==start && mutation==Constants.EMPTY_AA)) {
              muts.add(new Mutation(mutationPosition, endPosition+1, original, mutation));
            }
            else {
              //muts.add(new Mutation(mutationPosition, endPosition+1, original, mutation, true));
            }
          }
        }
      }
      
      // evaluate the possibility of an insertion
      for (AminoAcid aa : AminoAcid.getAminoAcids(-delta)) {
        //System.out.println("Adding " + aa.getResidue());
        muts.add(new Mutation(start, endPosition+1, Constants.EMPTY_AA, aa.getResidue()));  
      }
      
    }
  }
  
  
  /**
   * This is the helper method that matching the given mass by allowing a single
   * mutation. If the mass of the targetMass is equal to the mass of any amino
   * acids, the insertion mutation will not be added. The methods reverses the
   * coordinates and matches in reverse
   * @param start the start coordinate in the database (inclusive).
   * @param targetMass the mass to match
   * @param muts the resulting mutations
   */
  private void matchDbWithMutationR(long start, int targetMass, ArrayList<Mutation> muts) {
    int currentCumMass = 0;
    
    // SPECIAL CASE: match target mass with an amino acids insertion to the DB
    if (targetMass <= Constants.MAX_AMINO_ACID_MASS) {
      for (AminoAcid aa : AminoAcid.getAminoAcids(targetMass)) {
        muts.add(new Mutation(start+1, start, Constants.EMPTY_AA, aa.getResidue()));
      }
    }
    
    // This is for optimization purposes. 
    // As we go through amino acids in the database, we keep track of the maximum seen mass. 
    // This will allow us to break once the cumulative database mass exceeds the target mass
    // by this maximum mass
    int maxMass = Integer.MIN_VALUE; // this is the upper bound for the database deletion
    
    int delta;
    
    // Try to extend until endPosition (inclusive)
    for (long endPosition=start; endPosition>=0; endPosition--) {
      
      // we are done if this is true
      if (!this.db.hasMass(endPosition)) break;
      
      int currentMass = this.db.getIntegerMass(endPosition);
      if (currentMass>maxMass) maxMass = currentMass;
      currentCumMass += currentMass;
      
      delta = currentCumMass - targetMass; 
      
      // no deletion can reconcile the masses difference
      if (delta > maxMass) break;
      
      // there is no mutation that reconcile this delta
      if (delta < -Constants.MAX_AMINO_ACID_MASS) continue;
      
      // do not allow mutations that do not change the mass
      if (delta == 0) continue;
      
      //System.out.println("Database mass " + currentCumMass + " delta " + delta);
      
      // we can try to mutate all the characters up to the deltaPosition
      for (long mutationPosition=start; mutationPosition>=endPosition; mutationPosition--) {
          
        char original = this.db.getCharAt(mutationPosition);
        ArrayList<Character> mutations = mutateTo(original, delta);
        if (mutations!=null) {
          // there are candidate mutations
          for (char mutation : mutations) {
            
            // SPECIAL CASE: Deleting the first amino in the db is the same as
            // deleting the last amino acid in the previous iteration, so we
            // do not add a duplicate to the matches
            if (!(mutationPosition==endPosition && mutation==Constants.EMPTY_AA)) {
              if (mutationPosition==start && mutation==Constants.EMPTY_AA) {
                // this adds an extra matched edge to the results
                muts.add(new Mutation(mutationPosition, endPosition-1, original, mutation, true));
              }
              else {
                muts.add(new Mutation(mutationPosition, endPosition-1, original, mutation));
              }
            }
            else {
              //muts.add(new Mutation(mutationPosition, endPosition+1, original, mutation, true));
            }
            
            //muts.add(new Mutation(mutationPosition, endPosition-1, original, mutation));
          }
        }
      }
      
      // evaluate the possibility of an insertion
      for (AminoAcid aa : AminoAcid.getAminoAcids(-delta)) {
        //System.out.println("Adding " + aa.getResidue());
        muts.add(new Mutation(start, endPosition-1, Constants.EMPTY_AA, aa.getResidue()));  
      }
      
    }
  }
  
  
  
  /**
   * Tries to match the db with the given mass.
   * @param start the start position to start the matching (include this position)
   * @param mass the mass to match
   * @return the position of the next place to start matching next, so we can
   *         chain this function or -1 for no match. The returning position is
   *         <b>inclusive</b>.
   */
  private long matchDb(long start, int mass) {
    int cumMass = 0;
    for (long position = start; position < this.db.getSize(); position++) {
      if (this.db.isTerminator(position)) break;
      cumMass += this.db.getIntegerMass(position);
      if (cumMass==mass) return position;
    }
    return -1;
  }
  
  /**
   * Tries to match the db with the given mass in reverse.
   * @param start the start position to start the matching (include this position)
   * @param mass the mass to match
   * @return the position of the next place to start matching next, so we can
   *         chain this function or -1 for no match. The returning position
   *         is <b>inclusive</b>.
   */
  private long matchDbR(long start, int mass) {
    int cumMass = 0;
    for (long position = start; position >= 0; position--) {
      if (this.db.isTerminator(position)) break;
      cumMass += this.db.getIntegerMass(position);
      if (cumMass==mass) return position;
      if (cumMass>mass) return -1;
    }
    return -1;
  }
  
  
  /**
   * Exhaustively extend to the left allowing for a mutation in the first
   * match
   * @param queryIndex the index of the query in question
   * @param lastMatchedEdgeIndex the index of the last matched edge
   * @param start the start index of the partial match (right most inclusive)
   * @param end the end index of the partial match (left most inclusive)
   * @param matches the match object to store the matches
   */
  private void extendLeft(int queryIndex, long start, long end, ArrayList<MutMatchObject> matches) {

    DirectedMassEdge query = this.queries[queryIndex];
    
    // figure out the partial match coordinates on the edge
    int matchedMass = this.db.getIntegerMass(start, end);
    int matchedIndex = query.getIndexFromEnd(matchedMass);
    
    for (long position=start; position<end; ) {
      //System.out.printf("matchedIndex: %d .minMatchLength: %d\n", matchedIndex, this.minMatchLength);
      
      if (end==975625) {
        System.out.println("Attempting to do an extension on " + this.queries[queryIndex] + " from mass " + query.getLabelAt(matchedIndex));  
      }
      
      if (matchedIndex > this.minMatchLength) break;
      
      // do a mutation match of for the next edge
      ArrayList<Mutation> muts = new ArrayList<Mutation>();
      matchDbWithMutationR(position-1, query.getLabelAt(matchedIndex-1), muts);
      //System.out.println("Manually extend to the left starting with mass (mutating) " + query.getLabelAt(matchedIndex-1));
      
      if (end==975625) {
        System.out.println("Mutation candidates: " + muts.size());  
      }
      
      int massToShove = query.getLabelAt(matchedIndex);
      int cumMass = 0;
      while (cumMass<massToShove) {        
        cumMass += this.db.getIntegerMass(position++);
      }
      
      
      for (Mutation m : muts) {
        
        if (!m.isInEdge() && matchedIndex + 1 > this.minMatchLength) continue;
        
        // do a manual extension, exact match for the remaining edges
        int targetEdgeIndex = matchedIndex-2;

        //System.out.println(query.getLabelAt(matchEdgeIndex) + " " + this.db.getCharAt(m.nextStart));
        
        long nextStart = m.getNextStart();
        while (targetEdgeIndex>=0) {
          nextStart = matchDbR(nextStart, query.getLabelAt(targetEdgeIndex));
          if (nextStart >= 0) {
            targetEdgeIndex--;
            nextStart--;
          }
          else {
            break;
          }
        }
          
        if (targetEdgeIndex<=0 && nextStart>=0) {
          // do not allow a deletion of the last amino acid if no additional matches
          //if (!(targetEdgeIndex-1<=0 && m.mutation==Constants.EMPTY_AA)) {
            // create the match objects
            MutMatchObject smmo = new MutMatchObject(nextStart+1,
                                                                           end, 
                                                                           m,
                                                                           this.db, 
                                                                           query.getLabels(), 
                                                                           queryIndex);
            matches.add(smmo);
          //}
        }
      } 
              
      
      //System.out.printf("%d vs %d\n", cumMass, massToShove);
      matchedIndex++;
    }
    
    return;
    
  }

   

  /**
   * Helper method to collect prefix-suffix matches.
   * @param node the node to recurse on, in the prefix tree
   * @param prefixes the TreeMap for (mass->prefixCoordinates) pairs
   * @param matches store the matches in this object
   */
  private void collectMatchesWithOneMutation(ComplexInternalNode node,
                                             TreeMap<Integer,ArrayList<Long>> prefixes,
                                             ArrayList<MutMatchObject> matches) {
    
    // add all prefix matches collected by taking the current node
    int minPrefixMass = Integer.MAX_VALUE;  // register this value for tree popping at the end
    
    for (int i=0; i<node.getPrefixMatchCount(); i++) {
      long start = node.getPrefixStartAtIndex(i), end = start+node.getPrefixExtendAtIndex(i);
      int cumMass = this.db.getIntegerMass(start, end);
      if (!prefixes.containsKey(cumMass)) {
        // put the empty entry
        prefixes.put(cumMass, new ArrayList<Long>());  
      }
      prefixes.get(cumMass).add(node.getCoordinatesAtIndex(i));
     
      // register the minimum mass so we can pop the tree accordingly later
      if (minPrefixMass > cumMass) minPrefixMass = cumMass;
    }
    
    
    // base case because there is one or many query that ends in here
    if (node.getPositions().length>0) {
      
      HashMap<Integer,HashMap<Long,Long>> suffixes = new HashMap<Integer,HashMap<Long,Long>>();
      
      /**
       * All prefixes have been collected and stored in the TreeMap at this point
       **/
      
      // cycle through all the queries ending in this node (this is normally 1 query)  
      for (int queryIndex : node.getPositions()) {
        
        DirectedMassEdge query = this.queries[queryIndex];
        int queryTotalMass = query.getTotalMass();
        
        
        
        // collect all suffixes into a suitable data structure
        // do manual extension to the left if necessary 
        ComplexInternalNode suffixMatchNode = this.traceBacks[queryIndex];
        suffixes.clear();
        while (true) {
          
          // cycle through all the possible coordinates
          for (int matchIndex=0; matchIndex<suffixMatchNode.getPrefixMatchCount(); matchIndex++) {
            long start = suffixMatchNode.getPrefixStartAtIndex(matchIndex);
            long end = suffixMatchNode.getPrefixExtendAtIndex(matchIndex) + start;
            int suffixMass = this.db.getIntegerMass(start, end);
            if (suffixMass < queryTotalMass) {
              if (!suffixes.containsKey(suffixMass)) {
                suffixes.put(suffixMass, new HashMap<Long,Long>());
              }
              suffixes.get(suffixMass).put(start, end);
              
             
              if (end==709488) {
                System.out.println("Retrieving " + query);
                System.out.println("--- Partial suffix match " + this.db.getCharAt(start-1) + " " + start + "-" + this.db.getSubsequence(start, end) + "-" + end);
                System.out.println("Index " + matchIndex);
              }
              
              if (end==709488) {
                System.out.println("Query: " + query + ". MatchIndex " + matchIndex);  
              }
              
              // check whether a manual extension is required
              extendLeft(queryIndex, start, end, matches);
              
            }
          }
          suffixMatchNode = suffixMatchNode.getParentNode();
          if (suffixMatchNode==null) break;
        }
        
        
        ArrayList<Long> allCoors = new ArrayList<Long>();
        for (int prefixMass : prefixes.keySet()) {
          allCoors.addAll(prefixes.get(prefixMass));
        }
        
        // recursively get all ending position too from here
        node.retrieveAllCoordinates(allCoors);
        

        // build all prefix matches
        //for (int prefixMass : prefixes.keySet()) {
        //  for (long coors : prefixes.get(prefixMass)) {
        for (long coors : allCoors) {
            long start = ComplexInternalNode.decodeStartPosition(coors);
            long end = ComplexInternalNode.decodeExtension(coors)+start;
            
            if (start==975619) {
              System.out.println("Found the correct start. The query is " + query);  
            }
            
            int prefixCumMass = 0;
            int targetEdgeIndex = 0;
            int targetEdgeMass = query.getLabelAt(targetEdgeIndex++);
            for (long position=start; position<=end && prefixCumMass<queryTotalMass; position++) {
              if (prefixCumMass==targetEdgeMass) {
                
                //System.out.println("--- Partial match " + start + "-" + this.db.getSubsequence(start, position) + " " + this.db.getCharAt(position));
                
                ArrayList<Mutation> muts = new ArrayList<Mutation>();
                int currentEdgeMass = query.getLabelAt(targetEdgeIndex);
                int leftOverMass = queryTotalMass - prefixCumMass - currentEdgeMass;
  
                // match the next edge
                matchDbWithMutation(position, currentEdgeMass, muts);
                
                //System.out.printf("Examining %d position with letter %c edge %d\n", position, this.db.getCharAt(position), currentEdgeMass);
                for (Mutation m : muts) {
                  
                  //System.out.printf("Looking for %d mass after mutation from %c to %c\n", leftOverMass, mutation.original, mutation.mutation);
                  
                  int leftOverEdges = query.size() - targetEdgeIndex - 1;
                  if (!m.isInEdge()) leftOverEdges++; // this suffix edge was complete matched
                  long matchEnd = -1;
                  
                  if (leftOverEdges>=this.minMatchLength) {
                    
                    if (suffixes.containsKey(leftOverMass)) {
                      // do a suffix look up
                      HashMap<Long,Long> suffixCoors = suffixes.get(leftOverMass); 
                      if (suffixCoors.containsKey(m.getNextStart())) {
                        matchEnd = suffixCoors.get(m.getNextStart());
                      }
                    }
                  }
                  else {
                    
                    // do a manual extension
                    int matchEdgeIndex = targetEdgeIndex + 1;
                    long nextStart = m.getNextStart();
                    while (matchEdgeIndex<query.size()) {
                      nextStart = matchDb(nextStart, query.getLabelAt(matchEdgeIndex));
                      if (nextStart >= 0) {
                        matchEdgeIndex++;
                        nextStart++;
                      }
                      else {
                        break;
                      }
                    }
                      
                    if (matchEdgeIndex>=query.size() && nextStart>=0) {
                      // do not allow a deletion of the last amino acid if no additional matches
                      if (!(targetEdgeIndex+1>=query.size() && m.getMutation()==Constants.EMPTY_AA)) {
                        matchEnd = nextStart;
                      }
                    }
                  }
                    
                  // System.out.println(this.db.getSubsequence(suffixStart, suffixCoors.get(suffixStart)));  
                  if (matchEnd>=0) {  
                    // create the match objects
                    MutMatchObject smmo = new MutMatchObject(start,
                                                                                matchEnd, 
                                                                                m,
                                                                                this.db, 
                                                                                query.getLabels(), 
                                                                                queryIndex);
                    matches.add(smmo);
                  }
                    
                  //System.out.println(mutation);
                }
                
                // update the next target mass
                targetEdgeMass += query.getLabelAt(targetEdgeIndex++);
              }
              prefixCumMass += this.db.getIntegerMass(position);
            }
          //}
        }
      } // end of for loop
    } 
    
    // navigate the tree, by recursing on the next node
    for (int i=0; i<node.getDegree(); i++) {
      collectMatchesWithOneMutation((ComplexInternalNode)node.getEdgeAt(i).getSink(), prefixes, matches);
    }
      
    // pop the ends of the tree until the minimum mass is reached
    if (minPrefixMass < Integer.MAX_VALUE) {
      prefixes.tailMap(minPrefixMass, true).clear();
    }
    
  }
  
  
  
  /**
   * Helper method to collect prefix-suffix matches.
   * @param node the node to recurse on, in the prefix tree
   * @param path the stack of the path taken so far
   * @param matches store the matches in this object
   */
  private void collectPrefixSuffixMatches(ComplexInternalNode node, 
                                          Stack<ComplexInternalNode> path, 
                                          ArrayList<MatchObject> matches) {
    path.push(node);
    
    HashMap<Integer,Long> mass2position = new HashMap<Integer,Long>();
    
    // base case
    if (node.getPositions().length>0) {     // this node has a matching query
      
      for (int queryIndex : node.getPositions()) {
        //System.out.println("Retrieve reverse path from query index " + queryIndex);

        // Retrieve all suffix matches in the reverse tree and store them
        TreeMap<Integer,ArrayList<Long>> suffixes = new TreeMap<Integer,ArrayList<Long>>();
        ComplexInternalNode currentNode = this.traceBacks[queryIndex];
        int parentMass = this.queries[queryIndex].getTotalMass();
        while (true) {
          // collect all the suffixes
          for (int i=0; i<currentNode.getPrefixMatchCount(); i++) {
            long start = currentNode.getPrefixStartAtIndex(i), end = start+currentNode.getPrefixExtendAtIndex(i);
            int thisMass = this.db.getIntegerMass(start, end);
            if (thisMass!=parentMass) {     // exclude exact matches
              if (!suffixes.containsKey(thisMass)) suffixes.put(thisMass, new ArrayList<Long>());
              suffixes.get(thisMass).add(currentNode.getCoordinatesAtIndex(i));
            }
          }
          
          if (currentNode.getParentNode()==null) 
            // we have reached the root
            break;
          
          currentNode = currentNode.getParentNode();
        }
        
        /**
         * At this point all potential matches are stored in suffixes
         */
        
        //System.out.println("Prefixes...");
        for (ComplexInternalNode n : path) {
          for (int i=0; i<n.getPrefixMatchCount(); i++) {
            long start = n.getPrefixStartAtIndex(i), end = start+n.getPrefixExtendAtIndex(i);
            int prefixMass = this.db.getIntegerMass(start, end);
           
            if (prefixMass<parentMass) {     // exclude exact matches
              // build all possible matching scenarios, these matches break at the boundary of a gap
              /*
              for (int suffixMass : suffixes.tailMap(parentMass-prefixMass, true).keySet()) {
                for (long suffixCoors : suffixes.get(suffixMass)) {
                  // this match meets the cut
                  MatchObject m = new PrefixSuffixMatchObject(n.getCoordinatesAtIndex(i), 
                                                                     suffixCoors, 
                                                                     prefixMass, 
                                                                     suffixMass,
                                                                     this.db,
                                                                     this.queries[queryIndex].getLabels(),
                                                                     parentMass,
                                                                     queryIndex);
                  matches.add(m);
                }
              }
              */
              
              // account for those matches that do not break at a boundary
              int nextMass = this.queries[queryIndex].getMassAtIndexForward(this.queries[queryIndex].getIndexFromStart(prefixMass)+1);
              if (suffixes.containsKey(parentMass-nextMass)) {
                for (long suffixCoors : suffixes.get(parentMass-nextMass)) {
                  // we need to check that the items meet
                  int gapMass = nextMass - prefixMass;
                  //System.out.println("Gap mass is " + gapMass);
                  long suffixStart = ComplexInternalNode.decodeStartPosition(suffixCoors);
                  //long suffixEnd = ComplexInternalNode.decodeExtension(suffixCoors)+suffixStart;
                  mass2position.clear();
                  int cumMass = 0;
                  while (true) {
                    suffixStart--;
                    if (db.isTerminator(suffixStart)) break;
                    cumMass += db.getIntegerMass(suffixStart);
                    if (cumMass >= gapMass) break;
                    mass2position.put(cumMass, suffixStart);
                  } // all possible suffix masses and their starting positions
                  
                  cumMass = 0;
                  long prefixEnd = end;
                  while (true) {
                    if (db.isTerminator(prefixEnd)) break;
                    cumMass += db.getIntegerMass(prefixEnd++);
                    if (cumMass >= gapMass) break;
                    if (mass2position.containsKey(gapMass-cumMass)) {
                      suffixStart = mass2position.get(gapMass-cumMass);
                      /*
                      MatchObject m = new PrefixSuffixMatchObject(ComplexInternalNode.encodePositions(start, prefixEnd), 
                                                                         ComplexInternalNode.encodePositions(suffixStart, suffixEnd), 
                                                                         prefixMass+cumMass, 
                                                                         (parentMass-nextMass)+(gapMass-cumMass),
                                                                         this.db,
                                                                         this.queries[queryIndex].getLabels(),
                                                                         parentMass,
                                                                         queryIndex);
                      matches.add(m);
                      */
                      /*
                      System.out.println(m);
                      if (!m.isCorrect()) {
                        System.out.println("Incorrect");
                        System.out.printf("%d\t%d\t%d\t%d\n", prefixMass+cumMass, (parentMass-nextMass)+(gapMass-cumMass), parentMass, prefixEnd);
                        System.out.println(m);
                        System.exit(-9);
                      }*/
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    // navigate the tree, by recursing on the next node
    for (int i=0; i<node.getDegree(); i++) {
      collectPrefixSuffixMatches((ComplexInternalNode)node.getEdgeAt(i).getSink(), path, matches);
    }
    
    path.pop();
  }
  
  
  /**
   * Reverse populating the tree
   */
  private void reversePopulate() {
    
    System.out.print("Reverse feeding of the tree ");
    long compCount = 0, searchCount = 0, nextStop = 0;
    ArrayList<Integer> dbMasses = new ArrayList<Integer>();
    for (long start=this.db.getSize()-1; start>=0;) {
      
      if (this.db.isTerminator(start)) {
        // shift to the next protein
        dbMasses.clear();
        start--;
        continue;
      }
      
      if (dbMasses.size()==0) {
        
        // create a new query
        dbMasses.add(this.db.getIntegerMass(start));
        long end;     // find the end
        for (end=start-1; end>=0; end--) {
          if (this.db.isTerminator(end)) break;
          dbMasses.add(this.db.getIntegerMass(end));
        }  
        
        // To avoid creating the query again and again for each protein,
        // reuse the query and use a shift index instead
        for (long subStart=end; subStart<start; subStart++) {
          int shift = (int)(subStart-end);
          //System.out.println("Searching with shift: " + shift);
          this.currentStart = start-shift+1;
          this.currentShift = shift;
          
          int cumMass = 0;
          for (int i=shift; i<dbMasses.size(); i++) {
            cumMass += dbMasses.get(i);
            
            // we have reached an impossible path
            if (cumMass > this.leaf.getMaximumEdge().getLabel()) break;

            int matchIndex = this.leaf.search(cumMass);
            compCount++;
            if (matchIndex>=0) {
              compCount += this.storeMatches(dbMasses, i+1, this.leaf.getEdgeAt(matchIndex), 1, 1, true);
            }
          }
          searchCount++;
        }
        start = end;
      }
      
      int percDone = (int)(100-(start*100/this.db.getSize()));
      if (percDone>=nextStop) {
        nextStop += 5;
        System.out.printf(" %d%%", percDone);
      }
      
    }
    System.out.println("\nAverage number of comparisons per position " + compCount / searchCount);
 
  }
  
 
  
  /**
   * Register all partial matches into the keyword tree
   */
  private void populate() {
  
    int rootMaxEdge = this.root.getMaximumEdge().getLabel();
    
    System.out.print("Initializing the prefix tree ");
    long compCount = 0, searchCount = 0, nextStop = 0;
    ArrayList<Integer> dbMasses = new ArrayList<Integer>();
    for (long start=0; start<this.db.getSize();) {
      
      if (!this.db.hasMass(start)) {
        // shift to the next protein or segment
        dbMasses.clear();
        start++;
        continue;
      }
      
      if (dbMasses.size()==0) {
        
        // create a new query as an arrayList of masses
        dbMasses.add(this.db.getIntegerMass(start));
        long end;
        for (end=start+1; end<this.db.getSize(); end++) {
          if (!this.db.hasMass(end)) break;
          dbMasses.add(this.db.getIntegerMass(end));
        }  
     
        // To avoid translating the db into masses again and again, translate
        // everything into an array and use a shift index to indicate where to start
        for (long subStart=start; subStart<end; subStart++) {
          int shift = (int)(subStart-start);
          //System.out.println("Searching with shift: " + shift);
          this.currentStart = subStart;
          this.currentShift = shift;
          
          //System.out.println("Looking up " + database.toString(subStart, end));
          int cumMass = 0;
          // try to match the first edge here
          for (int i=shift; i<dbMasses.size(); i++) {
            cumMass += dbMasses.get(i);
            
            // we have reached an impossible path
            if (cumMass > rootMaxEdge) break;

            int matchIndex = this.root.search(cumMass);
            compCount++;
            if (matchIndex>=0) {
              compCount += this.storeMatches(dbMasses, i+1, this.root.getEdgeAt(matchIndex), 1, 1, false);
            }
          }
          
          searchCount++;
        }
        start = end;
      }
      
      int percDone = (int)(start*100/this.db.getSize());
      if (percDone>=nextStop) {
        nextStop += 5;
        System.out.printf(" %d%%", percDone);
      }
     
    }
    System.out.println("\nAverage number of comparisons per position " + compCount / searchCount);

  }

  
  /**
   * Helper method to match an edge against an edge.
   * @param dbMasses the masses of the database
   * @param startIndex the next index to try to match
   * @param currentMatchingEdge the matching edge
   * @param currentEdgeIndex the next index of the matching edge. This item could be 
   *                  equal to e.size(), which means, we have matched the edge
   *                  completely, and query.getLabelAt(query) should be matched
   *                  to an outgoing edge of e.getSink()
   * @param matchedMass the matched mass so far
   * @param reversed indicates whether the coordinates are reversed
   * @return the number of comparisons done for this query.
   */
  private int storeMatches(ArrayList<Integer> dbMasses, 
                           int startIndex, 
                           Edge currentMatchingEdge, 
                           int currentEdgeIndex,
                           int matchedEdgeCount,
                           boolean reversed) {
    int compCount = 0;
    
    ComplexInternalNode lastNode = null;    // the node to store the match or prefixes
    boolean hasMatch = false;               // flag to see whether there was match
    
    // We have reached the leaf of this branch OR we ran out of database
    if ((currentEdgeIndex==currentMatchingEdge.size() && currentMatchingEdge.getSink().getDegree()==0) 
        ||
        startIndex >= dbMasses.size()) {
      // base case
      lastNode = (ComplexInternalNode)currentMatchingEdge.getSink();
    }
    else {
      
      int cumMass = 0;
      if (currentEdgeIndex==currentMatchingEdge.size()) {
        // we have consumed a complete edge from the compressed tree, and we 
        // arrived to a node, so we need to find a matching edge from the node
        Node sink = currentMatchingEdge.getSink();
        // optimization stuff
        int lower = sink.getMinimumEdge().getLabel();
        int upper = sink.getMaximumEdge().getLabel();
      
        // branch many possible gaps
        for (int i=startIndex; i<dbMasses.size(); i++) {
          cumMass += dbMasses.get(i);
  
          // optimization
          if (cumMass < lower) continue;
          if (cumMass > upper) break;
          
          int matchIndex = sink.search(cumMass);
          compCount++;
          if (matchIndex>=0) {
            // recurse
            compCount += storeMatches(dbMasses, i+1, sink.getEdgeAt(matchIndex), 1, matchedEdgeCount+1, reversed);
            hasMatch = true;
          }
        }
        lastNode = (ComplexInternalNode)sink;
      }
      else {
        int currentMass = currentMatchingEdge.getLabelAt(currentEdgeIndex);
        for (int i=startIndex; i<dbMasses.size(); i++) {
          cumMass += dbMasses.get(i);
          
          if (cumMass < currentMass) continue;
          if (cumMass > currentMass) {
            // we have reached an impossible path, by over shooting
            break;
          }
          
          // There is a match
          compCount += storeMatches(dbMasses, i+1, currentMatchingEdge, currentEdgeIndex+1, matchedEdgeCount+1, reversed);
          hasMatch = true;
        }
        //lastNode = ((ComplexInternalNode)currentMatchingEdge.getSink()).getParentNode();
        lastNode = ((ComplexInternalNode)currentMatchingEdge.getSink());
      }
      
    }
    
    long start = this.currentStart-(startIndex-this.currentShift);
    long end = start+startIndex-this.currentShift;
    if (end==975625) {
      System.out.println("End 975625 found. Extends " + (end-start) + ". " + this.db.getSubsequence(start, end) + ". Match count " + matchedEdgeCount);
    }
    
    // add only those with a the minimum number of matched edges
    if (!hasMatch && matchedEdgeCount >= this.minMatchLength) {
      if (!reversed) {
        if (this.currentStart==975619) {
          System.out.println("Start 975619 found. Extends " + (startIndex-this.currentShift) + ". " + this.db.getSubsequence(this.currentStart, this.currentStart+startIndex-this.currentShift));
        }
        //System.out.println("Added prefix " + this.db.toString(this.currentStart, queryIndex));
        lastNode.addPartialMatch(this.currentStart, startIndex-this.currentShift);
      }
      else {
        
        if (end==975625) {
          System.out.println("Storage granted. Match count " + matchedEdgeCount + ". Last matched mass " + currentMatchingEdge.getLabel());
        }
        
        //System.out.println("Added " + this.db.toString(this.currentStart-(queryIndex-this.currentShift), this.currentStart));
        lastNode.addPartialMatch(this.currentStart-(startIndex-this.currentShift), startIndex-this.currentShift);
      }
    }
    return compCount;
  }
  
  
  @Override
  public String toString() {
    StringBuffer sb = new StringBuffer();
    for (int i=0; i<root.getDegree(); i++) {
      sb.append(root.getEdgeAt(i)+"\n");
    }
    return sb.toString();
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
  public static ArrayList<Character> mutateTo(char from, int mass) {
    if (mutTable.get(from).containsKey(mass)) return mutTable.get(from).get(mass);
    return null;
  }
  
  private static HashMap<Character,HashMap<Integer,ArrayList<Character>>> mutTable;
  static {
    mutTable = new HashMap<Character,HashMap<Integer,ArrayList<Character>>>();
    AminoAcid[] aas = AminoAcid.getStandardAminoAcids();
    for (AminoAcid aa : aas) {
      mutTable.put(aa.getResidue(), new HashMap<Integer,ArrayList<Character>>());
      for (AminoAcid mutant : aas) {
        if (mutant.getResidueStr()!=aa.getResidueStr()) {
          int massDiff = aa.getNominalMass() - mutant.getNominalMass();
          if (!mutTable.get(aa.getResidueStr()).containsKey(massDiff)) {
            mutTable.get(aa.getResidueStr()).put(massDiff, new ArrayList<Character>());
          }
          mutTable.get(aa.getResidueStr()).get(massDiff).add(mutant.getResidue());
        }
      }
      // add the deletion possibility 
      if (!mutTable.get(aa.getResidueStr()).containsKey(aa.getNominalMass())) {
        mutTable.get(aa.getResidueStr()).put(aa.getNominalMass(), new ArrayList<Character>());
      }
      mutTable.get(aa.getResidueStr()).get(aa.getNominalMass()).add(Constants.EMPTY_AA);
    }
  }
}
