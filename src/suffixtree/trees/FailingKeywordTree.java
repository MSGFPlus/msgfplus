package suffixtree.trees;

import java.util.ArrayList;
import java.util.HashSet;

import sequences.ProteinFastaSequence;
import suffixtree.Constants;
import suffixtree.edges.Edge;
import suffixtree.edges.MassEdge;
import suffixtree.matches.ExactMatchObject;
import suffixtree.nodes.FailingNode;
import suffixtree.nodes.InternalNode;



/**
 * This is the keyword tree with failure links
 * @author jung
 *
 */
public class FailingKeywordTree {

  private FailingNode root;
  private MassEdge[] queries;
  
  /**
   * Constructor taking a list of list of integer masses.
   * @param queries the array of array of masses to initialize the Keyword tree.
   */
  public FailingKeywordTree(ArrayList<ArrayList<Integer>> queries) {
    // initialize
    this.root = new FailingNode();
    this.root.setFail(this.root);
    this.queries = new MassEdge[queries.size()];
    
    int count = 0;
    for (ArrayList<Integer> iArray : queries) {
      // create individual edges for each mass
      /*
      FailingNode currentNode = null;
      MassEdge edge = null;
      for (int mass : iArray) {
        int[] massArray = {mass};
        FailingNode sinkNode = new FailingNode();
        MassEdge subEdge = new MassEdge(massArray, sinkNode);
        if (currentNode!=null) currentNode.insert(subEdge, 0);
        else                   edge = subEdge;
        currentNode = sinkNode;
      }
      currentNode.addPosition(count);
      */
      this.queries[count] = new MassEdge(iArray, new FailingNode(count));
      //System.out.println("Inserting " + this.queries[count]);
      this.insert(this.queries[count].duplicate());
      if (++count%100000 == 0) {
        System.out.println("Tree built for " + count + " sequences");
      }
    }

    // convert this into a trie 
    System.out.println("Converting to a trie");
    breakAllEdges(this.root);

    System.out.println("Starting to make failure links");
    makeFailureLinks();
    
    // check that everything was inserted correctly
    /*
    for (ArrayList<Integer> query : queries) {
      FailingNode current = this.root;
      for (int mass : query) {
        int match = current.search(mass);
        if (match<0) {
          System.err.println("Not found!");
        }
        else {
          current = (FailingNode)current.getEdgeAt(match).getSink();
        }
      }
    }*/
  }

  
  private void breakAllEdges(FailingNode n) {
    for (int i=0; i<n.getDegree(); i++) {
      FailingNode fn = new FailingNode((InternalNode)n.getEdgeAt(i).getSink());
      n.getEdgeAt(i).setSink(fn);
      breakAllEdges(fn);
      ((MassEdge)n.getEdgeAt(i)).divideAll();
    }
  }
  
  /**
   * Helper method for the insertion of a mass edge.
   * @param edge the edge to insert to this object.
   */
  private void insert(MassEdge edge) {
    this.root.insert(edge);
  }
  
  
  /**
   * Use recursive method to build the failing keyword tree
   */
  private void makeFailureLinks() {
    
    ArrayList<FailingNode> currentStack = new ArrayList<FailingNode>();
    for (int i=0; i<this.root.getDegree(); i++) {
      FailingNode nextFn = new FailingNode(((InternalNode)this.root.getEdgeAt(i).getSink()));
      this.root.getEdgeAt(i).setSink(nextFn);
      nextFn.setFail(this.root);
      currentStack.add(nextFn);
    }
    
    //System.out.println("Root degree " + currentStack.size());
    
    while (currentStack.size() > 0) {
      //System.out.println("At level " + ++level + " Node count: " + currentStack.size());
      ArrayList<FailingNode> nextStack = new ArrayList<FailingNode>();
      for (FailingNode n : currentStack) {
        for (int i=0; i<n.getDegree(); i++) {
          Edge outEdge = n.getEdgeAt(i);
          FailingNode nextFn = new FailingNode((InternalNode)outEdge.getSink());
          outEdge.setSink(nextFn);
          FailingNode fn = (FailingNode)n.getFail();
          int matchIndex = fn.search(outEdge);
          if (matchIndex>=0) {
            nextFn.setFail((FailingNode)fn.getEdgeAt(matchIndex).getSink());
          }
          else {
            nextFn.setFail(this.root);
          }
          nextStack.add(nextFn);
        }
      }
      currentStack = nextStack;
    }
  }
  
  
  /**
   * Feed in the database and retrieve the matches
   * @param database the amino acid fasta database
   * @param results the set of match objects representing the matching
   */
  public void feed(ProteinFastaSequence sequence, HashSet<ExactMatchObject> matches) {
    
    // Make this a circular rotating list, so we don't have to recreate the hashSets
    ArrayList<HashSet<FailingNode>> hashList = new ArrayList<HashSet<FailingNode>>();
    
    int hashListIndex = 0;
    int hashListLength = 1000;
    for (int i=0; i<hashListLength; i++) hashList.add(new HashSet<FailingNode>());
    
    hashList.get(hashListIndex).add(this.root);
    
    //byte b;
    
    // cycle through all the start positions
    int nodeCount = 0;
    for (int start=0; start<sequence.getSize(); start++) {
      
      if (start%10000==3) {
        System.out.println("At " + start + ". Averaging: " + nodeCount / start + ". Matches: " + matches.size());
      }
      
      if (sequence.isTerminator(start)) {
        // clear all the hashes in the list
        do {
          hashList.get(hashListIndex).clear();
          hashListIndex = (hashListIndex+1) % hashListLength;
        } while(hashList.get(hashListIndex).size()>0);
        hashList.get(hashListIndex).add(this.root);
        continue;
      }
      
      //hashList.get(hashListIndex).add(this.root);
      
      int cumMass = 0;
      for (int next=start; next<sequence.getSize(); next++) {
        if (sequence.isTerminator(next)) break;
        cumMass += sequence.getIntegerMass(next);
        if (cumMass > Constants.MAX_GAP_MASS) break;
        
        int nextHashIndex = (hashListIndex+next-start+1)%hashListLength;
        nodeCount += hashList.get(hashListIndex).size();
        for (FailingNode fn : hashList.get(hashListIndex)) {
          
          // register the matches
          for (int pos : fn.getPositions()) {
            MassEdge query = this.queries[pos];
            //System.out.println("Position: " + next +" Match: " + query);
            
            // only reasonable way to trace back the starting position is to go back and accumulate the mass
            int totalMass = query.getTotalMass();
            int backMass = 0, begin = next-1;
            while (totalMass>backMass) {
              backMass += sequence.getIntegerMass(begin--); 
            }
            //System.out.println(sequence.toString(begin+1, next));
            matches.add(new ExactMatchObject(sequence, begin, next, query.getLabels(), pos));
          }
          
          int matchIndex = fn.search(cumMass);
          if (matchIndex >= 0) {
            // add the node to the next hashPosition
            hashList.get(nextHashIndex).add((FailingNode)fn.getEdgeAt(matchIndex).getSink());
          }
          else {
            // add the failing nodes until we reach the root
            FailingNode testNode = (FailingNode)fn.getFail();
            while (testNode!=this.root) {
              matchIndex = testNode.search(cumMass);
              if (matchIndex>=0) {
                testNode = (FailingNode)testNode.getEdgeAt(matchIndex).getSink();
                break;
              }
              else {
                testNode = (FailingNode)testNode.getFail();
              }
            }
            
            if (testNode!=this.root) {
              hashList.get(nextHashIndex).add(testNode);
            }
            else {
              int nextIndex = (hashListIndex+1) % hashListLength;;
              hashList.get(nextIndex).add(this.root);
            }
          }
        }
      }
      hashList.get(hashListIndex).clear();
      hashListIndex = (hashListIndex+1) % hashListLength;
    }
    
  }
}
