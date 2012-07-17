package sequences;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Set;

public class FastaSequences implements Sequence {

  
  // the (path) name of the read files
  private ArrayList<String> files;
  
  // the end positions for each sequence (exclusive)
  private ArrayList<Long> positions;
  
  // the sequences, in case we need random access
  private ArrayList<FastaSequence> sequences;
  
  // the sequence currently loaded in memory, for sequencial access
  private FastaSequence current;
  private int currentIndex;
  
  // the alphabet specification
  private String aaSpec;
  
  private int id;
  
  private static final String metafileName = "sequences.ginfo";
  
  
  
  /**
   * Constructor create an object using the standard 20 amino acids (18 unique masses)
   * @param directory the directory where the fasta files are located
   * @param randomAccess flag to indicate loading of all sequences into memory
   */
  public FastaSequences(String directory, boolean randomAccess) {
    this(directory, Constants.AMINO_ACIDS_18, randomAccess);
  }
  
  
  /**
   * Constructor using an specific amino acid alphabet specification
   * @param directory the directory of the fasta files
   * @param aaSpec the amino acid alphabet specification
   * @param randomAccess flag to indicate loading of all sequences into memory
   */
  @SuppressWarnings("unchecked")
  public FastaSequences(String directory, String aaSpec, boolean randomAccess) {
    
    File dir = new File(directory);

    this.aaSpec = aaSpec;
    
    if (randomAccess) {
     // load all the sequences
      this.sequences = new ArrayList<FastaSequence>();  
    }
    
    // check whether the meta file exists
    if (new File(dir, metafileName).exists()) {
      // read the initialization parameters
      try {
        ObjectInputStream in = new ObjectInputStream(new FileInputStream(new File(dir, metafileName).getPath()));
        files = (ArrayList<String>)in.readObject();
        positions = (ArrayList<Long>)in.readObject();
        in.close();
      }
      catch (ClassNotFoundException e) {
        e.printStackTrace();
      }
      catch(FileNotFoundException e) {
        e.printStackTrace();
      } 
      catch (IOException e) {
        e.printStackTrace();
      }
      
      if (randomAccess) {
        for (String fileName : this.files) {
          sequences.add(new ProteinFastaSequence(fileName, aaSpec));
        }
      }
    }
    else {
      this.files = new ArrayList<String>();
      this.positions = new ArrayList<Long>();
      long cumPos = 0;
      // initialize the files and positions
      for (String file : dir.list()) {
        if (file.endsWith(".fasta")) {
          ProteinFastaSequence seq = new ProteinFastaSequence(new File(dir, file).getPath(), aaSpec);
          cumPos += seq.getSize();
          System.out.println("Loaded " + file);
          files.add(new File(dir, file).getPath());
          positions.add(cumPos);
          
          if (randomAccess) {
            sequences.add(seq);
          }
        }
      }
      // write the items to the file
      try {
        ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(new File(dir, metafileName).getPath()));
        out.writeObject(this.files);
        out.writeObject(this.positions);
        out.close();
      } 
      catch (FileNotFoundException e) {
        e.printStackTrace();
      } 
      catch (IOException e) {
        e.printStackTrace();
      }
    }
    
    // initialize the current items to the first sequence
    this.currentIndex = -1;
    this.current = getSequence(0);
    this.id = this.current.getId();
  }

  
  /**
   * Helper class that loads or retrieves the sequence at the given index.
   * @param index the index of the sequence to look up
   * @return the sequence object
   */
  private FastaSequence getSequence(int index) {
    if (this.sequences==null) {
      if (index!=this.currentIndex) {
        // load it
        this.current = new FastaSequence(this.files.get(index), this.aaSpec);
        this.currentIndex = index;
      }
      return current;
    }
    return this.sequences.get(index);
  }
  
  /**
   * Gets the array of individual protein sequences
   * @return the list of proteins sequences
   */
  public ArrayList<FastaSequence> getSequences() {
    return sequences;
  }
  
  /**
   * Helps translate the give position to a pair composed of a index of the
   * fasta sequence and the subindex in that sequence.  
   * @param position the absolute position
   * @return the relative position with the upper 32 bits as the sequence index
   *         and the lower 32 bits as the index in the sequence.
   */
  private long translate(long position) {
    int matchIndex = Collections.binarySearch(this.positions, position);
    
    long offset = 0;
    int sequenceIndex = 0;
    if (matchIndex < 0) {
      sequenceIndex = -matchIndex-1-1;
    }
    else {
      sequenceIndex = matchIndex-1;
    }
    if (sequenceIndex>=0) {
      offset = this.positions.get(sequenceIndex);
    }
    sequenceIndex++;
    return (((long)sequenceIndex)<<32) | ((int)(position-offset));
  }
  
  public int getAlphabetSize() {
    return current.getAlphabetSize();
  }

  public String getAnnotation(long position) {
    long pair = translate(position);
    //System.out.println(this.files.get((int)(pair>>>32)) + " ");
    return getSequence((int)(pair>>>32)).getAnnotation((int)pair);
  }

  public byte getByteAt(long position) {
    long pair = translate(position);
    return getSequence((int)(pair>>>32)).getByteAt((int)pair);
  }

  public int getId() {
    return this.id;
  }

  public String getMatchingEntry(long position) {
    long pair = translate(position);
    return getSequence((int)(pair>>>32)).getMatchingEntry((int)pair);
  }

  public String getMatchingEntry(String name) {
    for (FastaSequence sequence : this.sequences) {
      String match = sequence.getMatchingEntry(name);
      if (match!=null) return match;
    }
    return null;
  }
  
  public long getSize() {
    return this.positions.get(this.positions.size()-1);
  }

  public char toChar(byte b) {
    return current.toChar(b);
  }
  
  public String toString(byte[] sequence) {
    return current.toString(sequence);
  }

  public char getCharAt(long position) {
    long pair = translate(position);
    return getSequence((int)(pair>>>32)).getCharAt((int)pair);
  }
    
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    String directory = userHome+"/Data/Databases/Scerv/gen";
    FastaSequences pfs = new FastaSequences(directory, false);
    
    System.out.println("Total number of bases: " + pfs.getSize());
    for (long start = 0; start < pfs.getSize(); start++) {
      if (start%1000000==0) {
        if (pfs.isTerminator(start))
          System.out.println(pfs.getAnnotation(start));
      }
      pfs.getByteAt(start);
    }
  }

  public byte[] getBytes(int start, int end) {
    long pair1 = translate(start);
    long pair2 = translate(end);
    int seqIndex = (int)(pair1>>>32);
    return getSequence(seqIndex).getBytes((int)pair1, (int)pair2);
  }

  public boolean isInAlphabet(char c) {
    return current.isInAlphabet(c);
  }

  public boolean isTerminator(long position) {
    long pair = translate(position);
    return getSequence((int)(pair>>>32)).isTerminator((int)pair);
  }

  public boolean isValid(long position) {
    long pair = translate(position);
    return getSequence((int)(pair>>>32)).isValid((int)pair);
  }

  public byte toByte(char c) {
    return current.toByte(c);
  }

  public Collection<Character> getAlphabet() {
    return current.getAlphabet();
  }

  public Set<Byte> getAlphabetAsBytes() {
    return current.getAlphabetAsBytes();
  }

  public String getSubsequence(long start, long end) {
    long pair1 = translate(start);
    long pair2 = translate(end);
    int seqIndex = (int)(pair1>>>32);
    return getSequence(seqIndex).getSubsequence((int)pair1, (int)pair2);
  }

  public long getStartPosition(long position) {
    long pair = translate(position);
    long subStart = getSequence((int)(pair>>>32)).getStartPosition((int)pair);
    return position-(((int)pair)-subStart);
  }
  
}
