package edu.ucsd.msjava.suffixarray;

import java.io.*;
import java.nio.*;
import java.nio.channels.FileChannel;
import java.util.Arrays;
import java.util.TreeSet;


import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.Peptide;


/**
 * The MassArray class is for fast localization of a mass given a database of
 * amino acids.
 * @author jung
 *
 */
public class MassArray {

  // constants
  private static final String MASS_FILE_EXTENSION = ".marray";
  private static final float MAX_GAP_MASS = 1000.0f;
  
  // the size of an int primitive type in bytes
  private static final int INT_BYTE_SIZE = Integer.SIZE / Byte.SIZE;
  private static final int FLOAT_BYTE_SIZE = Float.SIZE / Byte.SIZE;
  
  /***** MEMBERS HERE *****/  
  private int size;
  private float[] masses;    // sorted array by mass
  private IntBuffer starts;
  private IntBuffer ends;
  
  
  /**
   * Constructor that creates a MassArray from the amino acid sequence in fasta
   * format and the name of the output file. This class uses a FastaSequence
   * object because the amino acid nomenclature and mass is needed to calculate
   * the mass.
   * @param sequence the FastaSequence object.
   * @param massArrayFile the path of the output file.
   */
  public MassArray(SuffixArraySequence sequence, String massArrayFile) {
    
    // if the does not exists create it 
    if(!new File(massArrayFile).exists())        
      createMassArrayFile(sequence, massArrayFile);
    
    // read the file in
    int id = readMassArrayFile(massArrayFile);
    
    // check that the files are consistent
    if(id != sequence.getId()) {
      System.err.println(massArrayFile + " was not created from the sequence " + sequence.getBaseFilepath());
      System.err.println("Please recreate the suffix array file.");
      System.exit(-1);
    }
  }
  
  
  /**
   * Constructor that takes the FastaSequence object as the only input. It uses
   * the default file name to check or create the pre-processed MassArray object.
   * @param sequence the FastaSequence object.
   */
  public MassArray(SuffixArraySequence sequence) {
    this(sequence, sequence.getBaseFilepath()+MASS_FILE_EXTENSION);
  }
  
  
  /**
   * This is the general matching routine that returns all matching elements 
   * in the MassArray. It can also return a null MatchSet if nothing is found.
   * @param lower the lower mass limit to match.
   * @param upper the upper mass limit to match.
   * @return a MatchSet representing the matches.
   */
  public MatchSet findAll(float lower, float upper) {
    // System.out.println("Size of masses: " + size);
    // System.out.println("First items: " + masses[0] + ", " + masses[1] + ", " + masses[size-1]);
    int matchIndex = Arrays.binarySearch(masses, lower);
    
    if(matchIndex < 0) {
      // no exact match
      matchIndex = -(matchIndex + 1);
    }
    
    // rewind the matchIndex until we can find the lowest insertion point
    for(int i = Math.min(matchIndex, this.size-1); i >= 0; i--) {
      if(this.masses[i] < lower) {
        matchIndex = i + 1;
        break;
      }
    }
    
    MatchSet ms = new MatchSet();
    
    // add the matches
    for(int i = matchIndex; i < this.size; i++) {
      if(masses[i] > upper) {
        break;
      }
      ms.add(this.starts.get(i), this.ends.get(i));
    }
    return ms;
  }
  
  
  
  /**
   * Debugging subroutine.
   */
  private static void debug() {
    long time = System.currentTimeMillis();
    String fastaFile;
    String userHome = System.getProperty("user.home");
    
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    fastaFile = userHome+"/Data/Databases/small.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    
    System.out.println("File name: "+fastaFile);
    SuffixArraySequence fs = new SuffixArraySequence(fastaFile);
    System.out.println("Total number of characters: " + fs.getSize());
    System.out.println("Alphabet size: " + fs.getAlphabetSize());
    System.out.println("Time to complete: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    
    // create the mass array and suffix array
    MassArray ma = new MassArray(fs);
    
    //SuffixArray sa = new SuffixArray(fs);
    // MatchSet exactMathSet = sa.findAll("M");
    // System.out.println(exactMathSet.getSize());
    
    MatchSet ms = ma.findAll(500f, 500.1f);
    
    for(int i = 0; i < ms.getSize(); i++) {
      int start = ms.getStart(i);
      int end = ms.getEnd(i);
      String match = fs.getSubsequence(start, end);
      System.out.println("Start: " + start + ". End: " + end + ". Sequence: " + match + " .Mass: " + Peptide.getMassFromString(match));
    }
    
  }
  
  
  /**
   * @param args
   */
  public static void main(String[] args) {
    debug();
  }
  
  
  private static class MassObject implements Comparable<MassObject>{
    float mass;
    int start, end;
 
    /**
     * Defaul constructor.
     * @param mass mass of this fragment.
     * @param start start position.
     * @param end end position.
     */
    public MassObject(float mass, int start, int end) {
      this.mass = mass;
      this.start = start;
      this.end = end;
    }
    
    
    public int compareTo(MassObject other) {
      if(this.mass > other.mass)       return 1;
      if(this.mass < other.mass)       return -1;
      
      if(this.start > other.start)     return 1;
      if(this.start < other.start)     return -1;
      
      if(this.end > other.end)       return 1;
      if(this.end < other.end)       return -1;       
      return 0;
    }  
  }
  
  
  /**
   * Reads in the file assuming the same format that was written.
   * @param massArrayFile the path to the file.
   * @return returns the id of this file for consistency check.
   */
  private int readMassArrayFile(String massArrayFile) {
    try {
      // read the size and id
      DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(massArrayFile)));
      this.size = in.readInt();
      int id = in.readInt();
      
      // use the FileChannel to read the chunks of data
      FileChannel fc = new FileInputStream(massArrayFile).getChannel();
      
      int startPos = 2*FLOAT_BYTE_SIZE;
      int sizeOfMasses = this.size*FLOAT_BYTE_SIZE;
      FloatBuffer fb = fc.map(FileChannel.MapMode.READ_ONLY, startPos, sizeOfMasses).asFloatBuffer();
      if(fb.hasArray()) {
        // has a backing array
        this.masses = fb.array();
      }
      else {
        // create the array
        //System.out.println("Creating the float array");
        this.masses = new float[this.size];
        for(int i = 0; i < this.size; i++) {
          this.masses[i] = fb.get(i);
        }
      }
      
      // get the other buffers
      startPos += sizeOfMasses;
      int startsSize = this.size*INT_BYTE_SIZE;
      this.starts = fc.map(FileChannel.MapMode.READ_ONLY, startPos, startsSize).asIntBuffer();
      
      startPos += startsSize;
      int endsSize = this.size*INT_BYTE_SIZE;
      this.ends = fc.map(FileChannel.MapMode.READ_ONLY, startPos, endsSize).asIntBuffer();
      fc.close();
      
      return id;
    }
    catch(IOException e) {
      e.printStackTrace();
      System.exit(-1);
    }
    
    return 0;
  }
  
  
  /**
   * Helper method that creates the massArrayFile.This routine requires a FastaSequence
   * object because we need to deduct the mass of the amino acids according to
   * the amino acid letters. The current format writes all masses, all starts and
   * all ends.Of course also write the size and id as first and second integers.
   * @param sequence the FastaSequence object that represents the database (text).
   * @param massFile the output file.
   */
  private void createMassArrayFile(SuffixArraySequence sequence, String massArrayFile) {
    System.out.println("Creating the mass array sorted file of size: " + sequence.getSize() + " characters");
    
    TreeSet<MassObject> sortedMasses = new TreeSet<MassObject>();
    
    int sequenceLen = (int)sequence.getSize(); 
    for(int index = 0; index < sequenceLen; index++) {
      float cumMass = 0f;
      int newIndex = index;
      while(cumMass < MAX_GAP_MASS && newIndex < sequenceLen) {
        AminoAcid aa = AminoAcid.getStandardAminoAcid(sequence.getCharAt(newIndex));
        if(aa != null) {
          // additional amino acid, record this.
          cumMass += aa.getMass();
          newIndex++;
          // note the the end is open
          sortedMasses.add(new MassArray.MassObject(cumMass, index, newIndex));
        }
        else {
          break;
        }
      }
    }
    
    System.out.println("Length of mass array " + sortedMasses.size());
    
    // write out the numbers
    try {
      DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(massArrayFile)));
      out.writeInt(sortedMasses.size());
      out.writeInt(sequence.getId());
      
      // write the masses
      for(MassObject mo : sortedMasses) {
        out.writeFloat(mo.mass);
      }
      
      // write in start indexes
      for(MassObject mo : sortedMasses) {
        out.writeInt(mo.start);
      }
      
      // write in start indexes
      for(MassObject mo : sortedMasses) {
        out.writeInt(mo.end);
      }
    
      out.close();
    }
    catch(IOException e) {
      e.printStackTrace();
      System.exit(-1);
    }
    
    
  }

}
