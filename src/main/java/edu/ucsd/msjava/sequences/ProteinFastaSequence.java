package edu.ucsd.msjava.sequences;

import java.util.HashSet;

import edu.ucsd.msjava.msutil.AminoAcidSet;

/**
 * This class is a wrapper to the FastaSequence that uses Amino Acids as the
 * alphabet by default.
 * to amino acid masses.
 * @author jung
 *
 */
public class ProteinFastaSequence extends FastaSequence implements MassSequence {

  private AminoAcidSet alpha = Constants.AA;
  private byte[] masses;               // the translated masses 
  private HashSet<Long> invalids;      // positions that are invalid
  
  
/***** Helpers *****/
  private void initialize() {
    this.invalids = new HashSet<Long>();
    this.masses = new byte[(int)this.getSize()];
    for (long position=0; position < getSize(); position++) {
      if (isTerminator(position) || !alpha.contains(getCharAt(position))) {
        this.invalids.add(position);
        this.masses[(int)position] = (byte)0;
      }
      else {
        // we scale it back, so all amino acids fit in a byte from -127 to 128
        this.masses[(int)position] = (byte)(alpha.getAminoAcid(getCharAt(position)).getNominalMass()-100);
      }
    }
  }
  
  
/***** Constructors *****/  
   /**
    * Constructor using all (standard) letters in the fasta file as amino acids
    * @param filepath the path to the fasta file
    */
   public ProteinFastaSequence(String filepath) {
     super(filepath, edu.ucsd.msjava.sequences.Constants.AMINO_ACIDS_20, edu.ucsd.msjava.sequences.Constants.PROTEIN_FILE_EXTENSION);
     initialize();
   }
   
   /**
    * Constructor using a customized alphabet. See FastaSequence, for the syntax
    * of the alphabet argument.
    * @param filepath the path to the fasta file
    * @param alphabet the alphabet specification
    */
   public ProteinFastaSequence(String filepath, String alphabet) {
     super(filepath, alphabet, edu.ucsd.msjava.sequences.Constants.PROTEIN_FILE_EXTENSION);
     initialize();
   }
   
   /**
    * Constructor using a customized alphabet. See FastaSequence, for the syntax
    * of the alphabet argument.
    * @param filepath the path to the fasta file
    * @param alphabet the alphabet specification
    * @param aaSet the amino acid set to use
    */
   public ProteinFastaSequence(String filepath, String alphabet, AminoAcidSet aaSet) {
     super(filepath, alphabet, edu.ucsd.msjava.sequences.Constants.PROTEIN_FILE_EXTENSION);
     this.alpha = aaSet;
     initialize();
   }
   
   
/***** Member methods *****/
   public int getIntegerMass(long index) {
     return this.masses[(int)index]+100;
     /*
     AminoAcid aa = alpha.getAminoAcid(getCharAt(index));
     if (aa!=null) return aa.getNominalMass();
     return 0;*/
   }
   
   public int getIntegerMass(long start, long end) {
     int cumMass = 0;
     for (long i=start; i<end; i++) cumMass += getIntegerMass(i);
     return cumMass;
   }
   
   public boolean hasMass(long position) {
     return !invalids.contains(position) && position<this.getSize() && position>=0;
   }
   
   
/***** Main method to test the size of a sequence *****/
   public static void main(String[] args) {
     ProteinFastaSequence s = new ProteinFastaSequence(System.getProperty("user.home")+"/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants.fasta");
     //ProteinFastaSequence s = new ProteinFastaSequence(System.getProperty("user.home")+"/Data/Databases/Asp/pro/translated.fasta");
     System.out.println("Size of database: " + s.getSize());
   }
}