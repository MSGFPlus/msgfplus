package edu.ucsd.msjava.sequences;

import edu.ucsd.msjava.msutil.AminoAcidSet;


/**
 * This class contains the hardcode or preset values for the sequence classes
 * @author jung
 *
 */
public class Constants {
  
  /**
   * This string contains the 20 standard amino acids. 
   */
  public static final String AMINO_ACIDS_20 = "A:C:D:E:F:G:H:I:K:L:M:N:P:Q:R:S:T:V:W:Y";
  
  /**
   * <p>This string contains the 19 standard amino acids where the L is replaced by I.
   * The syntax for alphabet encoding is very simple. All amino acids that are
   * grouped together by the token separator ":" are considered equivalent when
   * doing the mapping. When doing the reverse mapping the first letter in the 
   * group is treated as the representative of the group.</p>
   * <p>For example, the contents of this String are "A:C:D:E:F:G:H:IL:K:M:N:P:Q:R:S:T:V:W:Y".</p>
   */
  public static final String AMINO_ACIDS_19 = "A:C:D:E:F:G:H:IL:K:M:N:P:Q:R:S:T:V:W:Y";
  
  /**
   * This string contains the 18 standard amino acids where the L is replaced by I and the Q by K.
   */
  public static final String AMINO_ACIDS_18 = "A:C:D:E:F:G:H:IL:KQ:M:N:P:R:S:T:V:W:Y";
  
  /**
   * This string contains the 18 standard amino acids where the L is replaced by I and the Q by K.
   */
  public static final String AMINO_ACIDS_18_X = "A:C:D:E:F:G:H:IL:KQ:M:N:P:R:S:T:V:W:X:Y";
  
  /**
   * The extension of the permanent storage files for the regular FastaSequence objects
   */
  public static final String FILE_EXTENSION = ".seq";
  
  /**
   * The extension of the permanent storage files for the ProteinFastaSequence objects 
   */
  public static final String PROTEIN_FILE_EXTENSION = ".pseq";
  
  /**
   * Add this suffix to the file extension for the annotation files
   */
  public static final String ANNO_FILE_SUFFIX = "anno";
  
  /**
   * The terminator byte representation.
   */
  public static final byte TERMINATOR = 0;
  
  /**
   * The terminator byte representation.
   */
  public static final char TERMINATOR_CHAR = '_';
  
  /**
   * The byte representation of the invalid character.
   */
  public static final byte INVALID_CHAR_CODE = 1;
  
  /**
   * The character representation of the invalid character.
   */
  public static final char INVALID_CHAR = '?';
  
  /**
   * Minumum number of peaks per spectrum.
   */
  public static final int MIN_NUM_PEAKS_PER_SPECTRUM = 20;
 
  public static final AminoAcidSet AA = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCysWithTerm();
}
