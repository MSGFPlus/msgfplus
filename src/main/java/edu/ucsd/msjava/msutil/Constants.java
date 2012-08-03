package edu.ucsd.msjava.msutil;



import java.text.DecimalFormat;


public class Constants {

  public static final float EPSILON = 1E-6f;     // very small number for float comparisons 
  
  public static final float MILLION = 1000000.0f;
  
  public static final float INTEGER_MASS_SCALER = 0.999497f;
  public static final float INTEGER_MASS_SCALER_HIGH_PRECISION = 274.335215f;
  
  public static final float ANALYSIS_VERSION = 1.0f;

  public static boolean   COMPARE_WITH_MASCOT = false;

  public static boolean   PRINT_PEAK_ERROR = false;

  public static boolean     PARAMETER_OPTIMIZER = false;

  public static boolean   COMPARE_WITH_INSPECT = false;

  public static boolean   RANDOM_SPEC_SELECT = false;

  public static int     RANDOM_SPEC_SPELECT_SIZE = 1000;

  

  public static final float UNIT_MASS = 1.f;

  public static final float B_ION_OFFSET = UNIT_MASS;

  public static final float Y_ION_OFFSET = UNIT_MASS*19;

  public static float   offsetMinPerGap = -100.f;

  public static float   offsetMaxPerGap = 500.f;

  public static float   offsetMaxPerPeptide = 500.f;

  public static float   offsetMinPerPeptide = -150.f;

  public static float   massTolerance = 0.5f;

  public static float   precursorTolerance = 1.5f;

  public static float   selectionWindowSize = 70;

  public static int     minNumOfPeaksInWindow = 2;

  public static int     maxNumOfPeaksInWindow = 100;  // currently not defined

  public static float   minPeptideMass = 400.f;

  public static float   maxPeptideMass = 4000.f;

  public static int     minTagLength = 2;

  public static int     minTagLengthPeptideShouldContain = 3;

  public static float   tagChainPruningRate = 0.5f;



  public static String    IDENTIFIER = "Ewha_HSP27";

  public static int     MiscleavageForProteinID = 1;

  public static int     MiscleavageForPTMSearch = 5;

  public static String    PROTEIN_DB_NAME = "hsp27.fasta";

  public static String    SPECTRUM_FILE_NAME = "";

//  public static SpectraFileType SPECTRA_FILE_TYPE = SpectraFileType.DTA;

  public static String    INSTRUMENTS_NAME = "QTOF";

  public static String    PTM_FILE_NAME = "PTMDB.xml";

  

  public static final int MAX_TAG_SIZE = 400;

  public static final int MAX_PEPTIDE_LENGTH = 50;

  // should add XML form

  public static float minNormIntensity = 0.1f;

  

  // for Peptide DB

  public static final int proteinIDModeSeqLength  = 3;

  public static final String SOURCE_PROTEIN_FILE_NAME = "sourceProtein.mprot";  

  // for PTM DB

  public static final int maxPTMSearchLength    = 12;

  public static final int maxPTMSizePerGap    = 5;

  // public static final int maxPTMOccurrence   = 5;



  public static final String  SPECTRUM_EXTENSION = ".unidta";

  public static final String  ANALYSIS_EXTENSION = ".unidrawing";

  public static final int   ThresholdForCompression = 1000000000;

  

  public static final String  UNIMOD_FILE_NAME = "unimod.xml";

  

  // for mother mass correction for LTQ/LCQ

  public static final float MINIMUM_PRECURSOR_MASS_ERROR = -1.5f;

  public static final float MAXIMIM_PRECURSOR_MASS_ERROR = 1.5f;

  // public static final boolean  IS_TRYPTIC = true; // not used currently



  // if true, write unidrawing only tag chains whose all gaps are annotated

  public static final boolean writeAnnotatedTagChainOnly = false;

  

  public static final int   MINIMUM_SHARED_PEAK_COUNT = 2; 

  

  // for offset

  public static final int newLineCharSize = new String("\r\n").getBytes().length;

  

  public static int getMaxPTMOccurrence( int seqLength )

  {

    if (seqLength>6) return 1;

    else if (seqLength>4) return 2;

    else return seqLength;

  }

  

  public static boolean equal(float v1, float v2)

  {

    if(Math.abs(v1-v2) < massTolerance)

      return true;

    else

      return false;

  }



  public static boolean equal(float v1, float v2, float tolerance)

  {

    if(Math.abs(v1-v2) <= tolerance)

      return true;

    else

      return false;

  }

  

  

  public static String  getString(float value)

  {

    return new DecimalFormat("#.###").format(value).toString();

  }

  

  public static float MASS_CAL_STD_THRESHOLD = 0.1f;

  public static float PTM_ADD_PENALTY = 0.2f;

  

  // Can use Wysocki paper results?
  /*
  public static float getMissingPenaltyWeight(PeakProperty property)

  {

    if(property == PeakProperty.Y_ION)

      return 0.4f;  

    else if(property == PeakProperty.Y_MINUS_NH3_ION)

      return 0;

    else if(property == PeakProperty.Y_MINUS_H2O_ION)

      return 0;

    else if(property == PeakProperty.B_ION)

      return 0;

    else if(property == PeakProperty.A_ION)

      return 0;

    else

      return 0;

  }
  */

  public static float getNotExplainedPenaltyWeight()

  {

    return 0.15f;

  }

}

