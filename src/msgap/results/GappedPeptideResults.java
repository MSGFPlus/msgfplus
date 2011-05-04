package msgap.results;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeMap;

/**
 * This class is to provide accessibility to the queries read from the result
 * file for a complete directory.
 * @author jung
 *
 */
public class GappedPeptideResults {

  /**
   * Holds the essential information of a given spectrum
   * @author jung
   *
   */
  private class SpecMetaInfo {
    private float pm;
    private int charge;
    private int scanNum;
    private String filename;
    private String ident;
    private String actMethod;
    
    private SpecMetaInfo(String filename, int scanNum, float pm, int charge, String actMethod, String ident) {
      this.filename = filename;
      this.pm = pm;
      this.charge = charge;
      this.scanNum = scanNum;
      this.ident = ident;
      
      if (actMethod!=null) this.actMethod = actMethod;
      else                 this.actMethod = "unknown_act_method";
      
    }
    
    @Override
    public String toString() {
      return String.format("%d\t%s\t%.3f\t%d\t%s", scanNum, filename, pm, charge, ident);
    }
    
  }
  
  // maps the spec id to the meta info of the spectrum
  private TreeMap<Integer,SpecMetaInfo> ids;
  // maps the spec id to the list of reconstructions by MSGap (indexes in the sequence file)
  private HashMap<Integer,ArrayList<Integer>> ids2sequences;
  // maps the sequence index to the spec id, multiple sequences might come from the same spectrum id
  private HashMap<Integer,Integer> sequences2ids;
  // the series of sequences as a list of ArrayList of integers
  private ArrayList<ArrayList<Integer>> sequences;
  // the largest specId stored by this object
  private int maxId;
  // the smallest specId stored by this object 
  private int minId;

  
  
  /**
   * Renders this object useless by clearing all the fields
   */
  public void clear() {
    ids.clear(); ids2sequences.clear(); sequences2ids.clear(); sequences.clear();  
  }
  
  /**
   * Default constructor.
   */
  public GappedPeptideResults() {
    this.ids = new TreeMap<Integer,SpecMetaInfo>();
    this.ids2sequences = new HashMap<Integer,ArrayList<Integer>>();
    this.sequences2ids = new HashMap<Integer,Integer>();
    this.sequences = new ArrayList<ArrayList<Integer>>();
    this.maxId = 0;
    this.minId = Integer.MAX_VALUE;
  }
  
  
  
  /**
   * Add a spectrum annotation 
   * @param specId the unique identifier of this item
   * @param sequence the sequence as a list of masses
   */
  public void addSequence(int specId, ArrayList<Integer> sequence) {
    
    if (!ids2sequences.containsKey(specId)) {
      ids2sequences.put(specId, new ArrayList<Integer>());  
    }
    
    ids2sequences.get(specId).add(this.sequences.size());
    sequences2ids.put(this.sequences.size(), specId);
    this.sequences.add(sequence);
  }
  
  
  
  /**
   * Adds the meta information of an spectrum and gets an unique identifier for
   * the given spectrum
   * @param specId the unique integer identifier of this item
   * @param scanNum the scan number of the spectrum
   * @param filename the filename of the spectrum
   * @param pm the parent mass
   * @param charge the charge
   * @param ident the annotation (identification)
   */
  public void addSpectrum(int specId, int scanNum, String filename, float pm, int charge, String actMethod, String ident) {
    SpecMetaInfo smi = new SpecMetaInfo(filename, scanNum, pm, charge, actMethod, ident);
    this.addSpectrum(specId, smi);
  }
  
  
  
  /**
   * Helper method to add the spectrum metainfo into the current object
   * @param specId the id of the item
   * @param smi the meta information
   */
  private void addSpectrum(int specId, SpecMetaInfo smi) {
    if (specId > this.maxId) this.maxId = specId;
    if (specId < this.minId) this.minId = specId;
    ids.put(specId, smi);
  }
  
  
  
  /**
   * Get the total number of sequences bundled by this object.
   * @return the total number of sequences.
   */
  public int getSequenceCount() {
    return this.sequences.size();
  }
  
  
  
  /**
   * Get the sequence at a specific index
   * @param index the index of the sequence
   * @return the sequence
   */
  public ArrayList<Integer> getSequenceAt(int index) {
    return this.sequences.get(index);
  }
  
  
  
  /**
   * Get the list of sequences stored by this object
   * @return the list of sequences
   */
  public ArrayList<ArrayList<Integer>> getSequences() {
    return this.sequences;
  }
  
  
  
  /**
   * The number of distinct spectra grouped in this object
   * @return the number of distinct spectra.
   */
  public int getSpectrumCount() {
    return this.ids.size();
  }
  
  public int getMaxId() {
    return this.maxId; 
  }
  
  public int getMinId() {
    return this.minId;
  }

  public Set<Integer> getIdList() {
    return this.ids.keySet();
  }
  
  /**
   * Create a new GPR object that does not include the id's specified from 
   * the parameters
   * @param exclude the id's not to include
   * @return the new GPR object
   */
  public GappedPeptideResults generateGPR(HashSet<Integer> exclude) {
    GappedPeptideResults ret = new GappedPeptideResults();
    for (int specId : this.ids2sequences.keySet()) {
      if (!exclude.contains(specId)) {
        ret.addSpectrum(specId, this.ids.get(specId));
        for (int seqId : this.ids2sequences.get(specId)) {
          ret.addSequence(specId, this.sequences.get(seqId));
        }
      }
    }
    return ret;
  }
  
  /**
   * Create a new GPR object that only contains the given id's
   * @param include include these id's
   * @return the new GPR object
   */
  public GappedPeptideResults retain(HashSet<Integer> include) {
    GappedPeptideResults ret = new GappedPeptideResults();
    for (int specId : this.ids2sequences.keySet()) {
      if (include.contains(specId)) {
        ret.addSpectrum(specId, this.ids.get(specId));
        for (int seqId : this.ids2sequences.get(specId)) {
          ret.addSequence(specId, this.sequences.get(seqId));
        }
      }
    }
    return ret;
  }
  
  
  /**
   * Given the information do a reverse lookup on the id of the spectrum identified
   * by the filename and scanNumebr given. This is inplemented as a linear search
   * of the metainformation stored.
   * @param filename the filename
   * @param scanNumber the scan number
   * @return the id of the spectrum or -1 if not found
   */
  public int getSpecId(String filename, int scanNumber) {
    for (int id : this.ids.keySet()) {
      SpecMetaInfo smi = this.ids.get(id);
      if (smi.filename.equals(filename) && scanNumber==smi.scanNum) return id;
    }
    return -1;
  }
  
  
  /**
   * Return the spectrum id of the given query index
   * @param sequenceIndex the index of the query in to inquire about
   * @return the scan number associated with this queryIndex
   */
  public int getSpecId(int sequenceIndex) { return sequences2ids.get(sequenceIndex); }
  
  /**
   * Return the scanNumber of the spectrum associated with given id
   * @param specId the id of the spectrum
   * @return the scan number associated with this queryIndex
   */
  public int getScanNumber(int specId) { return this.ids.get(specId).scanNum; }

  /**
   * Get the charge of the given item
   * @param specId the item to look up
   * @return the charge of the item
   */
  public int getCharge(int specId) { return this.ids.get(specId).charge; }
  
  /**
   * Get the spectrum precursor mass
   * @param specId the item to loop up
   * @return the precursor mass
   */
  public float getPrecursorMass(int specId) { return this.ids.get(specId).pm; }
  
  /**
   * Get the identification of the given spectrum
   * @param specId the item to look up
   * @return the identification/annotation if it exist
   */
  public String getIdent(int specId) { return this.ids.get(specId).ident; }
  
  /**
   * Get the activation method of the given spectrum
   * @param specId the item to look up
   * @return the activation method information
   */
  public String getActmethod(int specId) { return this.ids.get(specId).actMethod; }
  
  /**
   * Return the fileName of the spectrum associated with the given id
   * @param specId the item to look up
   * @return the filename of the given id
   */
  public String getFileName(int specId) {
    return this.ids.get(specId).filename;
  }
  
  /**
   * Check whether an specId is contained in this set of results.
   * @param specId the specId to check against
   * @return true if it is contained, false otherwise.
   */
  public boolean containsSpecId(int specId) {
    return this.ids.containsKey(specId);  
  }
  
  /**
   * Return the query specified by the index
   * @param index the index of the query to retrieve
   * @return the list of masses representing the query at the given index
   */
  public ArrayList<Integer> getQueryAt(int index) {
    return this.sequences.get(index);
  }
  
 
  
  /**
   * Get the spectrum id's contained in this object
   * @return the list of spectrum id's.
   */
  public Set<Integer> getSpecIds() {
    return this.ids.keySet();
  }

 
  /**
   * Write this object into a file  
   * @param out the out stream
   */
  public void toFile(PrintWriter out) {
    for (int specId : this.ids.keySet()) {
      // print the header
      out.printf("#%d\t%s\n", specId, this.ids.get(specId));
      for (int seqId : this.ids2sequences.get(specId)) {
        out.println(this.getSequenceAt(seqId));
      }
    }
  }
  
}
