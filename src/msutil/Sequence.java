package msutil;

import java.util.*;
import msgf.MassListComparator;
import msgf.Tolerance;


/**
 * Superclass for a list of masses. Peptide, GappedPeptide, Tag should extend
 * this class.
 * @author jung
 *
 */
public class Sequence<T extends Matter> extends ArrayList<T>  {

  //this is recommended for Serializable objects
  static final private long serialVersionUID = 1L;
  
  
  /**
   * Sums up the masses of this Sequence.
   * @return the mass in Daltons of the mono isotopic masses.
   */
  public float getMass() {
    return getMass(0, this.size());
  }
  
  /**
   * Sums up the masses of this Sequence (double-precision).
   * @return the mass in Daltons of the mono isotopic masses (double-precision).
   */
  public double getAccurateMass() {
    return getMass(0, this.size());
  }
  
  /**
   * Sums up the masses of the specified range of masses (half open
   * intervals).
   * @param from the index of the starting mass (inclusive).
   * @param to the end index of the mass (exclusive).
   * @return the mass in Daltons.
   */
  public float getMass(int from, int to) {
    from = java.lang.Math.max(from, 0);
    to = java.lang.Math.min(to, this.size());
    float sum = 0.f;
    for(int i=from; i<to; i++)
      sum += this.get(i).getMass();
    return sum;
  }

  /**
   * Similar to getMass(), but returns double.
   * @param from the index of the starting mass (inclusive).
   * @param to the end index of the mass (exclusive).
   * @return the mass in Daltons (double).
   */
  public double getAccurateMass(int from, int to) {
    from = java.lang.Math.max(from, 0);
    to = java.lang.Math.min(to, this.size());
    double sum = 0;
    for(int i=from; i<to; i++)
      sum += this.get(i).getAccurateMass();
    return sum;
  }
  
  /**
   * Returns a subsequence of the specified range (half open intervals).
   * @param fromIndex the index of the starting subsequence (inclusive).
   * @param toIndex the end index of the subsequence (exclusive)
   * @return a subsequence of specified range
   */
  public Sequence<T> subSequence(int fromIndex, int toIndex) {
	  return (Sequence<T>)super.subList(fromIndex, toIndex);
  }
  
  /**
   * String representation of this sequence.
   * @return the String representing the amino acid letters in this sequence.
   */
  public String toString() {
    StringBuffer output = new StringBuffer();
    for(T matter : this) {
      output.append(matter.toString()+" "); 
    }
    return output.toString();
  }
  
  /**
   * Returns the union of two input sequences.
   * @param seq1 the first sequence
   * @param seq2 the second sequence
   * @return the union of seq1 and seq2
   */
  public static <T extends Matter> Sequence<T> getIntersection(Sequence<T> seq1, Sequence<T> seq2)
  {
	  Sequence<T> union = new Sequence<T>();
	  HashSet<T> set = new HashSet<T>();
	  for(T m : seq1)
		  set.add(m);
	  for(T m : seq2)
		  if(set.contains(m))
			  union.add(m);
	  return union;
  }
  
  /**
   * Checks if this sequence matches to the specified peptide within the input tolerance
   * @param peptide	Peptide.
   * @param tolerance Tolerance.
   * @return	True if matches, false otherwise.
   */
  public boolean isMatchedTo(Peptide peptide, Tolerance tolerance, boolean isPrefix)
  {
	  ArrayList<Mass> pepMassList = new ArrayList<Mass>();
	  float mass = 0;
	  for(int i=0; i<peptide.size(); i++)
	  {
		  if(isPrefix)
			  mass += peptide.get(i).getMass();
		  else
			  mass += peptide.get(peptide.size()-1-i).getMass();
		  pepMassList.add(new Mass(mass));
	  }
	  ArrayList<Mass> massList = new ArrayList<Mass>();
	  for(int i=0; i<this.size(); i++)
		  massList.add(new Mass(this.get(i).getMass()));
	  MassListComparator<Mass> comparator = new MassListComparator<Mass>(pepMassList, massList);
	  int matchSize = comparator.getMatchedList(tolerance).length;
	  return (matchSize == this.size());
  }
  
  /**
   * Checks if this sequence matches to the specified peptide. Use nominal masses.
   * @param peptide	Peptide.
   * @param isTolerancePPM	Tolerance is interpreted as ppm is true. If false, Tolerance is interpreted as Da.
   * @return	True if matches, false otherwise.
   */
  public boolean isMatchedToNominalMasses(Peptide peptide, boolean isPrefix)
  {
	  HashSet<Integer> massList = new HashSet<Integer>();
	  int mass = 0;
	  for(int i=0; i<peptide.size(); i++)
	  {
		  if(isPrefix)
			  mass += peptide.get(i).getNominalMass();
		  else
			  mass += peptide.get(peptide.size()-1-i).getNominalMass();
		 massList.add(mass);
	  }
	  for(Matter m : this)
	  {
		  if(!massList.contains(m.getNominalMass()))
			  return false;
	  }
	  return true;
  }
  
  /**
   * Converts this sequence into an array of masses.
   * @return a mass array of this object. null if 
   */
  public float[] toMassArray()
  {
	  float[] massArr = new float[this.size()];
	  int index = 0;
	  for(T m : this)
		  massArr[index++] = m.getMass();
	  return massArr;
  }
  
  /*
  public static Sequence<NominalMass> getIntMassGappedPeptide(ArrayList<Peptide> dictionary, float minProbability, boolean prefix)
  {
	  GappedPeptide<IntMass> gp = new GappedPeptide<IntMass>();
	  
	  Hashtable<IntMass, Integer> hist = new Hashtable<IntMass, Integer>();

	  for(Peptide peptide : dictionary)
	  {
		  IntMass mass = new IntMass(0);
		  for(int i=0; i<peptide.size(); i++)
		  {
			  AminoAcid aa;
			  if(prefix)
				  aa = peptide.get(i);
			  else
				  aa = peptide.get(peptide.size()-1-i);
			  mass = mass.getAddition(IntMass.getIntMass(aa.getMass()));
			  Integer occ = hist.get(mass);
			  if(occ == null)
				  hist.put(mass, 1);
			  else
				  hist.put(mass, occ+1);
		  }
	  }

	  Iterator<Entry<IntMass, Integer>> itr = hist.entrySet().iterator();
	  while(itr.hasNext())
	  {
		  Entry<IntMass, Integer> entry = itr.next();
		  if(entry.getValue() > dictionary.size()*minProbability)
			  gp.add(entry.getKey());
	  }
	  Collections.sort(gp);
	  return gp;	  
  }
   */
  /*
   * Returns the gapped peptide as an array of compositions/ 
   * The spectral profile of dictionary is generated and prm compositions exceeding minProbability are selected.
   * @param dictionary	the spectral dictionary
   * @param minProbability	the threshold of profile probability
   * @param prefix	gapped peptide is a set of prefixes if true, suffixes if false
   * @returns the gapped peptide of compositions 
  public static GappedPeptide<Composition> getCompositionGappedPeptide(ArrayList<Peptide> dictionary, float minProbability, boolean prefix)
  {
	  GappedPeptide<Composition> gp = new GappedPeptide<Composition>();
	  
	  Hashtable<Composition, Integer> hist = new Hashtable<Composition, Integer>();

	  for(Peptide peptide : dictionary)
	  {
		  Composition composition = new Composition(0);
		  for(int i=0; i<peptide.size(); i++)
		  {
			  AminoAcid aa;
			  if(prefix)
				  aa = peptide.get(i);
			  else
				  aa = peptide.get(peptide.size()-1-i);
			  composition = composition.getAddition(aa.getComposition());
			  Integer occ = hist.get(composition);
			  if(occ == null)
				  hist.put(composition, 1);
			  else
				  hist.put(composition, occ+1);
		  }
	  }

	  Iterator<Entry<Composition, Integer>> itr = hist.entrySet().iterator();
	  while(itr.hasNext())
	  {
		  Entry<Composition, Integer> entry = itr.next();
		  if(entry.getValue() > dictionary.size()*minProbability)
			  gp.add(entry.getKey());
	  }
	  Collections.sort(gp);
	  return gp;
  }  
  */
}
