package suffixtree.matches;

import java.util.ArrayList;

import msutil.Peptide;

import sequences.MassSequence;

public class MegaPrefixSuffixMatchObject extends MatchObject {

  private ArrayList<Coor> prefixes, suffixes;
  
  public MegaPrefixSuffixMatchObject(ArrayList<Coor> prefixes, ArrayList<Coor> suffixes, MassSequence db, ArrayList<Integer> query, int queryIndex) {
    this.prefixes = prefixes;
    this.suffixes = suffixes;
    setQuery(query);
    setQueryIndex(queryIndex);
  }

  @Override
  public Peptide getPeptide() {
    // generate 1 peptide
    return new Peptide(getMatchAsString());
  }

  @Override
  public Peptide getUnmodifiedPeptide() {
    return getPeptide();
  }

  @Override
  public int getStart() {
    return (int)prefixes.get(0).getStart();
  }

  @Override
  public int getEnd() {
    return (int)suffixes.get(0).getEnd();
  }


  @Override
  public String getSummaryLine(String filename, int scanNum, String actMethod,
      float pm, int charge, float offset) {
    // TODO Auto-generated method stub
    return null;
  }
  
  
}
