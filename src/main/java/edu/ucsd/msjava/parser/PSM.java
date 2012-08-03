package edu.ucsd.msjava.parser;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeMap;

import edu.ucsd.msjava.msutil.Peptide;

/** 
 * A class represents Peptide-Spectrum Match (PSM)
 * @author sangtaekim
 */
public class PSM {
	public String getSpecFileName() {
		return specFileName;
	}
	public PSM specFileName(String fileName) {
		this.specFileName = fileName;
		return this;
	}
	public int getSpecIndex() {
		return specIndex;
	}
	public int getScanNum() {
		return scanNum;
	}
	public PSM scanNum(int scanNum) {
		this.scanNum = scanNum;
		return this;
	}
	public PSM specIndex(int specIndex) {
		this.specIndex = specIndex;
		return this;
	}
	public float getPrecursorMz() {
		return precursorMz;
	}
	public PSM precursorMz(float precursorMz) {
		this.precursorMz = precursorMz;
		return this;
	}
	public int getCharge() {
		return charge;
	}
	public PSM charge(int charge) {
		this.charge = charge;
		return this;
	}
	public String getTitle() {
		return title;
	}
	public PSM title(String title) {
		this.title = title;
		return this;
	}
	public Peptide getPeptide()		{ return peptide; }
	public String getPeptideStr()	{ return peptide.toString(); }
	public PSM peptide(Peptide peptide) {
		this.peptide = peptide;
		return this;
	}
	public String getPTM() {
		return ptm;
	}
	public PSM precedingResidue(char precedingResidue) {
		this.precedingResidue = precedingResidue;
		return this;
	}
	public char getPrecedingResidue() {
		return precedingResidue;
	}
	public PSM succeedingResidue(char succeedingResidue) {
		this.succeedingResidue = succeedingResidue;
		return this;
	}
	public char getSucceedingResidue() {
		return succeedingResidue;
	}
	public PSM ptm(String ptm) {
		this.ptm = ptm;
		return this;
	}
	public String getProtein() {
		return protein;
	}
	public PSM protein(String protein) {
		this.protein = protein;
		return this;
	}
	public float getProbScore() {
		return probScore;
	}
	public PSM probScore(float probScore) {
		this.probScore = probScore;
		return this;
	}
	public float getRawScore() {
		return rawScore;
	}
	public PSM rawScore(float rawScore) {
		this.rawScore = rawScore;
		return this;
	}
	public float getScore(String scoreName) {
		Float score = optionalScores.get(scoreName);
		if(score == null)
			return Float.NaN;
		else
			return score;
	}
	public ArrayList<String> getScoreNames() {
		return new ArrayList<String>(optionalScores.keySet());
	}
	
	public PSM score(String name, float score) {
		if(optionalScores == null)
			optionalScores = new TreeMap<String,Float>();
		optionalScores.put(name, score);
		
		return this;
	}
	
	public static class PSMSpecFileAndScanNumComparator implements Comparator<PSM> {
		public int compare(PSM o1, PSM o2) {
			int fileNameComparison = o1.specFileName.compareTo(o2.specFileName);
			if(fileNameComparison != 0)
				return fileNameComparison;
			else
			{
				if(o1.scanNum > o2.scanNum)
					return 1;
				else if(o1.scanNum == o2.scanNum)
					return 0;
				else
					return -1;
			}
		}
	}
	
	static class PSMSpecNumComparator implements Comparator<PSM> {
		public int compare(PSM o1, PSM o2) {
			if(o1.specIndex > o2.specIndex)
				return 1;
			else if(o1.specIndex == o2.specIndex)
				return 0;
			else
				return -1;
		}
	}

	static class PSMScanNumComparator implements Comparator<PSM> {
		public int compare(PSM o1, PSM o2) {
			if(o1.scanNum > o2.scanNum)
				return 1;
			else if(o1.scanNum == o2.scanNum)
				return 0;
			else
				return -1;
		}
	}
	
	static class PSMProbScoreComparator implements Comparator<PSM> {
		public int compare(PSM o1, PSM o2) {
			if(o1.probScore > o2.probScore)
				return 1;
			else if(o1.probScore == o2.probScore)
				return 0;
			else
				return -1;
		}
	}
	
	// Spectrum File name
	private String specFileName = null;
	
	// Spectrum
	private int specIndex = -1;	// sequencial number of the spectrum in the file, starting from 0
	private int scanNum = -1;	// scan number
	private String title = null;	// Title
	private float precursorMz = 0;
	private int charge = 0;
	
	// Peptide
	private String protein = null;	// protein annotation
	private char precedingResidue = '*';
	private char succeedingResidue = '*'; 
	private String ptm = null;	// TODO: need to be refined later
	private Peptide peptide;
	
	// Score
	private float probScore;	// E-value, P-value or FDR, the smaller the better
	private float rawScore;		// Raw score, the greater the better
	
	// optional score
	private TreeMap<String, Float> optionalScores = null;
	
	public String toString()
	{
		return specIndex+"\t"+scanNum+"\t"+peptide+"\t"+charge+"\t"+probScore;
	}
	
	public String toStringAllFields() 
	{
		return (specIndex > 0 ? specIndex : "") + "\t"+ (scanNum > 0 ? scanNum : "") + "\t" +
		(title != null ? title : "") + "\t" + (charge > 0 ? charge : "") + "\t" +
		(peptide != null ? peptide : "") + "\t" + (protein != null ? protein : "") + "\t" +
		(ptm != null ? ptm : "") + "\t" + probScore + "\t" + rawScore;
	}
	
}
