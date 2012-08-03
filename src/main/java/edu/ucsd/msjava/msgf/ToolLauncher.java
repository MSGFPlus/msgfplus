package edu.ucsd.msjava.msgf;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Iterator;

import edu.ucsd.msjava.msscorer.NewAdditiveScorer;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Spectrum;

public abstract class ToolLauncher {
	// Essential parameters, set by the constructor
	protected final Iterator<Spectrum> specIterator;
	protected final NewAdditiveScorer scorer;
	
	// Optional parameters set by builders. 
	
	protected float specProb = 1e-9f;
	
	protected boolean trypticOnly = true;

	// Tolerance
	protected Tolerance pmTolerance = new Tolerance(30, true);
	protected Tolerance fragTolerance = new Tolerance(30, true);

	protected float minParentMass = 400f;
	protected float maxParentMass = 2000f;
	protected int msgfScoreThreshold = 0;
	
	// Amino acid set, default: standard + Carboamidomethyl C
	protected AminoAcidSet aaSet;
	
	// output
	protected PrintStream out;
	
	/**
	 * A constructor specifies spectral file name and database file name. Database must be "fasta" format.
	 * @param specIterator spectra iterator.
	 * @param scorer a scorer object.
	 */
	protected ToolLauncher(Iterator<Spectrum> specIterator, NewAdditiveScorer scorer)
	{
		this.specIterator = specIterator;
		this.scorer = scorer;
		this.out = System.out;
		this.aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
	}
	
	/**
	 * A builder method to set spectral probability.
	 * @param specProb spectral probability
	 * @return this object.
	 */
	public ToolLauncher specProb(float specProb)	{ this.specProb = specProb; return this; }
	
	/**
	 * If this method is called, non-tryptic peptides are generated. 
	 * Otherwise, only peptides ends with 'K' or 'R' are generated. 
	 * @return this object.
	 */
	public ToolLauncher allowNonTryptic()	{ this.trypticOnly = false; return this; }
	
	
	/**
	 * Set parent mass tolerance.
	 * @param tolerance tolerance.
	 * @return this object.
	 */
	public ToolLauncher pmTolerance(Tolerance pmTolerance)
	{
		this.pmTolerance = pmTolerance;
		return this;
	}

	/**
	 * Set fragment mass tolerance.
	 * @param tolerance tolerance.
	 * @return this object.
	 */
	public ToolLauncher fragTolerance(Tolerance fragTolerance)
	{
		this.fragTolerance = fragTolerance;
		return this;
	}
	
	/**
	 * Set minimum parent mass.
	 * @param minParentMass minimum parent mass.
	 * @return this object.
	 */
	public ToolLauncher minParentMass(float minParentMass)
	{
		this.minParentMass = minParentMass;
		return this;
	}
	
	/**
	 * Set maximum parent mass.
	 * @param maxParentMass maximum parent mass.
	 * @return this object.
	 */
	public ToolLauncher maxParentMass(float maxParentMass)
	{
		this.maxParentMass = maxParentMass;
		return this;
	}

	/**
	 * Set max MSGF score threshold. Ignore all spectra whose best de novo scores are below thresholdScore. 
	 * @param thresholdScore max MS-GF score threshold.
	 * @return this object.
	 */
	public ToolLauncher msgfScoreThreshold(int thresholdScore)
	{
		this.msgfScoreThreshold = thresholdScore;
		return this;
	}
	
	/**
	 * Set the amino acid set.
	 * @param aaSet amino acid set.
	 * @return this object.
	 */
	public ToolLauncher aminoAcidSet(AminoAcidSet aaSet)
	{
		this.aaSet = aaSet;
		return this;
	}
	
	/**
	 * Set the output.
	 * @param outputFileName output file name.
	 * @return this object.
	 */
	public ToolLauncher outputFileName(String outputFileName)
	{
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return this;
	}
}
