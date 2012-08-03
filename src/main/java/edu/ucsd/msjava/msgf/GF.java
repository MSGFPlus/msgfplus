package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.Annotation;
import edu.ucsd.msjava.msutil.Matter;

public interface GF<T extends Matter> {
	public boolean computeGeneratingFunction();
	public int getScore(Annotation annotation);
	public double getSpectralProbability(int score);
	public int getMaxScore();
	public ScoreDist getScoreDist();
}
