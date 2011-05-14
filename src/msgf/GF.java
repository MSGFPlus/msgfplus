package msgf;

import msutil.Annotation;
import msutil.Matter;

public interface GF<T extends Matter> {
	public boolean computeGeneratingFunction();
	public int getScore(Annotation annotation);
	public double getSpectralProbability(int score);
	public int getMaxScore();
	public ScoreDist getScoreDist();
}
