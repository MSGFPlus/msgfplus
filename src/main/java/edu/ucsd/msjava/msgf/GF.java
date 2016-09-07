package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.Annotation;
import edu.ucsd.msjava.msutil.Matter;

public interface GF<T extends Matter> {
    boolean computeGeneratingFunction();

    int getScore(Annotation annotation);

    double getSpectralProbability(int score);

    int getMaxScore();

    ScoreDist getScoreDist();
}
