package msscorer;

import msutil.IonType;

public interface NewAdditiveScorer {
	// for scoring nodes
	public float getNodeScore(Partition part, IonType ionType, int rank);
	public float getMissingIonScore(Partition part, IonType ionType);
	
	// for scoring edges
	public float getErrorScore(Partition part, float error);
	// index => nn:0, ny:1, yn:2, yy:3
	public float getIonExistenceScore(Partition part, int index, float probPeak);
}