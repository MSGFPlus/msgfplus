package msscorer;

import java.util.List;

import msgf.NominalMass;
import msgf.ScoredSpectrum;
import msgf.ScoredSpectrumSum;

public class DBScanScorerSum extends ScoredSpectrumSum<NominalMass> implements SimpleDBSearchScorer<NominalMass> {

	private DBScanScorer[] scorerArr;
	public DBScanScorerSum(List<ScoredSpectrum<NominalMass>> scoredSpecList, int peptideMass) 
	{
		super(scoredSpecList);
		scorerArr = new DBScanScorer[scoredSpecList.size()];
		for(int i=0; i<scoredSpecList.size(); i++)
		{
			NewScoredSpectrum<NominalMass> scoredSpec = (NewScoredSpectrum<NominalMass>)scoredSpecList.get(i);
			scorerArr[i] = new DBScanScorer(scoredSpec, peptideMass);
		}
	}

	public int getScore(double[] prefixMassArr, int[] nominalPrefixMassArr,	int fromIndex, int toIndex, int numMods) {
		int sum = 0;
		for(DBScanScorer scorer : scorerArr)
			sum += scorer.getScore(prefixMassArr, nominalPrefixMassArr, fromIndex, toIndex, numMods);
		return sum;
	}
	
	@Override
	public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
		int sum = 0;
		for(DBScanScorer scoredSpec : scorerArr)
			sum += scoredSpec.getEdgeScore(curNode, prevNode, theoMass);
		return sum;
	}
	
}
