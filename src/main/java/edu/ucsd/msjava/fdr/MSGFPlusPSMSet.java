package edu.ucsd.msjava.fdr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import edu.ucsd.msjava.msdbsearch.CompactSuffixArray;
import edu.ucsd.msjava.msdbsearch.DatabaseMatch;
import edu.ucsd.msjava.msdbsearch.MSGFPlusMatch;
import edu.ucsd.msjava.ui.MSGFPlus;

public class MSGFPlusPSMSet extends PSMSet {

	private final List<MSGFPlusMatch> msgfPlusPSMList;
	private final boolean isDecoy;
	private final CompactSuffixArray sa;
	
	public MSGFPlusPSMSet(List<MSGFPlusMatch> msgfPlusPSMList, boolean isDecoy, CompactSuffixArray sa)
	{
		this.msgfPlusPSMList = msgfPlusPSMList;
		this.isDecoy = isDecoy;
		this.sa = sa;
	}
	
	@Override
	public boolean isGreaterBetter() 
	{
		return false;
	}

	// set-up 	ArrayList<ScoredString> psmList and HashMap<String,Float> peptideScoreTable
	@Override
	public void read() 
	{
		psmList = new ArrayList<ScoredString>();		
		peptideScoreTable = new HashMap<String,Float>();
		
		for(MSGFPlusMatch match : msgfPlusPSMList)
		{
			DatabaseMatch bestMatch = match.getBestDBMatch();
			String pepSeq = bestMatch.getPepSeq();
			
			boolean isDecoy = true; 
			for(int index : bestMatch.getIndices())
			{
				String protAcc = sa.getSequence().getAnnotation(index);
				if(!protAcc.startsWith(MSGFPlus.DECOY_PROTEIN_PREFIX))
				{
					isDecoy = false;
					break;
				}
			}
			
			if(this.isDecoy != isDecoy)
				continue;

			float specEValue = (float)bestMatch.getSpecEValue();
			psmList.add(new ScoredString(pepSeq, specEValue));
			Float prevSpecEValue = peptideScoreTable.get(pepSeq);
			if(prevSpecEValue == null || specEValue < prevSpecEValue)
				peptideScoreTable.put(pepSeq, specEValue);
		}
	}

}
