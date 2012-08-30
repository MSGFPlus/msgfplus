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
	
	private boolean considerBestMatchOnly = false;
	
	public MSGFPlusPSMSet(List<MSGFPlusMatch> msgfPlusPSMList, boolean isDecoy, CompactSuffixArray sa)
	{
		this.msgfPlusPSMList = msgfPlusPSMList;
		this.isDecoy = isDecoy;
		this.sa = sa;
	}
	
	public MSGFPlusPSMSet setConsiderBestMatchOnly(boolean considerBestMatchOnly)
	{
		this.considerBestMatchOnly = considerBestMatchOnly;
		return this;
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
			List<DatabaseMatch> dbMatchList;
			if(considerBestMatchOnly)
			{
				dbMatchList = new ArrayList<DatabaseMatch>();
				dbMatchList.add(match.getBestDBMatch());
			}
			else
				dbMatchList = match.getMatchList();
			
			for(DatabaseMatch m : dbMatchList)
			{
				String pepSeq = m.getPepSeq();
				
				boolean isDecoy = true; 
				for(int index : m.getIndices())
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

				float specEValue = (float)m.getSpecEValue();
				psmList.add(new ScoredString(pepSeq, specEValue));
				Float prevSpecEValue = peptideScoreTable.get(pepSeq);
				if(prevSpecEValue == null || specEValue < prevSpecEValue)
					peptideScoreTable.put(pepSeq, specEValue);
			}
		}
	}

}
