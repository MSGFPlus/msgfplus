package edu.ucsd.msjava.fdr;

import edu.ucsd.msjava.msdbsearch.CompactSuffixArray;
import edu.ucsd.msjava.msdbsearch.DatabaseMatch;
import edu.ucsd.msjava.msdbsearch.MSGFPlusMatch;
import edu.ucsd.msjava.ui.MSGFPlus;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class MSGFPlusPSMSet extends PSMSet {

    private final List<MSGFPlusMatch> msgfPlusPSMList;
    private final boolean isDecoy;
    private final CompactSuffixArray sa;
    private final String decoyProteinPrefix;

    private boolean considerBestMatchOnly = false;

    public MSGFPlusPSMSet(
            List<MSGFPlusMatch> msgfPlusPSMList,
            boolean isDecoy,
            CompactSuffixArray sa,
            String decoyProteinPrefix) {

        this.msgfPlusPSMList = msgfPlusPSMList;
        this.isDecoy = isDecoy;
        this.sa = sa;

        if (decoyProteinPrefix == null || decoyProteinPrefix.trim().isEmpty())
            this.decoyProteinPrefix = MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX;
        else
            this.decoyProteinPrefix = decoyProteinPrefix;
    }

    public MSGFPlusPSMSet setConsiderBestMatchOnly(boolean considerBestMatchOnly) {
        this.considerBestMatchOnly = considerBestMatchOnly;
        return this;
    }

    @Override
    public boolean isGreaterBetter() {
        return false;
    }

    // set-up ArrayList<ScoredString> psmList and HashMap<String,Float> peptideScoreTable
    @Override
    public void read() {
        psmList = new ArrayList<ScoredString>();
        peptideScoreTable = new HashMap<String, Float>();

        for (MSGFPlusMatch match : msgfPlusPSMList) {
            List<DatabaseMatch> dbMatchList;
            if (considerBestMatchOnly) {
                dbMatchList = new ArrayList<DatabaseMatch>();
                dbMatchList.add(match.getBestDBMatch());
            } else
                dbMatchList = match.getMatchList();

            for (DatabaseMatch m : dbMatchList) {
                String pepSeq = m.getPepSeq();

                boolean isDecoy = true;
                for (int index : m.getIndices()) {
                    String protAcc = sa.getSequence().getAnnotation(index);

                    // Note: By default, decoyProteinPrefix will not end in an underscore
                    // However, if the user defines a custom decoy prefix and they include an underscore, this test will still be valid
                    if (!protAcc.startsWith(decoyProteinPrefix)) {
                        isDecoy = false;
                        break;
                    }
                }

                if (this.isDecoy != isDecoy)
                    continue;

                float specEValue = (float) m.getSpecEValue();
                psmList.add(new ScoredString(pepSeq, specEValue));
                Float prevSpecEValue = peptideScoreTable.get(pepSeq);
                if (prevSpecEValue == null || specEValue < prevSpecEValue)
                    peptideScoreTable.put(pepSeq, specEValue);
            }
        }
    }

}
