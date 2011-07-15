package msdbsearch;

import java.util.List;

import msutil.SpecKey;

public class ConcurrentMSGFDB {
	public static class PreProcessSpectra implements Runnable {

		private final ScoredSpectraMap specMap;
		private final List<SpecKey> specKeyList;
		public PreProcessSpectra(final ScoredSpectraMap specMap, final List<SpecKey> specKeyList)
		{
			this.specMap = specMap;
			this.specKeyList = specKeyList;
		}
		
		@Override
		public void run() 
		{
			specMap.preProcessSpectra(specKeyList);
		}
	}
}
