package msdbsearch;

import msutil.Enzyme;

public class ConcurrentMSGFDB {
	public static class PreProcessSpectra implements Runnable {
		private final ScoredSpectraMap specMap;
		private final int fromIndex;
		private final int toIndex;
		
		public PreProcessSpectra(final ScoredSpectraMap specMap, final int fromIndex, final int toIndex)
		{
			this.specMap = specMap;
			this.fromIndex = fromIndex;
			this.toIndex = toIndex;
		}
		
		@Override
		public void run() 
		{
			specMap.preProcessSpectra(fromIndex, toIndex);
		}
	}

	public static class ComputeSpecProb implements Runnable {
		private final DBScanner scanner;
		private final int fromIndex;
		private final int toIndex;
		private final boolean storeScoreDist;
		
		public ComputeSpecProb(final DBScanner scanner, boolean storeScoreDist, final int fromIndex, final int toIndex)
		{
			this.scanner = scanner;
			this.fromIndex = fromIndex;
			this.toIndex = toIndex;
			this.storeScoreDist = storeScoreDist;
		}
		
		@Override
		public void run() 
		{
			scanner.computeSpecProb(storeScoreDist, fromIndex, toIndex);
		}
	}
	
	public static class RunDBSearch implements Runnable {
		private final DBScanner scanner;
		private final int numberOfAllowableNonEnzymaticTermini;
		private final int fromIndex;
		private final int toIndex;
		private final int searchMode;
		
		public RunDBSearch(final DBScanner scanner, final int numberOfAllowableNonEnzymaticTermini, final int searchMode, final int fromIndex, final int toIndex)
		{
			this.scanner = scanner;
			this.numberOfAllowableNonEnzymaticTermini = numberOfAllowableNonEnzymaticTermini;
			this.fromIndex = fromIndex;
			this.toIndex = toIndex;
			this.searchMode = searchMode;
		}
		
		@Override
		public void run() 
		{
			if(searchMode == 1)
				scanner.dbSearchNoEnzyme(fromIndex, toIndex, true);	// currently not supported
			else if(searchMode == 2)
				scanner.dbSearchCTermEnzymeNoMod(numberOfAllowableNonEnzymaticTermini, fromIndex, toIndex, true);
			else if(searchMode == 3)
				scanner.dbSearchNTermEnzyme(numberOfAllowableNonEnzymaticTermini, fromIndex, toIndex, true);
			else
				scanner.dbSearchCTermEnzyme(numberOfAllowableNonEnzymaticTermini, fromIndex, toIndex, true);
		}
		
	}
}
