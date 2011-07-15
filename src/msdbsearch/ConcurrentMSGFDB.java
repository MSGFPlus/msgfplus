package msdbsearch;

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
}
