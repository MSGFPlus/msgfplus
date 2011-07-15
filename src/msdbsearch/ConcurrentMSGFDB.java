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
}
