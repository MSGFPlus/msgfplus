package msdbsearch;

import msgf.MSGFDBResultGenerator;
import msutil.AminoAcidSet;
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
	
	public static class RunMSGFDB implements Runnable {
		private final ScoredSpectraMap specScanner;
		private final DBScanner scanner;
		private final int numberOfAllowableNonEnzymaticTermini;
		private final int searchMode;
		private final boolean storeScoreDist;
		private final String specFileName;
		private final MSGFDBResultGenerator gen;
		
		public RunMSGFDB(
				ScoredSpectraMap specScanner,
				SuffixArrayForMSGFDB sa,
				Enzyme enzyme,
				AminoAcidSet aaSet,
				int numPeptidesPerSpec,
				int minPeptideLength,
				int maxPeptideLength,
				int numberOfAllowableNonEnzymaticTermini, 
				boolean storeScoreDist,
				MSGFDBResultGenerator gen, 
				String specFileName
				)
		{
			this.specScanner = specScanner;
			this.scanner = new DBScanner(specScanner, sa, enzyme, aaSet, numPeptidesPerSpec, minPeptideLength, maxPeptideLength);
			this.numberOfAllowableNonEnzymaticTermini = numberOfAllowableNonEnzymaticTermini;
			this.storeScoreDist = storeScoreDist;
			this.specFileName = specFileName;
			this.gen = gen;
			
			int searchMode = 0;
			if(enzyme == null || enzyme.getResidues() == null)
				searchMode = 1;	
			else if(enzyme.isCTerm())
			{
				if(!aaSet.containsModification())
					searchMode = 2;
				else
					searchMode = 0;
			}
			else
				searchMode = 3;
			this.searchMode = searchMode;
			
		}
		
		@Override
		public void run() 
		{
			String threadName = Thread.currentThread().getName();
			
			// Pre-process spectra
			long time = System.currentTimeMillis();
			if(specScanner.getPepMassSpecKeyMap().size() == 0)
				specScanner.makePepMassSpecKeyMap();
			specScanner.preProcessSpectra();
			System.out.print(threadName+": Preprocessing spectra... ");
			System.out.format("%.3f sec\n", (float)((System.currentTimeMillis()-time)/1000));
			
			time = System.currentTimeMillis();
			// DB search
			if(searchMode == 1)
				scanner.dbSearchNoEnzyme(true);
			else if(searchMode == 2)
				scanner.dbSearchCTermEnzymeNoMod(numberOfAllowableNonEnzymaticTermini, true);
			else if(searchMode == 3)
				scanner.dbSearchNTermEnzyme(numberOfAllowableNonEnzymaticTermini, true);
			else
				scanner.dbSearchCTermEnzyme(numberOfAllowableNonEnzymaticTermini, true);
			System.out.print(threadName+": Database search... ");
			System.out.format("%.3f sec\n", (float)((System.currentTimeMillis()-time)/1000));
			
			time = System.currentTimeMillis();
			scanner.computeSpecProb(storeScoreDist);
			System.out.print(threadName+": Computing spectral probabilities... ");
			System.out.format("%.3f sec\n", (float)((System.currentTimeMillis()-time)/1000));
			
			time = System.currentTimeMillis();
			scanner.addDBSearchResults(gen, specFileName);
			System.out.print(threadName+": Generating results... ");
			System.out.format("%.3f sec\n", (float)((System.currentTimeMillis()-time)/1000));
		}
	}	
}
