package edu.ucsd.msjava.msdbsearch;

import java.util.List;


import edu.ucsd.msjava.msgf.MSGFDBResultGenerator;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;

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
		
		public void run() 
		{
			if(searchMode == 1)
				scanner.dbSearch(2, fromIndex, toIndex, true);
			else if(searchMode == 2)
				scanner.dbSearch(numberOfAllowableNonEnzymaticTermini, fromIndex, toIndex, true);
			else if(searchMode == 3)
				scanner.dbSearch(numberOfAllowableNonEnzymaticTermini, fromIndex, toIndex, true);
			else
				scanner.dbSearch(numberOfAllowableNonEnzymaticTermini, fromIndex, toIndex, true);
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
		
		public void run() 
		{
			scanner.computeSpecEValue(storeScoreDist, fromIndex, toIndex);
		}
	}
	
	public static class RunMSGFDB implements Runnable {
		private final ScoredSpectraMap specScanner;
		private final DBScanner scanner;
		private final int numberOfAllowableNonEnzymaticTermini;
		private final int searchMode;
		private final boolean storeScoreDist;
		private final String specFileName;
		private final List<MSGFDBResultGenerator.DBMatch> gen;
		private final boolean replicateMergedResults;
		
		public RunMSGFDB(
				ScoredSpectraMap specScanner,
				CompactSuffixArray sa,
				Enzyme enzyme,
				AminoAcidSet aaSet,
				int numPeptidesPerSpec,
				int minPeptideLength,
				int maxPeptideLength,
				int numberOfAllowableNonEnzymaticTermini, 
				boolean storeScoreDist,
				List<MSGFDBResultGenerator.DBMatch> gen, 
				String specFileName,
				boolean replicateMergedResults
				)
		{
			this.specScanner = specScanner;
			this.scanner = new DBScanner(specScanner, sa, enzyme, aaSet, numPeptidesPerSpec, minPeptideLength, maxPeptideLength);
			this.numberOfAllowableNonEnzymaticTermini = numberOfAllowableNonEnzymaticTermini;
			this.storeScoreDist = storeScoreDist;
			this.specFileName = specFileName;
			this.gen = gen;
			this.replicateMergedResults = replicateMergedResults;
			
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
		
		public void run() 
		{
			String threadName = Thread.currentThread().getName();
			
			// Pre-process spectra
			long time = System.currentTimeMillis();
			if(specScanner.getPepMassSpecKeyMap().size() == 0)
				specScanner.makePepMassSpecKeyMap();
			System.out.println(threadName+": Preprocessing spectra...");
			specScanner.preProcessSpectra();
			System.out.print(threadName+": Preprocessing spectra finished ");
			System.out.format("(elapsed time: %.2f sec)\n", (float)((System.currentTimeMillis()-time)/1000));
			
			time = System.currentTimeMillis();
			// DB search
			System.out.println(threadName+": Database search...");
			scanner.setThreadName(threadName);
			if(searchMode == 1)
				scanner.dbSearchNoEnzyme(true);
			else if(searchMode == 2)
				scanner.dbSearchCTermEnzymeNoMod(numberOfAllowableNonEnzymaticTermini, true);
			else if(searchMode == 3)
				scanner.dbSearchNTermEnzyme(numberOfAllowableNonEnzymaticTermini, true);
			else
				scanner.dbSearchCTermEnzyme(numberOfAllowableNonEnzymaticTermini, true);
			System.out.print(threadName+": Database search finished ");
			System.out.format("(elapsed time: %.2f sec)\n", (float)((System.currentTimeMillis()-time)/1000));
			
			time = System.currentTimeMillis();
			System.out.println(threadName+": Computing spectral probabilities...");
			scanner.computeSpecEValue(storeScoreDist);
			System.out.print(threadName+": Computing spectral probabilities finished ");
			System.out.format("(elapsed time: %.2f sec)\n", (float)((System.currentTimeMillis()-time)/1000));
			
			scanner.addDBSearchResults(gen, specFileName, replicateMergedResults);
		}
	}	
	
	public static class RunMSGFDBLib implements Runnable {
		private final ScoredSpectraMap specScanner;
		private final LibraryScanner scanner;
		private final String specFileName;
		private final List<MSGFDBResultGenerator.DBMatch> gen;
		private final String libraryFileName;
		
		public RunMSGFDBLib(
				ScoredSpectraMap specScanner,
				int numPeptidesPerSpec,
				List<MSGFDBResultGenerator.DBMatch> gen, 
				String specFileName,
				String libraryFileName
				)
		{
			this.specScanner = specScanner;
			this.scanner = new LibraryScanner(specScanner, numPeptidesPerSpec);
			this.specFileName = specFileName;
			this.gen = gen;
			this.libraryFileName = libraryFileName;
		}
		
		public void run() 
		{
			String threadName = Thread.currentThread().getName();
			
			// Pre-process spectra
			long time = System.currentTimeMillis();
			if(specScanner.getPepMassSpecKeyMap().size() == 0)
				specScanner.makePepMassSpecKeyMap();
			System.out.println(threadName+": Preprocessing spectra...");
			specScanner.preProcessSpectra();
			System.out.print(threadName+": Preprocessing spectra finished ");
			System.out.format("(elapsed time: %.2f sec)\n", (float)((System.currentTimeMillis()-time)/1000));
			
			time = System.currentTimeMillis();
			
			// Library search
			System.out.println(threadName+": Library search...");
			scanner.setThreadName(threadName);
			scanner.libSearch(libraryFileName, true);
			System.out.print(threadName+": Library search finished ");
			System.out.format("(elapsed time: %.2f sec)\n", (float)((System.currentTimeMillis()-time)/1000));

			// Computing spectral probabilities
			time = System.currentTimeMillis();
			System.out.println(threadName+": Computing spectral probabilities...");
			scanner.computeSpecProb();
			System.out.print(threadName+": Computing spectral probabilities finished ");
			System.out.format("(elapsed time: %.2f sec)\n", (float)((System.currentTimeMillis()-time)/1000));
			
			scanner.addLibSearchResults(gen, specFileName);
		}
	}		
}
