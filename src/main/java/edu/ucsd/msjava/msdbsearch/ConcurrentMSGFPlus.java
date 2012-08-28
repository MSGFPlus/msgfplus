package edu.ucsd.msjava.msdbsearch;

import java.util.List;

public class ConcurrentMSGFPlus {
	public static class RunMSGFPlus implements Runnable {
		private final ScoredSpectraMap specScanner;
		private final DBScanner scanner;
		SearchParams params;
		List<MSGFPlusMatch> resultList;
		
		public RunMSGFPlus(
				ScoredSpectraMap specScanner,
				CompactSuffixArray sa,
				SearchParams params,
				List<MSGFPlusMatch> resultList
				)
		{
			this.specScanner = specScanner;
			this.params = params;
			this.scanner = new DBScanner(
					specScanner, 
					sa, 
					params.getEnzyme(), 
					params.getAASet(), 
					params.getNumMatchesPerSpec(), 
					params.getMinPeptideLength(), 
					params.getMaxPeptideLength());
			this.resultList = resultList;
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
			
			int ntt = params.getNumTolerableTermini();
			if(params.getEnzyme() == null)
				ntt = 2;
			scanner.dbSearch(ntt);
			System.out.print(threadName+": Database search finished ");
			System.out.format("(elapsed time: %.2f sec)\n", (float)((System.currentTimeMillis()-time)/1000));
			
			time = System.currentTimeMillis();
			System.out.println(threadName+": Computing spectral E-values...");
			scanner.computeSpecEValue(false);
			System.out.print(threadName+": Computing spectral E-values finished ");
			System.out.format("(elapsed time: %.2f sec)\n", (float)((System.currentTimeMillis()-time)/1000));
			
			scanner.generateSpecIndexDBMatchMap();
			
			if(params.outputAdditionalFeatures())
				scanner.addAdditionalFeatures();
			
			scanner.addResultsToList(resultList);
//			gen.addSpectrumIdentificationResults(scanner.getSpecIndexDBMatchMap());
		}
	}	
}
