package edu.ucsd.msjava.msdbsearch;

import java.util.List;


import edu.ucsd.msjava.mzid.MZIdentMLGen;

public class ConcurrentMSGFPlus {
	public static class RunMSGFPlus implements Runnable {
		private final ScoredSpectraMap specScanner;
		private final DBScanner scanner;
		SearchParams params;
		MZIdentMLGen gen;
		
		public RunMSGFPlus(
				ScoredSpectraMap specScanner,
				CompactSuffixArray sa,
				SearchParams params,
				MZIdentMLGen gen 
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
			this.gen = gen;
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
			scanner.computeSpecProb(false);
			System.out.print(threadName+": Computing spectral E-values finished ");
			System.out.format("(elapsed time: %.2f sec)\n", (float)((System.currentTimeMillis()-time)/1000));
			
			gen.addSpectrumIdentificationResults(scanner.getSpecKeyDBMatchMap());
		}
	}	
}
