package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.misc.ProgressData;

import java.util.List;

public class ConcurrentMSGFPlus {
    public static class RunMSGFPlus implements Runnable {
        private final ScoredSpectraMap specScanner;
        private final DBScanner scanner;
        SearchParams params;
        List<MSGFPlusMatch> resultList;
        private final int taskNum;
        private ProgressData progress;

        public void setProgressData(ProgressData data) {
            progress = data;
        }

        public ProgressData getProgressData() {
            return progress;
        }

        public RunMSGFPlus(
                ScoredSpectraMap specScanner,
                CompactSuffixArray sa,
                SearchParams params,
                List<MSGFPlusMatch> resultList,
                int taskNum
        ) {
            this.specScanner = specScanner;
            this.params = params;
            this.scanner = new DBScanner(
                    specScanner,
                    sa,
                    params.getEnzyme(),
                    params.getAASet(),
                    params.getNumMatchesPerSpec(),
                    params.getMinPeptideLength(),
                    params.getMaxPeptideLength(),
                    params.getMaxNumVariatsPerPeptide(),
                    params.getMinDeNovoScore(),
                    params.ignoreMetCleavage()
            );
            this.resultList = resultList;
            this.taskNum = taskNum;
            progress = null;
        }

        @Override
        public void run() {
            if (progress == null) {
                progress = new ProgressData();
            }
            progress.stepRange(5.0);
            String threadName = Thread.currentThread().getName();
            System.out.println(threadName + ": Starting task " + taskNum);

            specScanner.setProgressObj(new ProgressData(progress));

            // Pre-process spectra
            long time = System.currentTimeMillis();
            if (specScanner.getPepMassSpecKeyMap().size() == 0)
                specScanner.makePepMassSpecKeyMap();
            System.out.println(threadName + ": Preprocessing spectra...");
            specScanner.preProcessSpectra();
            System.out.print(threadName + ": Preprocessing spectra finished ");
            System.out.format("(elapsed time: %.2f sec)\n", (float) ((System.currentTimeMillis() - time) / 1000));

            specScanner.getProgressObj().setParentProgressObj(null);
            progress.report(5.0);
            progress.stepRange(80.0);
            scanner.setProgressObj(new ProgressData(progress));

            time = System.currentTimeMillis();
            // DB search
            System.out.println(threadName + ": Database search...");
            scanner.setThreadName(threadName);

            int ntt = params.getNumTolerableTermini();
            if (params.getEnzyme() == null)
                ntt = 0;
            int nnet = 2 - ntt;
            scanner.dbSearch(nnet);
            System.out.print(threadName + ": Database search finished ");
            System.out.format("(elapsed time: %.2f sec)\n", (float) ((System.currentTimeMillis() - time) / 1000));

            progress.stepRange(95.0);

            time = System.currentTimeMillis();
            System.out.println(threadName + ": Computing spectral E-values...");
            scanner.computeSpecEValue(false);
            System.out.print(threadName + ": Computing spectral E-values finished ");
            System.out.format("(elapsed time: %.2f sec)\n", (float) ((System.currentTimeMillis() - time) / 1000));

            scanner.getProgressObj().setParentProgressObj(null);
            progress.stepRange(100);

            scanner.generateSpecIndexDBMatchMap();

            progress.report(30.0);

            if (params.outputAdditionalFeatures())
                scanner.addAdditionalFeatures();

            progress.report(60.0);

            scanner.addResultsToList(resultList);

            progress.report(100.0);
//			gen.addSpectrumIdentificationResults(scanner.getSpecIndexDBMatchMap());
            System.out.println(threadName + ": Task " + taskNum + " completed.");
        }
    }
}
