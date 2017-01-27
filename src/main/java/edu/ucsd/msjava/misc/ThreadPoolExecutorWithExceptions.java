package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.msdbsearch.ConcurrentMSGFPlus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.*;

/**
 * @author Bryson Gibbons
 */
public class ThreadPoolExecutorWithExceptions extends ThreadPoolExecutor {

    private Throwable thrownData;
    private boolean hasThrownData;
    private long startTime;
    private final ScheduledExecutorService statusExecutor = Executors.newSingleThreadScheduledExecutor();
    private final Runnable progressReportRunnable = new Runnable() {
        @Override
        public void run() {
            outputProgressReport();
        }
    };
    private ScheduledFuture<?> currentProgressReportFuture;

    private final List<ProgressData> progressObjects;

    public static ThreadPoolExecutorWithExceptions newFixedThreadPool(int nThreads) {
        return new ThreadPoolExecutorWithExceptions(nThreads, nThreads,
                0L, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<Runnable>());
    }

    public static ThreadPoolExecutorWithExceptions newFixedThreadPool(int nThreads, ThreadFactory threadFactory) {
        return new ThreadPoolExecutorWithExceptions(nThreads, nThreads,
                0L, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<Runnable>(),
                threadFactory);
    }

    private ThreadPoolExecutorWithExceptions(int corePoolSize, int maximumPoolSize, long keepAliveTime, TimeUnit unit, BlockingQueue<Runnable> workQueue) {
        super(corePoolSize, maximumPoolSize, keepAliveTime, unit, workQueue, Executors.defaultThreadFactory());
        thrownData = null;
        hasThrownData = false;
        progressObjects = Collections.synchronizedList(new ArrayList<ProgressData>(maximumPoolSize));
        startTime = -1;
    }

    private ThreadPoolExecutorWithExceptions(int corePoolSize, int maximumPoolSize, long keepAliveTime, TimeUnit unit, BlockingQueue<Runnable> workQueue, ThreadFactory threadFactory) {
        super(corePoolSize, maximumPoolSize, keepAliveTime, unit, workQueue, threadFactory);
        thrownData = null;
        hasThrownData = false;
        progressObjects = Collections.synchronizedList(new ArrayList<ProgressData>(maximumPoolSize));
        startTime = -1;
    }
    
    @Override
    public void execute(Runnable command) {
        if (startTime < 0) {
            startTime = System.currentTimeMillis();
            currentProgressReportFuture = statusExecutor.scheduleAtFixedRate(progressReportRunnable, 1, 1, TimeUnit.MINUTES);
        }
        super.execute(command);
    }

    @Override
    protected void afterExecute(Runnable r, Throwable t) {
        super.afterExecute(r, t);
        outputProgressReport();
        if (r instanceof ConcurrentMSGFPlus.RunMSGFPlus) {
            ConcurrentMSGFPlus.RunMSGFPlus run = (ConcurrentMSGFPlus.RunMSGFPlus) r;
            progressObjects.remove(run.getProgressData());
        }
        if (t != null && thrownData == null) {
            // store the throwable, to get meaningful data.
            thrownData = t;
            hasThrownData = true;
            this.shutdownNow();
        }
    }

    @Override
    protected void beforeExecute(Thread t, Runnable r) {
        super.beforeExecute(t, r);
        if (r instanceof ConcurrentMSGFPlus.RunMSGFPlus) {
            ConcurrentMSGFPlus.RunMSGFPlus run = (ConcurrentMSGFPlus.RunMSGFPlus) r;
            run.setProgressData(new ProgressData());
            progressObjects.add(run.getProgressData());
        }
    }

    @Override
    public boolean awaitTermination(long timeout, TimeUnit unit) throws InterruptedException {
        boolean result = false;
        InterruptedException except = null;
        try {
            result = super.awaitTermination(timeout, unit);
        } catch (InterruptedException e) {
            except = e;
        }

        // Shutdown the progress reporting
        currentProgressReportFuture.cancel(true);
        statusExecutor.shutdown();

        // Return/throw the original result
        if (except != null)
        {
            throw except;
        }
        return result;
    }

    public boolean awaitTerminationWithExceptions(long timeout, TimeUnit unit) throws Throwable {
        boolean result = false;
        InterruptedException interrupted = null;
        try {
            result = this.awaitTermination(timeout, unit);
        } catch (InterruptedException e) {
            interrupted = e;
        }

        // If we have data thrown by a thread, throw that instead of the result of awaitTermination
        if (hasThrownData) {
            throw thrownData;
        }

        // No data thrown by a thread? Return/throw the original result
        if (interrupted != null) {
            throw interrupted;
        }
        return result;
    }

    public boolean HasThrownData() {
        return hasThrownData;
    }

    public Throwable getThrownData() {
        return thrownData;
    }

    /*
    * Get the adjustment value for progress reporting
    */
    public double getProgressAdjustment() {
        double count = 0.0;
        double progressSum = 0.0;
        synchronized (progressObjects) {
            for (ProgressData data : progressObjects) {
                count += 1;
                progressSum += data.getProgress();
            }
        }
        if (count < 1) {
            // No active tasks, prevent divide by zero
            return 0.0;
        }
        double progress = progressSum / count;
        double weight = count / this.getTaskCount();
        return progress * weight;
    }

    /*
    * Output a progress report to the console
    */
    public void outputProgressReport() {
        double completed = getCompletedTaskCount();
        double total = getTaskCount();
        if (total < 1) {
            // prevent divide by zero - should never be zero (unless someone rearranges code), but here just in case.
            total = 1;
        }
        double progress = (completed / total) * 100.0;
        
        double time = (System.currentTimeMillis() - startTime) / 1000.0;
        String units = "seconds";
        if (time > 3600) {
            time = time / 3600;
            units = "hours";
        } else if (time > 60) {
            time = time / 60;
            units = "minutes";
        }
        double totalProgress = progress + getProgressAdjustment();
        System.out.format("Search progress: %.0f / %.0f tasks, %.2f%%\t\t%.2f %s elapsed%n", completed, total, totalProgress, time, units);
    }
}
