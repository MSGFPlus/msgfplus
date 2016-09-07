package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.msdbsearch.ConcurrentMSGFPlus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * @author Bryson Gibbons
 */
public class ThreadPoolExecutorWithExceptions extends ThreadPoolExecutor {

    private Throwable thrownData;
    private boolean hasThrownData;

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
    }

    private ThreadPoolExecutorWithExceptions(int corePoolSize, int maximumPoolSize, long keepAliveTime, TimeUnit unit, BlockingQueue<Runnable> workQueue, ThreadFactory threadFactory) {
        super(corePoolSize, maximumPoolSize, keepAliveTime, unit, workQueue, threadFactory);
        thrownData = null;
        hasThrownData = false;
        progressObjects = Collections.synchronizedList(new ArrayList<ProgressData>(maximumPoolSize));
    }

    @Override
    public void afterExecute(Runnable r, Throwable t) {
        super.afterExecute(r, t);
        outputProgressReport();
        if (r instanceof ConcurrentMSGFPlus.RunMSGFPlus) {
            ConcurrentMSGFPlus.RunMSGFPlus run = (ConcurrentMSGFPlus.RunMSGFPlus) r;
            progressObjects.remove(run.getProgressData());
        }
        if (t != null) {
            // store the throwable, to get meaningful data.
            thrownData = t;
            hasThrownData = true;
        }
    }

    @Override
    public void beforeExecute(Thread t, Runnable r) {
        super.beforeExecute(t, r);
        if (r instanceof ConcurrentMSGFPlus.RunMSGFPlus) {
            ConcurrentMSGFPlus.RunMSGFPlus run = (ConcurrentMSGFPlus.RunMSGFPlus) r;
            run.setProgressData(new ProgressData());
            progressObjects.add(run.getProgressData());
        }
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
        System.out.format("Search progress: %.0f / %.0f tasks, %.2f%%%n", completed, total, progress + getProgressAdjustment());
    }
}
