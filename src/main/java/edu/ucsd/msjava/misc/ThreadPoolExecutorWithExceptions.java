package edu.ucsd.msjava.misc;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author Bryson Gibbons
 */
public class ThreadPoolExecutorWithExceptions extends ThreadPoolExecutor {

    private Throwable thrownData;
    private boolean hasThrownData;

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
    }

    private ThreadPoolExecutorWithExceptions(int corePoolSize, int maximumPoolSize, long keepAliveTime, TimeUnit unit, BlockingQueue<Runnable> workQueue, ThreadFactory threadFactory) {
        super(corePoolSize, maximumPoolSize, keepAliveTime, unit, workQueue, threadFactory);
        thrownData = null;
        hasThrownData = false;
    }

    @Override
    public void afterExecute(Runnable r, Throwable t) {
        super.afterExecute(r, t);
        if (t != null) {
            // store the throwable, to get meaningful data.
            thrownData = t;
            hasThrownData = true;
        }
    }

    public boolean HasThrownData() {
        return hasThrownData;
    }

    public Throwable getThrownData() {
        return thrownData;
    }
}
