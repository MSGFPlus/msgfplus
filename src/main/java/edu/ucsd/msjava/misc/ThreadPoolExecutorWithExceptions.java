package edu.ucsd.msjava.misc;

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
    private String taskName;
    private String progressTitle = "Progress";
    private long startTime;
    private final ScheduledExecutorService statusExecutor = Executors.newSingleThreadScheduledExecutor();
    private final Runnable progressReportRunnable = new Runnable() {
        @Override
        public void run() {
            outputProgressReport();
        }
    };
    private ScheduledFuture<?> currentProgressReportFuture;
    private int progressReportDelayNextChangeMinutes = 0;

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
            if (currentProgressReportFuture == null) {
                currentProgressReportFuture = statusExecutor.scheduleAtFixedRate(progressReportRunnable, 1, 1, TimeUnit.MINUTES);
            }
        }
        super.execute(command);
    }

    @Override
    protected void afterExecute(Runnable r, Throwable t) {
        super.afterExecute(r, t);
        outputProgressReport();
        if (r instanceof ProgressReporter) {
            ProgressReporter reporter = (ProgressReporter) r;
            progressObjects.remove(reporter.getProgressData());
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
        if (r instanceof ProgressReporter) {
            ProgressReporter reporter = (ProgressReporter) r;
            reporter.setProgressData(new ProgressData());
            progressObjects.add(reporter.getProgressData());
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
    
    public void setTaskName(String taskName) {
        this.taskName = taskName;
        this.progressTitle = taskName + " progress";
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
        double timeMinutes = time / 60;
        String units = "seconds";
        if (time > 3600) {
            time = time / 3600;
            units = "hours";
        } else if (time > 60) {
            time = time / 60;
            units = "minutes";
        }
        double totalProgress = progress + getProgressAdjustment();
        System.out.format("%s: %.0f / %.0f tasks, %.2f%%\t\t%.2f %s elapsed%n", this.progressTitle, completed, total, totalProgress, time, units);
        
        if (timeMinutes >= progressReportDelayNextChangeMinutes) {
            ChangeProgressReportDelay();
        }
    }
    
    private void ChangeProgressReportDelay() {
        int nextDelayValue;
        TimeUnit nextDelayUnits;
        switch (progressReportDelayNextChangeMinutes) {
            case 0:
                nextDelayValue = 1;
                nextDelayUnits = TimeUnit.MINUTES;
                progressReportDelayNextChangeMinutes = 60;
                break;
            case 60:
                nextDelayValue = 5;
                nextDelayUnits = TimeUnit.MINUTES;
                progressReportDelayNextChangeMinutes = 180;
                break;
            case 180:
                nextDelayValue = 15;
                nextDelayUnits = TimeUnit.MINUTES;
                progressReportDelayNextChangeMinutes = 600;
                break;
            case 600:
                nextDelayValue = 30;
                nextDelayUnits = TimeUnit.MINUTES;
                progressReportDelayNextChangeMinutes = Integer.MAX_VALUE;
                break;
            default:
                return;
        }
        if (currentProgressReportFuture != null) {
            currentProgressReportFuture.cancel(false);
        }
        currentProgressReportFuture = statusExecutor.scheduleAtFixedRate(progressReportRunnable, nextDelayValue, nextDelayValue, nextDelayUnits);
    }
}
