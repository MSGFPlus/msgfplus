package edu.ucsd.msjava.misc;

/**
 *
 * @author bryson
 */
public class ProgressData {
    private double progress;
    private double minPercent;
    private double maxPercent;
    private ProgressData parentProgress;
    
    public ProgressData() {
        progress = 0.0;
        minPercent = 0;
        maxPercent = 100;
        isPartialRange = false;
        parentProgress = null;
    }
    
    public ProgressData(ProgressData parent) {
        progress = 0.0;
        minPercent = 0;
        maxPercent = 100;
        isPartialRange = false;
        parentProgress = parent;
    }
    
    public void setParentProgressObj(ProgressData progressObj) {
        parentProgress = progressObj;
    }
    
    public ProgressData getParentProgressObj() {
        return parentProgress;
    }
    
    public void resetProgress() {
        progress = 0.0;
    }
    
    private void setProgress(double pct) {
        progress = pct;
    }
    
    public double getProgress() {
        if (isPartialRange) {
            return progress * ((maxPercent - minPercent) / 100) + minPercent;
        }
        return progress;
    }
    
    public boolean isPartialRange;
    
    public void setMinPercentage(double pct) {
        checkSetMinMaxRange(pct, maxPercent);
    }
    
    public double getMinPercentage(double pct) {
        return minPercent;
    }
    
    public void setMaxPercentage(double pct) {
        checkSetMinMaxRange(minPercent, pct);
    }
    
    public double getMaxPercentage(double pct) {
        return maxPercent;
    }
    
    public void stepRange(double newMaxPct) {
        if (!isPartialRange) {
            isPartialRange = true;
            
            minPercent = 0;
            if (maxPercent >= 100) {
                maxPercent = 0;
            }
        }
        checkSetMinMaxRange(maxPercent, newMaxPct);
    }

    private void checkSetMinMaxRange(double minPct, double maxPct) {
        boolean partial = isPartialRange;
        double pct = progress;
        progress = pct;
        isPartialRange = false;

        if (maxPct > minPct) {
            minPercent = minPct;
            maxPercent = maxPct;
        }
        if (minPercent < 0) {
            minPercent = 0;
        }
        if (maxPercent > 100.0) {
            maxPercent = 100;
        }

        isPartialRange = partial;

        if (partial) {
            // Trigger a report so the data is correct
            report(0.0);
        }
    }

    public void updateProgress(double pct) {
        setProgress(pct);
    }
    
    public void report(double pct) {
        setProgress(pct);
        if (parentProgress != null) {
            parentProgress.report(this.getProgress());
        }
        // report to callable?
    }
    
    public void reportDecimal(double pct) {
        report(pct * 100.0);
    }
    
    public void report(double count, double total) {
        reportDecimal(count / total);
    }
}
