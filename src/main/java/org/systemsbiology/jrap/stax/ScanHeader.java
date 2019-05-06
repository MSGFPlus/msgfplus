/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) ScanHeader.java * Author: * Mathijs Vogelzang
 * m_v@dds.nl
 * ****************************************************************************** * * *
 * This software is provided ``AS IS'' and any express or implied * *
 * warranties, including, but not limited to, the implied warranties of * *
 * merchantability and fitness for a particular purpose, are disclaimed. * * In
 * no event shall the authors or the Institute for Systems Biology * * liable
 * for any direct, indirect, incidental, special, exemplary, or * *
 * consequential damages (including, but not limited to, procurement of * *
 * substitute goods or services; loss of use, data, or profits; or * * business
 * interruption) however caused and on any theory of liability, * * whether in
 * contract, strict liability, or tort (including negligence * * or otherwise)
 * arising in any way out of the use of this software, even * * if advised of
 * the possibility of such damage. * * *
 * ******************************************************************************
 * 
 * ChangeLog
 * 
 * 10-05-2004 Added this header
 * 
 * Created on Jan 12, 2004
 *
 * add support for mzXML_schema_3.0 and S(M)RM by Ning Zhang
 *  
 ******************************************************************************/
package org.systemsbiology.jrap.stax;

import java.io.Serializable;

/**
 * ScanHeader is a class that contains all information
 * associated with a Scan, except for the actual peakList.
 * The separation between the peaklist and the other information
 * was made because parsing the peaklist costs a lot of time, and
 * in this way, programs can parse headers separately, and not parse
 * the peaklist when it's not needed.
 *
 * dhmay: rt and retentionTime are completely separate fields, which is horribly confusing.  Probably getRetentionTime()
 * should form a String around rt, and getDoubleRetentionTime() should just be a cover for getRT(), if both need to
 * exist. Noting this on 2009/03/10 but not touching it, in case there are unknown dependencies on this separation.
 * 
 * @author M. Vogelzang 
 */
public class ScanHeader implements Serializable
{
	//
	// Class Members and Defaults
	// 

	/** Scan Number */
	protected int num = -1;

	/** MS Scan Level */
	protected int msLevel = -1;

	/** Number of peaks in scan */
	protected int peaksCount = -1;

	/** TODO: Describe */
	protected String polarity = null;

	/** TODO: Describe */
	protected String scanType = null;

	/** TODO: Describe */
	protected int centroided = -1;

	/** TODO: Describe */
	protected int deisotoped = -1;

	/** TODO: Describe */
	protected int chargeDeconvoluted = -1;

	/** TODO: Describe */
	protected String retentionTime = null;

    /** for mzML */
    protected double rt = -1;

	/** TODO: Describe */
	protected float startMz = -1;

	/** TODO: Describe */
	protected float endMz = -1;

	/** TODO: Describe */
	protected float lowMz = -1;

	/** TODO: Describe */
	protected float highMz = -1;

	/** TODO: Describe */
	protected float basePeakMz = -1;

	/** TODO: Describe */
	protected float basePeakIntensity = -1;

	/** TODO: Describe */
	protected float totIonCurrent = -1;

	/** TODO: Describe */
	protected float precursorMz = -1;

	/** TODO: Describe */
	protected int precursorScanNum = -1;

	/** TODO: Describe */
	protected int precursorCharge = -1;

    protected float precursorIntensity = -1f;

    // Added by Sangtae
	protected String activationMethod = null;

	/** TODO: Describe */
	protected float collisionEnergy = -1;

	/** TODO: Describe */
	protected float ionisationEnergy = -1;

	/** TODO: Describe */
	protected int precision = -1;
    
    /** for S(M)RM */
    protected String filterLine = null;

    /**Peaks attribute for mzXML_3.0*/
    protected String byteOrder = null;

    protected String contentType = null;

    protected String compressionType = null;

    protected int compressedLen = -1;

    /** for mzML */
    protected int massPrecision = -1;
    protected String massCompressionType = null;
    protected int massCompressedLen = -1;

    protected int intenPrecision = -1;
    protected String intenCompressionType = null;
    protected int intenCompressedLen = -1;

    /**
     * Store the byte offset, within the mz(X)ML file, at which the binary data for this scan are found.
     * dhmay re-adding 20091021.  This was removed by in mid-2008, super important.  Note: this must be set explicitly
     * by calling code -- the offset won't be found in the scan XML itself */
    protected long scanOffset = -1;

    /**
	 * @return Returns the basePeakIntensity.
	 */
	public float getBasePeakIntensity()
	{
		return basePeakIntensity;
	}

	/**
	 * @param basePeakIntensity
	 *            The basePeakIntensity to set.
	 */
	public void setBasePeakIntensity(float basePeakIntensity)
	{
		this.basePeakIntensity = basePeakIntensity;
	}

	/**
	 * @return Returns the basePeakMz.
	 */
	public float getBasePeakMz()
	{
		return basePeakMz;
	}

	/**
	 * @param basePeakMz
	 *            The basePeakMz to set.
	 */
	public void setBasePeakMz(float basePeakMz)
	{
		this.basePeakMz = basePeakMz;
	}

    /**
     *@return returns the byteOrder
     */
    public String getByteOrder()
    {
	return byteOrder;
    }

    /**
     *@param byteOrder set the byteOrder
     */
    public void setByteOrder(String byteOrder)
    {
	this.byteOrder = byteOrder;
    }

  
	/**
	 * @return Returns the centroided.
	 */
	public int getCentroided()
	{
		return centroided;
	}

	/**
	 * @param centroided
	 *            The centroided to set.
	 */
	public void setCentroided(int centroided)
	{
		this.centroided = centroided;
	}

	/**
	 * @return Returns the chargeDeconvoluted.
	 */
	public int getChargeDeconvoluted()
	{
		return chargeDeconvoluted;
	}

	/**
	 * @param chargeDeconvoluted
	 *            The chargeDeconvoluted to set.
	 */
	public void setChargeDeconvoluted(int chargeDeconvoluted)
	{
		this.chargeDeconvoluted = chargeDeconvoluted;
	}

	/**
	 * @return Returns the collisionEnergy.
	 */
	public float getCollisionEnergy()
	{
		return collisionEnergy;
	}

	/**
	 * @param collisionEnergy
	 *            The collisionEnergy to set.
	 */
	public void setCollisionEnergy(float collisionEnergy)
	{
		this.collisionEnergy = collisionEnergy;
	}

    /**
     *@return returns compressionType
     */
    public String getCompressionType()
    {
	return compressionType;
    }

    /**
     *@param compressionType set compressionType
     */
    public void setCompressionType(String compressionType)
    {
	this.compressionType = compressionType;
    }

    /**
     *@return returns compressedLen
     */
    public int getCompressedLen()
    {
	return compressedLen;
    }

    /**
     *@param compressedLen set compressedLen
     */
    public void setCompressedLen(int compressedLen)
    {
	this.compressedLen = compressedLen;
    }

    /**
     *@return returns the contentType
     */
    public String getContentType()
    {
	return contentType;
    }

    /**
     *@param contentType set the contentType
     */
    public void setContentType(String contentType)
    {
	this.contentType = contentType;
    }

	/**
	 * @return Returns the deisotoped.
	 */
	public int getDeisotoped()
	{
		return deisotoped;
	}

	/**
	 * @param deisotoped
	 *            The deisotoped to set.
	 */
	public void setDeisotoped(int deisotoped)
	{
		this.deisotoped = deisotoped;
	}

	/**
	 * @return Returns the endMz.
	 */
	public float getEndMz()
	{
		return endMz;
	}

	/**
	 * @param endMz
	 *            The endMz to set.
	 */
	public void setEndMz(float endMz)
	{
		this.endMz = endMz;
	}

    /**
     *@return Returns the filterLine.
     */
    public String getFilterLine()
    {
	return filterLine;
    }

    /**
     *@param filterLine
     */

    public void setFilterLine(String filterLine)
    {
	this.filterLine = filterLine;
    }

	/**
	 * @return Returns the highMz.
	 */
	public float getHighMz()
	{
		return highMz;
	}

	/**
	 * @param highMz
	 *            The highMz to set.
	 */
	public void setHighMz(float highMz)
	{
		this.highMz = highMz;
	}

	/**
	 * @return Returns the ionisationEnergy.
	 */
	public float getIonisationEnergy()
	{
		return ionisationEnergy;
	}

	/**
	 * @param ionisationEnergy
	 *            The ionisationEnergy to set.
	 */
	public void setIonisationEnergy(float ionisationEnergy)
	{
		this.ionisationEnergy = ionisationEnergy;
	}

	/**
	 * @return Returns the lowMz.
	 */
	public float getLowMz()
	{
		return lowMz;
	}

	/**
	 * @param lowMz
	 *            The lowMz to set.
	 */
	public void setLowMz(float lowMz)
	{
		this.lowMz = lowMz;
	}

	/**
	 * @return Returns the msLevel.
	 */
	public int getMsLevel()
	{
		return msLevel;
	}

	/**
	 * @param msLevel
	 *            The msLevel to set.
	 */
	public void setMsLevel(int msLevel)
	{
		this.msLevel = msLevel;
	}

	/**
	 * @return Returns the num.
	 */
	public int getNum()
	{
		return num;
	}

	/**
	 * @param num
	 *            The num to set.
	 */
	public void setNum(int num)
	{
		this.num = num;
	}

	/**
	 * @return Returns the peaksCount.
	 */
	public int getPeaksCount()
	{
		return peaksCount;
	}

	/**
	 * @param peaksCount
	 *            The peaksCount to set.
	 */
	public void setPeaksCount(int peaksCount)
	{
		this.peaksCount = peaksCount;
	}

	/**
	 * @return Returns the polarity.
	 */
	public String getPolarity()
	{
		return polarity;
	}

	/**
	 * @param polarity
	 *            The polarity to set.
	 */
	public void setPolarity(String polarity)
	{
		this.polarity = polarity;
	}

	/**
	 * @return Returns the precision.
	 */
	public int getPrecision()
	{
		return precision;
	}

	/**
	 * @param precision
	 *            The precision to set.
	 */
	public void setPrecision(int precision)
	{
		this.precision = precision;
	}

	/**
	 * @return Returns the precursorCharge.
	 */
	public int getPrecursorCharge()
	{
		return precursorCharge;
	}

	/**
	 * @param precursorCharge
	 *            The precursorCharge to set.
	 */
	public void setPrecursorCharge(int precursorCharge)
	{
		this.precursorCharge = precursorCharge;
	}

	/**
	 * @return Returns the precursorMz.
	 */
	public float getPrecursorMz()
	{
		return precursorMz;
	}

	/**
	 * @param precursorMz
	 *            The precursorMz to set.
	 */
	public void setPrecursorMz(float precursorMz)
	{
		this.precursorMz = precursorMz;
	}

	/**
	 * @return Returns the precursorScanNum.
	 */
	public int getPrecursorScanNum()
	{
		return precursorScanNum;
	}

	/**
	 * @param precursorScanNum
	 *            The precursorScanNum to set.
	 */
	public void setPrecursorScanNum(int precursorScanNum)
	{
		this.precursorScanNum = precursorScanNum;
	}

    /**
     * @return Returns the  precursorIntensity
     *
     */
    public float getPrecursorIntensity()
    {
	return precursorIntensity;
    }

    /**
     * @param precursorIntensity
     *         The precursorIntensity to set
     */
    public void setPrecursorIntensity(float precursorIntensity)
    {
	this.precursorIntensity = precursorIntensity;
    }

    	/**
	 * @return Returns the activaion method
	 */
    	public String getActivationMethod()
	{
		return activationMethod;
	}

	/**
	 * @param activationMethod
	 * 	The activation method to be set
	 */
	public void setActivationMethod(String activationMethod)
	{
		this.activationMethod = activationMethod;
	}

	/**
	 * @return Returns the retentionTime.
	 */
	public String getRetentionTime()
	{
		return retentionTime;
	}

	/**
	 * @param retentionTime
	 *            The retentionTime to set.
	 */
	public void setRetentionTime(String retentionTime)
	{
		this.retentionTime = retentionTime;
	}


    /**
     * 
     *             The retentionTime for mzML.
     */
    public double getRT()
    {
	return rt;
    }

    public void setRT(double rt)
    {
	this.rt = rt;
    }

	/**
	 * @return Returns the scanType.
	 */
	public String getScanType()
	{
		return scanType;
	}

	/**
	 * @param scanType
	 *            The scanType to set.
	 */
	public void setScanType(String scanType)
	{
		this.scanType = scanType;
	}

	/**
	 * @return Returns the startMz.
	 */
	public float getStartMz()
	{
		return startMz;
	}

	/**
	 * @param startMz
	 *            The startMz to set.
	 */
	public void setStartMz(float startMz)
	{
		this.startMz = startMz;
	}

	/**
	 * @return Returns the totIonCurrent.
	 */
	public float getTotIonCurrent()
	{
		return totIonCurrent;
	}

	/**
	 * @param totIonCurrent
	 *            The totIonCurrent to set.
	 */
	public void setTotIonCurrent(float totIonCurrent)
	{
		this.totIonCurrent = totIonCurrent;
	}

	public double getDoubleRetentionTime()
	{
        // TODO: more robust ISO time conversion?
		if (retentionTime.charAt(0) != 'P'
			|| retentionTime.charAt(1) != 'T'
			|| retentionTime.charAt(retentionTime.length() - 1) != 'S')
		{
			throw new IllegalArgumentException(
				"Format of retentiontime is not PTxxxxS, don't know how to parse "
					+ retentionTime);
		}

		return Double.parseDouble(
			retentionTime.substring(2, retentionTime.length() - 1));
	}

    /**
     * @return massPrecision
     */
    public int getMassPrecision()
    {
	return massPrecision;
    }

    /**
     * @param massPrecision to set
     */

    public void setMassPrecision(int massPrecision)
    {
	this.massPrecision = massPrecision;
    }

    /**
     * @return massCompressionType
     */
    public String getMassCompressionType()
    {
	return massCompressionType;
    }

    /**
     * @param massCompressionType
     */
    public void setMassCompressionType(String massCompressionType)
    {
	this.massCompressionType = massCompressionType;
    }

    /**
     *@return massCompressionLen
     */

    public int getMassCompressedLen()
    {
	return massCompressedLen;
    }

    /**
     *@param massCompressedLen
     */
    public void setMassCompressedLen(int massCompressedLen)
    {
	this.massCompressedLen = massCompressedLen;
    }

    /**
     *@return intenPrecision
     */

    public int getIntenPrecision()
    {
	return intenPrecision;
    }

    /**
     *@param intenPrecision
     */
    public void setIntenPrecision(int intenPrecision)
    {
	this.intenPrecision = intenPrecision;
    }

    /**
     *@return intenCompressionType
     */

    public String getIntenCompressionType()
    {
	return intenCompressionType;
    }
    
    /**
     *@param intenCompressionType
     */

    public void setIntenCompressionType(String intenCompressionType)
    {
	this.intenCompressionType = intenCompressionType;
    }

    /**
     *@return intenCompressedLen
     */
    public int getIntenCompressedLen()
    {
	return intenCompressedLen;
    }

    /**
     *@param intenCompressedLen
     */
    public void setIntenCompressedLen(int intenCompressedLen)
    {
	this.intenCompressedLen = intenCompressedLen;
    }
	
	/**
	 * String respresentation of a ScanHeader object.
	 * 
	 * Note: This is most likely not an optimal way to build the string.
	 * Hopefully this method will only be used for testing.
	 */
	public String toString()
	{
		StringBuffer tmpStrBuffer = new StringBuffer(1000);
		tmpStrBuffer.append("SCANHEADER\n");
		tmpStrBuffer.append("==========\n");
		tmpStrBuffer.append("num = " + num + "\n");
		tmpStrBuffer.append("msLevel = " + msLevel + "\n");
		tmpStrBuffer.append("peaksCount = " + peaksCount + "\n");
		tmpStrBuffer.append("polarity = " + polarity + "\n");
		tmpStrBuffer.append("scanType = " + scanType + "\n");
		tmpStrBuffer.append("centroided = " + centroided + "\n");
		tmpStrBuffer.append("deisotoped = " + deisotoped + "\n");
		tmpStrBuffer.append(
				"chargeDeconvoluted = " + chargeDeconvoluted + "\n");
		tmpStrBuffer.append("retentionTime = " + retentionTime + "\n");
		tmpStrBuffer.append("startMz = " + startMz + "\n");
		tmpStrBuffer.append("endMz = " + endMz + "\n");
		tmpStrBuffer.append("lowMz = " + lowMz + "\n");
		tmpStrBuffer.append("highMz = " + highMz + "\n");
		tmpStrBuffer.append("basePeakMz = " + basePeakMz + "\n");
		tmpStrBuffer.append("basePeakIntensity = " + basePeakIntensity + "\n");
		tmpStrBuffer.append("totIonCurrent = " + totIonCurrent + "\n");
		tmpStrBuffer.append("precursorMz = " + precursorMz + "\n");
		tmpStrBuffer.append("precursorScanNum = " + precursorScanNum + "\n");
		tmpStrBuffer.append("precursorCharge = " + precursorCharge + "\n");
		tmpStrBuffer.append("precursorIntensity = " + precursorIntensity +"\n");
		tmpStrBuffer.append("collisionEnergy = " + collisionEnergy + "\n");
		tmpStrBuffer.append("ionisationEnergy = " + ionisationEnergy + "\n");
		tmpStrBuffer.append("precision = " + precision + "\n");
		//add for mzXML_3.0
		tmpStrBuffer.append("byteOrder = "+ byteOrder + "\n");
		tmpStrBuffer.append("contentType = " + contentType + "\n");
		tmpStrBuffer.append("compressionType = " + compressionType + "\n");
		tmpStrBuffer.append("compressedLen = " + compressedLen + "\n");

		//for mzML
		tmpStrBuffer.append("rt "+rt+"\n");
		tmpStrBuffer.append("massPrecision "+massPrecision+"\n");
		tmpStrBuffer.append("massCompressionType "+massCompressionType+"\n");
		tmpStrBuffer.append("massCompressedLen "+massCompressedLen+"\n");
		tmpStrBuffer.append("intenPrecision "+intenPrecision+"\n");
		tmpStrBuffer.append("intenCompressionType "+intenCompressionType+"\n");
		tmpStrBuffer.append("intenCompressedLen "+intenCompressedLen+"\n");

		//add for support S(M)RM
		tmpStrBuffer.append("filterLine = "+filterLine+"\n");

		return (tmpStrBuffer.toString());
	}

    public long getScanOffset()
    {
        return scanOffset;
    }

    public void setScanOffset(long scanOffset)
    {
        this.scanOffset = scanOffset;
    }
}
