/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) Scan.java * Author: * Ning Zhang
 * nzhang@systemsbiology.org
 * ****************************************************************************** * * *
 * This software is provided ``AS IS'' and any express or implied 
 * warranties, including, but not limited to, the implied warranties of 
 * merchantability and fitness for a particular purpose, are disclaimed.  In
 * no event shall the authors or the Institute for Systems Biology  liable
 * for any direct, indirect, incidental, special, exemplary, or 
 * consequential damages (including, but not limited to, procurement of 
 * substitute goods or services; loss of use, data, or profits; or business
 * interruption) however caused and on any theory of liability, whether in
 * contract, strict liability, or tort (including negligence or otherwise)
 * arising in any way out of the use of this software, even  if advised of
 * the possibility of such damage. * * *
 * ******************************************************************************
 * 
 * redesigned by Ning Zhang @Jan., 2008. 
 * 
 *  
 ******************************************************************************/
package org.systemsbiology.jrap.stax;

import java.io.Serializable;

/**
 * A simple class to hold the contents of a scan from a MSXML file.
 * 
 * This is a start. For those who want to get more fancy you should only have to
 * modify the SAX2ScanHandler to replace this.
 *
 * Note: doubleMassList and doubleIntensityList are separate entities from the components of massIntensityList.
 * dhmay noting this on 2009/03/10 but not touching it, in case there are unknown dependencies on this separation.
 *  
 */
public final class Scan implements Serializable
{
    
    public ScanHeader header;

	/**
	 * A 2-dimensional array, element 0 contains a list of masses of peaks,
	 * element 1 contains a list of intensities of peaks.
	 */
 
    /**
     * No matter 32-bit or 64-bit peak pair, return as double list.
     * Support mzXML
     */
    protected double[][] massIntensityList;
	
    /**
     * No matter 32-bit m/z or 64 bit m/z return as a double list
     * Support for mzML
     */
    public double[] doubleMassList= null;
    public double[] doubleIntensityList=null;
    
    
    public void setHeader(ScanHeader header)
    {
	this.header = header;
    }

    public ScanHeader getHeader()
    {
	return header;
    }

      

    public void setDoubleMassList(double[] newValue)
    {
	doubleMassList = newValue;
    }

    public double[] getDoubleMassList()
    {
	return doubleMassList;
    }

   

    public void setDoubleIntensityList(double[] newValue)
    {
	doubleIntensityList = newValue;
    }

    public double[] getDoubleIntensityList()
    {
	return doubleIntensityList;
    }
    
    //for support mzXML
    public void setMassIntensityList(double[][] massIntensityList)
    {
	this.massIntensityList = massIntensityList;
    }

    public double[][] getMassIntensityList()
    {
	return massIntensityList;
    }

	/**
	 * String respresentation of a Scan object.
	 * 
	 * Note: This is most likely not an optimal way to build the string.
	 * Hopefully this method will only be used for testing.
	 */
    public String toString()
    {
	StringBuffer tmpStrBuffer = new StringBuffer();
	tmpStrBuffer.append("================================================\n");
	tmpStrBuffer.append("peaks:\n");
	      
	if(doubleMassList != null)
	    {
		for(int i=0; i<doubleMassList.length; i++)
		    tmpStrBuffer.append("    mass="+String.format("%.6f",doubleMassList[i])
					+"    intensity="+String.format("%.6f",doubleIntensityList[i])+"\n");
	    }
	    
	else
	    {
		for(int i=0; i<massIntensityList[0].length; i++)
		    tmpStrBuffer.append("    mass="+String.format("%.6f",massIntensityList[0][i])
					+"   intensity="+String.format("%.6f",massIntensityList[1][i])+"\n");
			
	    }

	return (header.toString()+tmpStrBuffer.toString());
    }

}
