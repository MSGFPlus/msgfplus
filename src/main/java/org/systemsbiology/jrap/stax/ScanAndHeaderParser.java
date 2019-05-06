/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) ScanHeaderParser.java * Author: * Ning Zhang
 * nzhang@systemsbiology.org
 * ****************************************************************************** 
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
 * ******************************************************************************/
 
package org.systemsbiology.jrap.stax;
import java.io.*;
import javax.xml.stream.*;
import javax.xml.stream.events.*;
import java.util.*;
import java.nio.ByteBuffer;
import java.util.zip.*;

public class ScanAndHeaderParser{

    public ScanHeader tmpScanHeader;
    public Scan tmpScan;

    
    FileInputStream fileIN = null;

    boolean isScan = false;
 

    public void setIsScan(boolean isScan)
    {
	this.isScan = isScan;
    }

 
   
    public void setFileInputStream(FileInputStream in)
    {
	try{
	    this.fileIN = in;
	}
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
   }

    public ScanHeader getHeader()
    {
	return tmpScanHeader;
    }

    public Scan getScan()
    {
	return tmpScan;
    }

    public void parseScanAndHeader()
    {
        XMLStreamReader xmlSR = null;
	try{
	    XMLInputFactory inputFactory = XMLInputFactory.newInstance();
	    xmlSR = inputFactory.createXMLStreamReader(fileIN,"ISO-8859-1");

	    parseScanAndHeader(xmlSR);
	   
	}
	catch(Exception e)
	    {
		String exception1=e.getMessage();
		if(!exception1.equals("ScanHeaderEndFoundException"))
		    {
			if(!exception1.equals("ScanEndFoundException"))
			    e.printStackTrace();
		    }
	    }
	    finally {
	        if(xmlSR != null) {
	            try {
                    xmlSR.close();
                }
                catch (XMLStreamException e) {
                    e.printStackTrace();
                }
	        }
	    }
    }


    public void parseScanAndHeader(XMLStreamReader xmlSR)
            throws XMLStreamException {
        boolean inPrecursorMZ = false;
	    boolean inPeaks = false;
	    String elementName = null;
	    String attriName = null;
	    String attriValue = null;
	  

	    StringBuffer precursorBuffer = null;
	    StringBuffer peaksBuffer = null;

		boolean isActivationMethodSet = false;
				
	    while(xmlSR.hasNext())
		{
		    int event = xmlSR.next();
		    if(event == xmlSR.START_ELEMENT)
			{
			    elementName = xmlSR.getLocalName();
			    if(elementName.equals("scan"))
				{
				    tmpScanHeader = new ScanHeader();
				    tmpScanHeader.setNum(getIntValue(xmlSR,"num"));
				    tmpScanHeader.setMsLevel(getIntValue(xmlSR, "msLevel"));
				    tmpScanHeader.setPeaksCount(getIntValue(xmlSR, "peaksCount"));
				    tmpScanHeader.setPolarity(getStringValue(xmlSR,"polarity"));
				    tmpScanHeader.setScanType(getStringValue(xmlSR,"scanType"));
				    tmpScanHeader.setCentroided(getIntValue(xmlSR, "centroided"));
				    tmpScanHeader.setDeisotoped(getIntValue(xmlSR, "deisotoped"));
				    tmpScanHeader.setChargeDeconvoluted(getIntValue(xmlSR, "chargeDeconvoluted"));
				    tmpScanHeader.setRetentionTime(getStringValue(xmlSR,"retentionTime"));
				    tmpScanHeader.setStartMz(getFloatValue(xmlSR, "startMz"));
				    tmpScanHeader.setEndMz(getFloatValue(xmlSR, "endMz"));
				    tmpScanHeader.setLowMz(getFloatValue(xmlSR, "lowMz"));
				    tmpScanHeader.setHighMz(getFloatValue(xmlSR, "highMz"));
				    tmpScanHeader.setBasePeakMz(getFloatValue(xmlSR, "basePeakMz"));
				    tmpScanHeader.setBasePeakIntensity(getFloatValue(xmlSR, "basePeakIntensity"));
				    tmpScanHeader.setTotIonCurrent(getFloatValue(xmlSR, "totIonCurrent"));
				    //for S(M)RM
				    tmpScanHeader.setFilterLine(getStringValue(xmlSR,"filterLine"));
					
				    String actMethod = getStringValue(xmlSR,"activationMethod");
				    if(!isActivationMethodSet && actMethod != null && actMethod.length() > 0)
				    {
						isActivationMethodSet = true; 
				    	tmpScanHeader.setActivationMethod(getStringValue(xmlSR,"activationMethod"));
				    }
					
				}
			    if(elementName.equals("peaks"))
				{
				    tmpScanHeader.setPrecision(getIntValue(xmlSR, "precision"));
				    tmpScanHeader.setByteOrder(getStringValue(xmlSR,"byteOrder"));
				    tmpScanHeader.setContentType(getStringValue(xmlSR,"contentType"));
				    tmpScanHeader.setCompressionType(getStringValue(xmlSR,"compressionType"));
				    tmpScanHeader.setCompressedLen(getIntValue(xmlSR, "compressedLen"));

				    if(isScan)
					{
					    inPeaks = true;
					    peaksBuffer = new StringBuffer();
					    tmpScan = new Scan();
					    tmpScan.setHeader(tmpScanHeader);
					}
				    else
					throw new XMLStreamException("ScanHeaderEndFoundException");
				}

			    if(elementName.equals("precursorMz"))
				{
				    tmpScanHeader.setPrecursorScanNum(getIntValue(xmlSR,"precursorScanNum"));
				    tmpScanHeader.setPrecursorCharge(getIntValue(xmlSR,"precursorCharge"));
				    tmpScanHeader.setCollisionEnergy(getFloatValue(xmlSR,"collisionEnergy"));
				    tmpScanHeader.setIonisationEnergy(getFloatValue(xmlSR,"ionisationEnergy"));
				    tmpScanHeader.setPrecursorIntensity(getFloatValue(xmlSR,"precursorIntensity"));
					
				    String actMethod = getStringValue(xmlSR,"activationMethod");
				    if(!isActivationMethodSet && actMethod != null && actMethod.length() > 0)
				    {
					isActivationMethodSet = true; 
				    	tmpScanHeader.setActivationMethod(getStringValue(xmlSR,"activationMethod"));
				    }
					
				    precursorBuffer = new StringBuffer();
				    inPrecursorMZ = true;
				}
			}
		    if(event == xmlSR.CHARACTERS)
			{
			    if(inPrecursorMZ)
				precursorBuffer.append(xmlSR.getText());
			    if(inPeaks)
				peaksBuffer.append(xmlSR.getText());
			}
		    if(event ==xmlSR.END_ELEMENT)
			{
			    elementName = xmlSR.getLocalName();
			    if(elementName.equals("precursorMz"))
				{
				    tmpScanHeader.setPrecursorMz(Float.parseFloat(precursorBuffer.toString()));
                                
				    precursorBuffer = null; // make available for garbage collection

				    inPrecursorMZ = false;
				}
			    if(elementName.equals("peaks"))
				{
				    
				    //get peaks, this time use ByteBuffer
				    getPeaks(peaksBuffer.toString());
				    inPeaks = false;
				    peaksBuffer = null;
				    throw new XMLStreamException("ScanEndFoundException");
				    
				    
				}
			}
		}
    }

    public String getStringValue(XMLStreamReader xmlSR, String name)
    {
	String value="";
	try{
	    if(xmlSR.getAttributeValue(null,name) == null)
		value="";
	    else
		value=xmlSR.getAttributeValue(null,name);
	}
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
	return value;
    }

    public int getIntValue(XMLStreamReader xmlSR, String name)
    {
	int value=-1;
	try{
	    if(xmlSR.getAttributeValue(null,name) == null)
		value = -1;
	    else
		value = Integer.parseInt(xmlSR.getAttributeValue(null,name));
	}
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
	return value;
    }

    public float getFloatValue(XMLStreamReader xmlSR, String name)
    {
	float value=-1f;
	try{
	    if(xmlSR.getAttributeValue(null,name) == null)
		value= -1f;
	    else
		value=Float.parseFloat(xmlSR.getAttributeValue(null,name));
	}
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
	return value;
    }

    public void getPeaks(String peakData)
    {
	//support non-zlib
	byte[] peakArray = peakData.getBytes();
	byte[] outPeakArray = peakArray;
	int outpos = Base64.decode(peakArray,0,peakArray.length,outPeakArray);

	double[][] massIntenList = null;
	int arrayLen = -1;
	ByteBuffer peakBuffer = null;
	//check if it's compressed
	byte[] result=null;
	if((tmpScanHeader.getCompressionType()).equals("zlib"))
	    {
		try{
		    Inflater decompresser = new Inflater();
		    decompresser.setInput(outPeakArray, 0, outpos);
		    int unCompLen = (tmpScanHeader.getPeaksCount())*(tmpScanHeader.getPrecision()/4);
		    result = new byte[unCompLen];
		    decompresser.inflate(result);
		    decompresser.end();

		}
		catch(DataFormatException e)
		    {
			e.printStackTrace();
		    }

		arrayLen = result.length/(tmpScanHeader.getPrecision()/8)/2;
		massIntenList = new double[2][arrayLen];
		peakBuffer = ByteBuffer.wrap(result);
		
	    }
	else
	    {
		arrayLen = outpos/(tmpScanHeader.getPrecision()/8)/2;
		massIntenList = new double[2][arrayLen];
		peakBuffer = ByteBuffer.wrap(outPeakArray,0,outpos);
	    }
	
	int i=0;
	while(peakBuffer.hasRemaining())
	    {
		if(tmpScanHeader.getPrecision() == 32)
		    {
			massIntenList[0][i]=(double)peakBuffer.getFloat();
			massIntenList[1][i]=(double)peakBuffer.getFloat();
		    }
		else
		    {
			massIntenList[0][i] = peakBuffer.getDouble();
			massIntenList[1][i] = peakBuffer.getDouble();
		    }
		i++;
	    }
	    
		
	tmpScan.setMassIntensityList(massIntenList);
    }
}
