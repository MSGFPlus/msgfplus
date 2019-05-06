/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 ***************************************************************************/

/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) MLScanHeaderParser.java * Author: * Ning Zhang
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
 
//support mzML parsing
package org.systemsbiology.jrap.stax;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

/**
 * dhmay changing 2009/10/21, incorporating Vagisha's changes and rebuilding my changes from 03/2009.
 * mzML 1.1 changes the way scan IDs are stored.  They are now stored in
 * the "id" attribute of "spectrum", which is being used to contain multiple name-value pairs; the
 * name of the name-value pair containing the scan number is "scan", so I'm knocking off everything but that pair.
 * Also changing to cobble together Scan.massIntensityList from its components, which was missed earlier.
 * Also calling tmpScanHeader.setRetentionTime(), which was previously not set.
 *
 * F Levander changing 2010-02-11. Changed parsing of CVparams from names to accession numbers since those are stable.
 */
public class MLScanAndHeaderParser
{

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
        try
        {
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

    /**
     * dhmay: mzML 1.1 changes the way scan IDs are stored.  They are now stored in
     * the "id" attribute of "spectrum", which is being used to contain multiple name-value pairs; the
     * name of the name-value pair containing the scan number is "scan", so I'm knocking off everything but that pair.
     * @param idString
     * @return The scan number or if a numeric value couldn't be parsed.
     */
     protected int parseScanNumberFromSpectrumIdField(String idString)
     {
    	 int retval=-1;
         if (idString.contains("scan="))
             idString = idString.substring(idString.indexOf("scan=") + "scan=".length());
         if (idString.contains("scanId="))
             idString = idString.substring(idString.indexOf("scanId=") + "scanId=".length());
         if (idString.contains(" "))
             idString = idString.substring(0, idString.indexOf(" "));
         try
         {
        	 retval=Integer.parseInt(idString);
         }
         catch(Exception e)
         {
             e.printStackTrace();
         }
         return retval;
     }


    public void parseMLScanAndHeader()
    {
        XMLStreamReader xmlSR = null;
        try{
            XMLInputFactory inputFactory = XMLInputFactory.newInstance();
            xmlSR = inputFactory.createXMLStreamReader(fileIN,"ISO-8859-1");

            parseMLScanAndHeader(xmlSR);

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
        finally
        {
            if(xmlSR != null)
            {
                try
                {
                    xmlSR.close();
                }
                catch (XMLStreamException e)
                {
                    e.printStackTrace();
                }
            }
        }
    }



    public void parseMLScanAndHeader(XMLStreamReader xmlSR)
            throws XMLStreamException
    {
        boolean inSpectrum = false;
        boolean inPeaks = false;
        boolean isPeaks = false;
        String elementName = null;
        String attriName = null;
        String attriValue = null;

        StringBuffer peaksBuffer = null;
        int count = 0;

        while(xmlSR.hasNext())
        {
            int event = xmlSR.next();
            if(event == xmlSR.START_ELEMENT)
            {
                elementName = xmlSR.getLocalName();
                if(elementName.equals("spectrum"))
                {
                    inSpectrum = true;
                    count=0;
                    tmpScanHeader = new ScanHeader();
                    //dhmay changing 2009/03/09.  mzML 1.1 changes the way scan IDs are stored
                    tmpScanHeader.setNum(parseScanNumberFromSpectrumIdField(getStringValue(xmlSR, "id")));
                    // If scan number couldn't be parsed, fall back to index+1.
                    if (tmpScanHeader.getNum()==-1) 
                    	tmpScanHeader.setNum(getIntValue(xmlSR, "index")+1);
                    tmpScanHeader.setPeaksCount(getIntValue(xmlSR, "defaultArrayLength"));

                }
                if(elementName.equals("cvParam"))
                {
                    attriName = xmlSR.getAttributeValue(null,"name");
                    String attriAccession = xmlSR.getAttributeValue(null,"accession");
                    if(inSpectrum)
                    {

                        if(attriAccession.equals("MS:1000511"))
                            tmpScanHeader.setMsLevel(getIntValue(xmlSR, "value"));
                        if(attriAccession.equals("MS:1000127"))
                            tmpScanHeader.setCentroided(1);
                        if(attriAccession.equals("MS:1000504"))
                            tmpScanHeader.setBasePeakMz(getFloatValue(xmlSR, "value"));
                        if(attriAccession.equals("MS:1000505"))
                            tmpScanHeader.setBasePeakIntensity(getFloatValue(xmlSR, "value"));
                        if(attriAccession.equals("MS:1000285"))
                            tmpScanHeader.setTotIonCurrent(getFloatValue(xmlSR, "value"));
                        if(attriAccession.equals("MS:1000528"))
                            tmpScanHeader.setStartMz(getFloatValue(xmlSR, "value"));
                        if(attriAccession.equals("MS:1000527"))
                            tmpScanHeader.setEndMz(getFloatValue(xmlSR, "value"));
                        if(attriAccession.equals("MS:1000501"))
                            tmpScanHeader.setLowMz(getFloatValue(xmlSR, "value"));
                        if(attriAccession.equals("MS:1000500" ))
                            tmpScanHeader.setHighMz(getFloatValue(xmlSR, "value"));
                        if(attriAccession.equals("MS:1000512"))
                            tmpScanHeader.setFilterLine(getStringValue(xmlSR,"value"));
                        if(attriAccession.equals("MS:1000498"))
                            tmpScanHeader.setScanType("full scan");
                        if(attriAccession.equals("MS:1000130"))
                            tmpScanHeader.setPolarity("+");
                        //dhmay changed this for mzML 1.1.0RC5,a nd then again for RC6.
                        //Hopefully the name of this attribute will settle down.
                        if(attriAccession.equals("MS:1000016"))
                        {
                            String timeType = xmlSR.getAttributeValue(null,"unitName");
                            double rt = Double.parseDouble(xmlSR.getAttributeValue(null, "value"));
                            // flevander changed 
                            if(timeType.equals("minute"))
                                rt = rt * 60;
                            tmpScanHeader.setRT(rt);

                            //dhmay adding for backward compatibility.  Probably this should be rewired so that
                            //getDoubleRetentionTime just accesses the rt variable, but I don't want to sort out
                            //that tangle
                            tmpScanHeader.setRetentionTime("PT" + rt + "S");                            
                        }
                        // Precursor m/z from isolation window target m/z or selected ion m/z
                        // Selected m7z comes afterwards and will have precedence
                        if(attriAccession.equals("MS:1000827") || attriAccession.equals("MS:1000744"))
                            tmpScanHeader.setPrecursorMz(getFloatValue(xmlSR,"value"));
                        if(attriAccession.equals("MS:1000042"))
                            tmpScanHeader.setPrecursorIntensity(getFloatValue(xmlSR,"value"));
                        if(attriAccession.equals("MS:1000041"))
                            tmpScanHeader.setPrecursorCharge(getIntValue(xmlSR,"value"));
                        if(attriAccession.equals("MS:1000045"))
                            tmpScanHeader.setCollisionEnergy(getFloatValue(xmlSR,"value"));
                    }
                    if(inPeaks)
                    {
                        if(attriAccession.equals("MS:1000523"))
                        {
                            if(count == 1)
                                tmpScanHeader.setMassPrecision(64);
                            if(count == 2)
                                tmpScanHeader.setIntenPrecision(64);
                        }
                        if(attriAccession.equals("MS:1000521"))
                        {
                            if(count == 1)
                                tmpScanHeader.setMassPrecision(32);
                            if(count == 2)
                                tmpScanHeader.setIntenPrecision(32);
                        }
                        if(attriAccession.equals("MS:1000576"))
                        {
                            if(count == 1)
                            {
                                tmpScanHeader.setMassCompressionType("None");
                            }
                            if(count == 2)
                            {
                                tmpScanHeader.setIntenCompressionType("None");
                            }
                        }
                        if(attriAccession.equals("MS:1000574"))
                        {
                            if(count == 1)
                            {
                                tmpScanHeader.setMassCompressionType("zlib");
                            }
                            if(count == 2)
                            {
                                tmpScanHeader.setIntenCompressionType("zlib");
                            }
                        }
                    }


                }

                if(elementName.equals("binaryDataArrayList"))
                {
                    if(isScan)
                    {
                        tmpScan = new Scan();
                        tmpScan.setHeader(tmpScanHeader);
                    }
                    else
                        throw new XMLStreamException("ScanHeaderEndFoundException");
                }

                if(elementName.equals("binaryDataArray"))
                {
                    inPeaks = true;
                    count++;
                    //System.out.println("count "+count);

                    if(count == 1)
                        tmpScanHeader.setMassCompressedLen(getIntValue(xmlSR, "encodedLength"));
                    if(count == 2)
                        tmpScanHeader.setIntenCompressedLen(getIntValue(xmlSR, "encodedLength"));


                }
                if(elementName.equals("binary"))
                {
                    inPeaks = false;
                    isPeaks = true;
                    peaksBuffer = new StringBuffer();
                }
            }
            if(event == xmlSR.CHARACTERS)
            {
                if(isPeaks)
                    peaksBuffer.append(xmlSR.getText());
            }
            if(event ==xmlSR.END_ELEMENT)
            {
                elementName = xmlSR.getLocalName();
                if(elementName.equals("spectrumDescription"))
                {
                    inSpectrum = false;
                }

                if(elementName.equals("binary"))
                {
                    getPeaks(peaksBuffer.toString(),count);
                    isPeaks = false;
                    peaksBuffer = null;
                }
                if(elementName.equals("binaryDataArrayList"))
                {
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

    public void getPeaks(String peakData, int count)
    {
        int precision = -1;
        if(count == 1)
        {
            precision = tmpScanHeader.getMassPrecision();
        }
        if(count == 2)
        {
            precision = tmpScanHeader.getIntenPrecision();
        }
        //support non-zlib
        byte[] peakArray = peakData.getBytes();
        byte[] outPeakArray = peakArray;
        int outpos = Base64.decode(peakArray,0,peakArray.length,outPeakArray);


        double[] doubleMassList = null;

        double[] doubleIntenList = null;

        ByteBuffer peakBuffer = null;
        //check if it's compressed
        byte[] result=null;
        int unCompLen = outpos;
        String compressType = "None";
        if(count == 1)
            compressType = tmpScanHeader.getMassCompressionType();
        if(count == 2)
            compressType = tmpScanHeader.getIntenCompressionType();

        if(compressType.equals("zlib"))
        {
            try{
                Inflater decompresser = new Inflater();
                decompresser.setInput(outPeakArray, 0, outpos);
                unCompLen = (tmpScanHeader.getPeaksCount())*(precision/4);
                result = new byte[unCompLen];
                decompresser.inflate(result);
                decompresser.end();

            }
            catch(DataFormatException e)
            {
                e.printStackTrace();
            }

            peakBuffer = ByteBuffer.wrap(result);
            peakBuffer.order(ByteOrder.LITTLE_ENDIAN);

        }
        else
        {
            peakBuffer = ByteBuffer.wrap(outPeakArray,0,outpos);
            peakBuffer.order(ByteOrder.LITTLE_ENDIAN);
        }

        if(precision == 64)
        {
            if(count == 1)
            {
                doubleMassList = new double[unCompLen/8];
                int i=0;
                while(peakBuffer.hasRemaining())
                {
                    doubleMassList[i] = peakBuffer.getDouble();
                    i++;
                }
                tmpScan.setDoubleMassList(doubleMassList);
                //System.out.println("massList size "+tmpScan.getDoubleMassList().length);

            }

            if(count == 2)
            {
                doubleIntenList = new double[unCompLen/8];
                int i=0;
                while(peakBuffer.hasRemaining())
                {
                    doubleIntenList[i] = peakBuffer.getDouble();
                    i++;
                }
                tmpScan.setDoubleIntensityList(doubleIntenList);
            }
        }
        else
        {

            if(count == 1)
            {
                doubleMassList = new double[unCompLen/4];
                int i=0;
                while(peakBuffer.hasRemaining())
                {
                    doubleMassList[i] = (double)peakBuffer.getFloat();
                    i++;
                }
                tmpScan.setDoubleMassList(doubleMassList);
            }
            if(count == 2)
            {
                doubleIntenList = new double[unCompLen/4];
                int i=0;
                while(peakBuffer.hasRemaining())
                {
                    doubleIntenList[i] = (double)peakBuffer.getFloat();
                    i++;
                }
                tmpScan.setDoubleIntensityList(doubleIntenList);
                //System.out.println("intenList size "+tmpScan.getFloatIntensityList().length);
            }
        }
        //dhmay fixing up the massIntensityList, 2009/03/09.  This seems to have been missed initially
        if (count == 2)
        {
//System.err.println("****Setting mass-int list");
            double[][] massIntensityList = new double[2][];
            massIntensityList[0] = tmpScan.getDoubleMassList();
            massIntensityList[1] = tmpScan.getDoubleIntensityList();
            tmpScan.setMassIntensityList(massIntensityList);
        }        
    }


}
