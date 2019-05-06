/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) FileHeaderParser.java * Author: * Ning Zhang
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
 * arising in any way out of the use of this software, even if advised of
 * the possibility of such damage. * * *
 * **********************************************************************************/

package org.systemsbiology.jrap.stax;

import java.io.FileInputStream;
import java.util.ArrayList;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

public class FileHeaderParser{

    String inputMZXMLfile;
    MZXMLFileInfo info;

    ArrayList<ParentFile> parentFiles = new ArrayList<ParentFile>();
    ArrayList<SoftwareInfo> dataProcessingSoftware = new ArrayList<SoftwareInfo>();

    boolean isXML = false;
    boolean isML = false;

    FileHeaderParser(String inputMZXMLfile)
    {
	this.inputMZXMLfile = inputMZXMLfile;
	info = new MZXMLFileInfo();

	if(inputMZXMLfile.indexOf("mzXML") != -1)
	    isXML = true;
	else
	    isML = true;
    }

    public MZXMLFileInfo getInfo()
    {
	return info;
    }

    public void parseFileHeader()
    {
	if(isXML)
	    parseXMLFileHeader();
	else
	    parseMLFileHeader();
    }
    
    public void parseFileHeader(XMLStreamReader xmlSR) throws XMLStreamException
    {
    if(isXML)
        parseXMLFileHeader(xmlSR);
    else
        parseMLFileHeader(xmlSR);
    }

    private void parseXMLFileHeader()
    {
	
	try{
	    XMLInputFactory inputFactory = XMLInputFactory.newInstance();
	    XMLStreamReader xmlSR = inputFactory.createXMLStreamReader(new FileInputStream(inputMZXMLfile));

	    parseXMLFileHeader(xmlSR);
	    xmlSR.close();
	}
	catch(Exception e)
	    {
		if(!(e.getMessage()).equals("HeaderEndFoundException"))
		    e.printStackTrace(System.err);
	    }
    }

    private void parseXMLFileHeader(XMLStreamReader xmlSR)
            throws XMLStreamException {
        String elementName = null;
	    int event = -1;

	    boolean isInstrument = false;
	    boolean isDataProcess = false;

	    while(xmlSR.hasNext())
		{
		    event = xmlSR.next();
		    if(event == XMLStreamConstants.START_ELEMENT)
			{
			    elementName = xmlSR.getLocalName();
			    //System.out.println("elementName "+elementName);
			    if(elementName.equals("parentFile"))
				{
				    parentFiles.add(new ParentFile(xmlSR.getAttributeValue(0),
								   xmlSR.getAttributeValue(1),
								   xmlSR.getAttributeValue(2)));
				}
			    if(elementName.equals("msInstrument"))
				{
				    isInstrument = true;
				}
			    if(elementName.equals("msManufacturer"))
				{
				    info.instrumentInfo.manufacturer = xmlSR.getAttributeValue(1);
				}
			    if(elementName.equals("msModel"))
				{
				    info.instrumentInfo.model = xmlSR.getAttributeValue(null,"value");
				}
			    if(elementName.equals("msIonisation"))
				{
				    info.instrumentInfo.ionization = xmlSR.getAttributeValue(null,"value");
				}
			    if(elementName.equals("msMassAnalyzer"))
				{
				    info.instrumentInfo.massAnalyzer = xmlSR.getAttributeValue(null,"value");
				}
			    if(elementName.equals("msDetector"))
				{
				    info.instrumentInfo.detector = xmlSR.getAttributeValue(null,"value");
				}
			    if(elementName.equals("operator"))
				{
				    MSOperator operator = new MSOperator();
				    operator.firstName = xmlSR.getAttributeValue(null,"first");
				    operator.lastName = xmlSR.getAttributeValue(null,"last");
				    operator.phoneNumber = xmlSR.getAttributeValue(null,"phone");
				    operator.email = xmlSR.getAttributeValue(null,"email");
				    operator.URI = xmlSR.getAttributeValue(null,"URI");
				    //System.out.println("operator "+operator);

				    info.instrumentInfo.operator = operator;
				}
			    if(elementName.equals("software"))
				{
				    if(isInstrument)
					{
					    
					    info.instrumentInfo.softwareInfo = parseSoftware(xmlSR);
					}
				    else if(isDataProcess)
					{
					    dataProcessingSoftware.add(parseSoftware(xmlSR));
					}
				}
			    if(elementName.equals("dataProcessing"))
				{
				    isDataProcess = true;
				    
					String value;
					if ((value = xmlSR.getAttributeValue(null,"intensityCutoff")) != null)
					    info.dataProcessing.intensityCutoff = Double.parseDouble(value);

					if ((value = xmlSR.getAttributeValue(null,"centroided")) != null)
					    info.dataProcessing.centroided = Integer.parseInt(value);

					if ((value = xmlSR.getAttributeValue(null,"deisotoped")) != null)
					    info.dataProcessing.deisotoped = Integer.parseInt(value);

					if ((value = xmlSR.getAttributeValue(null,"chargeDeconvoluted")) != null)
					    info.dataProcessing.chargeDeconvoluted = Integer
						.parseInt(value);

					if ((value = xmlSR.getAttributeValue(null,"spotIntegration")) != null)
					    info.dataProcessing.spotIntegration = Integer.parseInt(value);	    
				}
			    
			}
		    if(event == XMLStreamConstants.END_ELEMENT)
			{
			    elementName = xmlSR.getLocalName();
			    if(elementName.equals("msInstrument"))
				isInstrument = false;
			    if(elementName.equals("dataProcessing"))
				{
				    //System.out.println("elementName "+elementName);
				    isDataProcess = false;
				    info.parentFiles = parentFiles;
				    info.dataProcessing.softwareUsed = dataProcessingSoftware;
				    throw new XMLStreamException("HeaderEndFoundException");
				}
			}
		}
    }

    private SoftwareInfo parseSoftware(XMLStreamReader xmlSR)
    {
	SoftwareInfo sInfo = new SoftwareInfo("","","");
	try{
	    if(isXML)
		sInfo.setType(xmlSR.getAttributeValue(null,"type"));
	    else
		sInfo.setType(xmlSR.getAttributeValue(null,"accession"));

	    sInfo.setName(xmlSR.getAttributeValue(null,"name"));
	    sInfo.setVersion(xmlSR.getAttributeValue(null,"version"));
	  
	}
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
	return sInfo;
    }

    private void parseMLFileHeader()
    {

	try{
	    XMLInputFactory inputFactory = XMLInputFactory.newInstance();
	    XMLStreamReader xmlSR = inputFactory.createXMLStreamReader(new FileInputStream(inputMZXMLfile));

	    parseMLFileHeader(xmlSR);
	    xmlSR.close();
	}
	catch(Exception e)
	    {
		if(!(e.getMessage()).equals("HeaderEndFoundException"))
		    e.printStackTrace(System.err);
	    }
    }

    private void parseMLFileHeader(XMLStreamReader xmlSR)
            throws XMLStreamException {
        String elementName = null;
	    int event = -1;

	    boolean isInstrument = false;
	    boolean isDataProcess = false;

	    boolean inSourceFile = false;
	    boolean inSource = false;
	    boolean inAnalyzer = false;
	    boolean inDetector = false;

	    String fileName = null;
	    String fileLocation = null;
	    String accession = null;
	    String fileType = null;
	    String sha1 = null;

	    
	    while(xmlSR.hasNext())
		{
		    event = xmlSR.next();
		    if(event == xmlSR.START_ELEMENT)
			{
			    elementName = xmlSR.getLocalName();
			    //System.out.println("elementName "+elementName);
			    if(elementName.equals("sourceFile"))
				{
				    fileName = xmlSR.getAttributeValue(null,"name");
				    fileLocation = xmlSR.getAttributeValue(null,"location");
				    inSourceFile = true;
				}
			    if(elementName.equals("cvParam"))
				{
				    accession = xmlSR.getAttributeValue(null,"accession");
				    if(inSourceFile)
					{
					    
					    if(accession.equals("MS:1000563"))
						{
						    fileType = xmlSR.getAttributeValue(null,"name");
						}
					    if(accession.equals("MS:1000569"))
						{
						    sha1= xmlSR.getAttributeValue(null,"value");
						}
					}
				    if(isInstrument)
					{
					    if(accession.equals("MS:1000554"))
						{
						   info.instrumentInfo.model = xmlSR.getAttributeValue(null,"name"); 
						}
					    if(accession.equals("MS:1000529"))
						{
						    info.instrumentInfo.manufacturer = xmlSR.getAttributeValue(null,"name")
							+" "+xmlSR.getAttributeValue(null,"value");
						}
					}
				    if(inSource)
					{
					    info.instrumentInfo.ionization = xmlSR.getAttributeValue(null,"name");
					}
				    if(inAnalyzer)
					{
					    info.instrumentInfo.massAnalyzer = xmlSR.getAttributeValue(null,"name");
					}
				    if(inDetector)
					{
					    info.instrumentInfo.detector = xmlSR.getAttributeValue(null,"name");
					}

				    if(isDataProcess)
					{
					    String name, value;
					    name = xmlSR.getAttributeValue(null, "name");
					    value = xmlSR.getAttributeValue(null, "value");
					    if(name.indexOf("deisotoping") != -1)
						{
						    if(value.equals("true"))
							info.dataProcessing.deisotoped = 1;
						    else
							info.dataProcessing.deisotoped = 0;
						}
					    if(name.indexOf("charge") != -1)
						{
						    if(value.equals("true"))
							info.dataProcessing.chargeDeconvoluted = 1;
						    else
							info.dataProcessing.chargeDeconvoluted = 0;
						}
					    if(name.indexOf("peak") != -1)
						{
						    if(value.equals("true"))
							info.dataProcessing.peakPicked = 1;
						    else
							info.dataProcessing.peakPicked = 0;
						}
					    if(name.indexOf("smoothing") != -1)
						{
						    if(value.equals("true"))
							info.dataProcessing.smoothed = 1;
						    else
							info.dataProcessing.smoothed = 0;
						}
					    if(name.indexOf("baseline") != -1)
						{
						    if(value.equals("true"))
							info.dataProcessing.baseLineReduced = 1;
						    else
							info.dataProcessing.baseLineReduced = 0;
						}
					    if(name.indexOf("low intensity") != -1)
						{
						    if(value.equals("true"))
							info.dataProcessing.lowIntensityDataRemoved = 1;
						    else
							info.dataProcessing.lowIntensityDataRemoved = 0;
						}
					}
						    
				}
				
			    if(elementName.equals("referenceableParamGroup"))
				{
				    if((xmlSR.getAttributeValue(null,"id")).indexOf("Instrument") != -1)
					isInstrument = true;
				}
			    if(elementName.equals("source"))
				{
				    inSource = true;
				    
				}
			    if(elementName.equals("analyzer"))
				{
				    inAnalyzer = true;
				    
				}
			    if(elementName.equals("detector"))
				{
				    inDetector = true;
				    
				}
			   
			    if(elementName.equals("softwareParam"))
				{
				    
					
				    dataProcessingSoftware.add(parseSoftware(xmlSR));
					    
				}
			    if(elementName.equals("dataProcessing"))
				{
				    isDataProcess = true;
				    
				}
			    
			}
		    if(event == xmlSR.END_ELEMENT)
			{
			    elementName = xmlSR.getLocalName();
			    if(elementName.equals("sourceFile"))
				{
				    parentFiles.add(new ParentFile(fileLocation+"/"+fileName,
								  fileType,
								  sha1));
				    inSourceFile = false;
				}
			    if(elementName.equals("referenceableParamGroup"))
				{
				    isInstrument = false;
				}
				
			    if(elementName.equals("source"))
				{
				    inSource = false;
				}
			    if(elementName.equals("analyzer"))
				{
				    inAnalyzer = false;
				}
			    if(elementName.equals("detector"))
				{
				    inDetector = false;
				}
			  
			    if(elementName.equals("dataProcessing"))
				{
				    //System.out.println("elementName "+elementName);
				    isDataProcess = false;
				    info.parentFiles = parentFiles;
				    info.dataProcessing.softwareUsed = dataProcessingSoftware;
				    throw new XMLStreamException("HeaderEndFoundException");
				}
			}
		}
    }

}