/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) ScanHeaderParser.java * Author: * Vagisha Sharma
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

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Map;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

/**
 * 
 */
public class MSXMLSequentialParser {

    /** The file we are in charge of reading */
    private String fileName = null;
    private XMLStreamReader xmlSR = null;
    private InputStream inputStr = null;


    private MZXMLFileInfo fileHeader = null;

    /** The indexes */
    private Map<Integer, Long> offsets;
    private int maxScan;
    private long chrogramIndex;


    private boolean isXML = false;
    private boolean isML = false;

    private int currentScan = 0; // current scan number being read
    
    public MSXMLSequentialParser() {}

    public void open(String fileName) throws FileNotFoundException, XMLStreamException {
        this.fileName = fileName;

        if(fileName.indexOf("mzXML") != -1)
            isXML = true;
        else {
            isML = true;
        }

        //using IndexParser get indexes
        IndexParser indexParser = new IndexParser(fileName); // this will open and close the file once. 
        indexParser.parseIndexes();
        offsets = indexParser.getOffsetMap();
        maxScan = indexParser.getMaxScan();
        chrogramIndex = indexParser.getChrogramIndex();

        inputStr = new FileInputStream(fileName);
        XMLInputFactory inputFactory = XMLInputFactory.newInstance();
        xmlSR = inputFactory.createXMLStreamReader(inputStr);

        // read the file header 
        try {
            readFileHeader(xmlSR);
        }
        catch(XMLStreamException e) {
            if(!(e.getMessage()).equals("HeaderEndFoundException")) {
                throw e;
            }
        }
    }

    public void close() {
        if(this.xmlSR != null) {
            try {xmlSR.close();}
            catch (XMLStreamException e) {}
        }
        if(this.inputStr != null) {
            try {inputStr.close();}
            catch(IOException e) {}
        }
    }

    /**this gives back the file header (info before scan)
     *@return the file header info (MZXMLFileInfo)
     * @throws XMLStreamException 
     */
    private void readFileHeader(XMLStreamReader reader) throws XMLStreamException {
        FileHeaderParser fileHeaderParser = new FileHeaderParser(fileName);
        fileHeaderParser.parseFileHeader(reader);
        this.fileHeader = fileHeaderParser.getInfo();
    }
    
    /**
     * This gives back the file header (info before scan)
     *@return the file header info (MZXMLFileInfo)
     */
    public MZXMLFileInfo getFileHeader() {
        return this.fileHeader;
    }


    /**
     * Returns true if there are more scans to be parsed in the file
     * @return
     */
    public boolean hasNextScan() {
        return currentScan != maxScan;
    }
    
    /**
     * Returns a Scan object with its peaks and header information
     * @return
     * @throws XMLStreamException 
     */
    public Scan getNextScan() throws XMLStreamException {
        if(isXML)
        {
            ScanAndHeaderParser scanParser = new ScanAndHeaderParser();
            scanParser.setIsScan(true);
            try {
                scanParser.parseScanAndHeader(xmlSR);
            }
            catch(XMLStreamException e) {
                if(!e.getMessage().equals("ScanHeaderEndFoundException") && 
                   !e.getMessage().equals("ScanEndFoundException")) {
                    throw e;
                }
            }
            this.currentScan = scanParser.getScan().getHeader().getNum();
            return scanParser.getScan();
        }
        else
        {
            MLScanAndHeaderParser scanParser = new MLScanAndHeaderParser();
            scanParser.setIsScan(true);
            try {
                scanParser.parseMLScanAndHeader(xmlSR);
            }
            catch(XMLStreamException e) {
                if(!e.getMessage().equals("ScanHeaderEndFoundException") && 
                   !e.getMessage().equals("ScanEndFoundException")) {
                    throw e;
                }
            }
            this.currentScan = scanParser.getScan().getHeader().getNum();
            return (scanParser.getScan());
        }
    }
    

    /**
     * Get the total number of scans in the mzXMLfile handled by this parser.
     *
     * @return The number of scans.
     */
    public int getScanCount()
    {
        return offsets.size();
    }

    public int getMaxScanNumber()
    {
        return maxScan;
    }
   
}