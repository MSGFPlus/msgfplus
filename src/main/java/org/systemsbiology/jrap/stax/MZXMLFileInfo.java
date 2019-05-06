/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) MZXMLFileInfo.java * Author: * Mathijs Vogelzang
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
 * Created on May 21, 2004
 *
 * Jan. 22, 2008 changed to apply to JDK1.6 by Ning Zhang
 *  
 ******************************************************************************/

package org.systemsbiology.jrap.stax;
import java.util.*;

/**
 * MZXMLFileInfo is a class that contains all information from the header of an
 * MzXML file that is constant for the entire file.
 * 
 * @author M. Vogelzang
 */
public class MZXMLFileInfo
{
    //protected ParentFile[] parentFiles;
    ArrayList<ParentFile> parentFiles = new ArrayList<ParentFile>();
	protected MSInstrumentInfo instrumentInfo;
	protected DataProcessingInfo dataProcessing;

	public MZXMLFileInfo() {
	    //parentFiles = null;
	    instrumentInfo = new MSInstrumentInfo();
		dataProcessing = new DataProcessingInfo();
	}

	/**
	 * Get information about parent files, chronologically ordered.
	 * 
	 * @return An array of information about parent files of an mzXML file.
	 */
	public ArrayList<ParentFile> getParentFiles()
	{
		return parentFiles;
	}

	/**
	 * Get information about the MS instrument used to extract data.
	 * 
	 * @return MS instrument information, or null when no information was
	 *         present in the file.
	 */
	public MSInstrumentInfo getInstrumentInfo()
	{
		return instrumentInfo;
	}

	/**
	 * Get data about how the data was processed.
	 * @return An instance of DataProcessingInfo.
	 */
	public DataProcessingInfo getDataProcessing()
	{
		return dataProcessing;
	}

    public String toString()
    {
	String outputLine="";
	ParentFile pFile = null;
	for(int i=0; i<parentFiles.size();i++)
	    {
		pFile = (parentFiles).get(i);
		outputLine += pFile.toString()+" ";
	    }
	outputLine += instrumentInfo.toString()+" "+dataProcessing.toString();

	return outputLine;
    }
}
