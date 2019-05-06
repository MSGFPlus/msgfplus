/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) MSInstrumentInfo.java * Author: * Mathijs Vogelzang
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
 ******************************************************************************/
package org.systemsbiology.jrap.stax;

/**
 * The MSInstrumentInfo class contains information about the
 * MS instrument used for a certain MzXML file.
 * 
 * @author M. Vogelzang
 */
public class MSInstrumentInfo
{
	protected String manufacturer, model, ionization, massAnalyzer, detector;
	protected SoftwareInfo softwareInfo;
	protected MSOperator operator;
	
	/**
	 * @return Returns the detector.
	 */
	public String getDetector()
	{
		return detector;
	}
	
	/**
	 * @return Returns the ionization method.
	 */
	public String getIonization()
	{
		return ionization;
	}
	
	/**
	 * @return Returns the manufacturer.
	 */
	public String getManufacturer()
	{
		return manufacturer;
	}
	/**
	 * @return Returns the mass analyzer.
	 */
	public String getMassAnalyzer()
	{
		return massAnalyzer;
	}
	/**
	 * @return Returns the model.
	 */
	public String getModel()
	{
		return model;
	}
	
	/**
	 * @return Returns the operator.
	 */
	public MSOperator getOperator()
	{
		return operator;
	}
	
	/**
	 * @return Returns the software information.
	 */
	public SoftwareInfo getSoftwareInfo()
	{
		return softwareInfo;
	}

    public String toString()
    {
	return ("msManufacturer "+manufacturer+" msModel "+model+" msIonization "+
		ionization+" msMassAnalyzer "+massAnalyzer+" detector "+detector);
		//+" SoftwareInfo "+softwareInfo.toString());
    }
}
