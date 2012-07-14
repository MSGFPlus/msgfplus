package jmzparser;

import java.io.File;
import java.util.List;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import uk.ac.ebi.jmzml.model.mzml.CV;
import uk.ac.ebi.jmzml.model.mzml.CVList;
import uk.ac.ebi.jmzml.model.mzml.FileDescription;
import uk.ac.ebi.jmzml.model.mzml.MzML;
import uk.ac.ebi.jmzml.xml.Constants;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

public class TestMzMLReader {
	public static void main(String argv[]) throws Exception
	{
		Logger.getRootLogger().setLevel(Level.FATAL);
		test();
	}
	
	public static void test() throws Exception
	{
		File xmlFile = new File("/Users/kims336/Research/Data/JMzReader/tiny.pwiz.mzML");
		MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(xmlFile);
		MzML completeMzML = unmarshaller.unmarshall();
		
		System.out.println("MzML Version: = " + unmarshaller.getMzMLVersion());
		System.out.println("MzML ID: = " + unmarshaller.getMzMLId());
		System.out.println("MzML Accession: = " + unmarshaller.getMzMLAccession());		
		
		CVList cvList = unmarshaller.unmarshalFromXpath("/cvList", CVList.class);
		//the object is now fully populated with the data from the XML file
		System.out.println("Number of defined CVs in this mzML: " + cvList.getCount());
		//retrieve the fileDescription element
		FileDescription fd = unmarshaller.unmarshalFromXpath("/fileDescription", FileDescription.class);
		System.out.println("Number of source files: " + fd.getSourceFileList().getCount());
		
		//supported XPath for indexed and non-indexed mzML 
		System.out.println("Supported XPath:" + Constants.XML_INDEXED_XPATHS);
	}
}
