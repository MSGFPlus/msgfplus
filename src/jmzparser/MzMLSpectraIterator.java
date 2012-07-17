//package jmzparser;
//
//import java.io.File;
//import java.io.FileNotFoundException;
//import java.util.Collections;
//import java.util.Iterator;
//import java.util.Map;
//import java.util.Map.Entry;
//
//import org.apache.log4j.Level;
//import org.apache.log4j.Logger;
//
//import msutil.ActivationMethod;
//import msutil.Peak;
//import msutil.SpecFileFormat;
//
//import uk.ac.ebi.jmzml.model.mzml.*;
//import uk.ac.ebi.jmzml.model.mzml.params.ActivationCVParam;
//import uk.ac.ebi.jmzml.xml.Constants;
//import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
//import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;
//import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
//import uk.ac.ebi.pride.tools.ms2_parser.Ms2File;
//import uk.ac.ebi.pride.tools.mzdata_parser.MzDataFile;
//import uk.ac.ebi.pride.tools.mzml_wrapper.MzMlWrapper;
//import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
//import uk.ac.ebi.pride.tools.pkl_parser.PklFile;
//
//public class MzMLSpectraIterator implements Iterator<msutil.Spectrum>, Iterable<msutil.Spectrum> {
//	private MzMLUnmarshaller unmarshaller;
//	private MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> itr;
//	private int specIndex = 0;
//	
//	public MzMLSpectraIterator(File specFile) throws FileNotFoundException
//	{
//		unmarshaller = new MzMLUnmarshaller(specFile);
//		itr = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", uk.ac.ebi.jmzml.model.mzml.Spectrum.class);
//	}
//	
//	@Override
//	public boolean hasNext() 
//	{
//		return itr.hasNext();
//	}
//
//	@Override
//	public msutil.Spectrum next() 
//	{
//		uk.ac.ebi.jmzml.model.mzml.Spectrum jmzSpec = itr.next();
//		return null;
////		msutil.Spectrum spec = JmzSpectrumConverter.convert(jmzSpec);
////		spec.setSpecIndex(++specIndex);
////		
////		return spec;
//	}
//
//	@Override
//	public void remove() 
//	{
//		throw new UnsupportedOperationException("SpectraIterator.remove() not implemented");
//	}
//	
//	@Override
//	public Iterator<msutil.Spectrum> iterator() {
//		return this;
//	}
//	
//	public static void main(String argv[]) throws Exception
//	{
//		Logger.getRootLogger().setLevel(Level.FATAL);
//		test();
//	}
//	
//	public static void test() throws Exception
//	{
//		File xmlFile = new File("/Users/kims336/Research/Data/JMzReader/example.mzML");
//		xmlFile = new File("/Users/kims336/Research/Data/JMzReader/small.pwiz.1.1.mzML");
//		
//		MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(xmlFile);
////		MzML completeMzML = unmarshaller.unmarshall();
//		
////		CVList cvList = unmarshaller.unmarshalFromXpath("/cvList", CVList.class);
////		//the object is now fully populated with the data from the XML file
////		System.out.println("Number of defined CVs in this mzML: " + cvList.getCount());
////		for(CV cv : cvList.getCv())
////		{
////			System.out.println(cv.getId()+"\t"+cv.getFullName());
////		}
////		
//		//supported XPath for indexed and non-indexed mzML 
////		System.out.println("Number of Supported XPath:" + Constants.XML_INDEXED_XPATHS.size());
////		for(String xpath : Constants.XML_INDEXED_XPATHS)
////			System.out.println(xpath);
//		
//		MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> itr = 
//				unmarshaller.unmarshalCollectionFromXpath(
//						"/run/spectrumList/spectrum", 
//						uk.ac.ebi.jmzml.model.mzml.Spectrum.class);
//		int numSpecs = 0;
//		while(itr.hasNext())
//		{
//			uk.ac.ebi.jmzml.model.mzml.Spectrum jmzSpec = itr.next();
//			System.out.println("ID: " + jmzSpec.getId());
//			System.out.println("Index: " + jmzSpec.getIndex());	// 0-based index
//			PrecursorList precursorList = jmzSpec.getPrecursorList();
//			if(precursorList != null)
//			{
//				for(Precursor precursor : precursorList.getPrecursor())
//				{
//					// set-up the activation method
////					ParamGroup paramGroup = precursor.getActivation();
////					System.out.println("CVParams");
////					for(CVParam param : paramGroup.getCvParam())
////					{
////						System.out.println("AC?: " + (param instanceof ActivationCVParam));
////						System.out.println(param.getName()+"\t|\t"+param.getValue()+"\t|\t"+param.getHid());
////						System.out.println(param.getAccession());
////						System.out.println(param.getCvRef());
////						System.out.println(param.getName());
////						System.out.println(param.getUnitAccession());
////						System.out.println(param.getUnitCvRef());
////						System.out.println(param.getUnitName());
////						System.out.println(param.getValue());
////					}
////					System.out.println("UserParams");
////					for(UserParam param : paramGroup.getUserParam())
////					{
////						System.out.println(param.getName()+"\t|\t"+param.getValue());
////					}
//					
//					for(ParamGroup paramGroup : precursor.getSelectedIonList().getSelectedIon())
//					{
//						for(CVParam param : paramGroup.getCvParam())
//						{
//							System.out.println(param.getAccession()+"\t"+param.getValue());
//						}
//						for(UserParam param : paramGroup.getUserParam())
//						{
//							System.out.println(param.getName()+"\t"+param.getValue());
//						}
//					}
//				}
//				System.out.println("NumPrecursors: " + precursorList.getPrecursor().size());
//				System.exit(1);
//			}
//		}
//	}
//	
//	public static msutil.Spectrum getSpectrumFromJMzMLSpec(uk.ac.ebi.jmzml.model.mzml.Spectrum jmzMLSpec)
//	{
//		msutil.Spectrum spec = new msutil.Spectrum();
//		
//		PrecursorList precursorList = jmzMLSpec.getPrecursorList();
//		if(precursorList != null && precursorList.getPrecursor().size() > 0)
//		{
//			Precursor precursor = precursorList.getPrecursor().get(0);	// consider only the first precursor
//			
//			// precursor mz, charge
//			float precursorMz = 0;
//			int precursorCharge = 0;
//			float precursorIntensity = 0;
//			
//			if(precursor.getSelectedIonList() != null)
//			{
//				for(ParamGroup paramGroup : precursor.getSelectedIonList().getSelectedIon())
//				{
//					for(CVParam param : paramGroup.getCvParam())
//					{
//						if(param.getAccession().equals("MS:1000744"))	// selected ion m/z
//						{
//							precursorMz = Float.parseFloat(param.getValue());	// assume that unit is m/z (MS:1000040)
//							break;
//						}
//						else if(param.getAccession().equals("MS:1000633"))	// possible charge state
//						{
//							if(precursorCharge == 0)
//								precursorCharge = Integer.parseInt(param.getValue());
//							else
//								precursorCharge = 0;	// if 2 or more charges are possible, set charge to zero
//						}
//					}
//					if(precursorMz != 0)
//						break;
//				}
//			}
//
//			// TODO: read precursor intensity
//			spec.setPrecursor(new Peak(precursorMz, precursorIntensity, precursorCharge));
//			
//			// activation method
//			ParamGroup actMethodParams = precursor.getActivation();
//			for(CVParam param : actMethodParams.getCvParam())
//			{
//				ActivationMethod am = ActivationMethod.getByCV(param.getAccession());
//				if(am != null)
//				{
//					spec.setActivationMethod(am);
//					break;
//				}
//			}
//			
//		}
//
//		// setup peak list
//		Map<Double, Double> peakList = jmzSpec.getPeakList();
//		Iterator<Entry<Double, Double>> itr = peakList.entrySet().iterator();
//		while(itr.hasNext())
//		{
//			Entry<Double, Double> entry = itr.next();
//			Double mz = entry.getValue();
//			Double intensity = entry.getKey();
//			if(mz != null && intensity != null)
//				spec.add(new Peak(mz.floatValue(), intensity.floatValue(), 1));
//		}
//
//		String id = jmzSpec.getId();
//		spec.setID(id);
//		
//		int msLevel = jmzSpec.getMsLevel();
//		spec.setMsLevel(msLevel);
//
//		// should set scanNum, specIndex, activationMethod
//		
//		// sort peaks by increasing order of m/z
//		Collections.sort(spec);
//		
//		return spec;
//		
////		return null;
//	}
//	
//}