package edu.ucsd.msjava.msutil.test;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;


/**
 * Test MzXMLSpectraIterator.
 * @author jung
 *
 */
public class MzXMLSpectraIteratorTest {

  public static void main(String argv[]) {
	  test2();
  }
  
  public static void test1()
  {
	    String fileName = System.getProperty("user.home")+"/Research/Data/HEK293/H293-total-try-a-200ug-2D34-081905-LTQ1-01.mzXML";
	    MzXMLSpectraIterator itr = new MzXMLSpectraIterator(fileName);
	    int specNum = 0;
	    while(itr.hasNext())
	    {
	      specNum++;
	      Spectrum spec = itr.next();
	      System.out.println(spec.getScanNum()+"\t"+spec.getPrecursorPeak().getMz() + "\t" + spec.getPrecursorPeak().getCharge()
	          +"\t"+spec.size());
	    }
	    System.out.println("NumSpecs: " + specNum);
  }
  
  public static void test2()
  {
	  String fileName = System.getProperty("user.home")+"/Research/Data/HeckWhole/090121_NM_Trypsin_22.mzXML";
	  MzXMLSpectraIterator itr = new MzXMLSpectraIterator(fileName);
	    int specNum = 0;
	    int prevScanNum = 0;
	    float prevPrecursorMz = 0;
	    ActivationMethod prevMethod = ActivationMethod.ETD;
	    while(itr.hasNext())
	    {
	      specNum++;
	      Spectrum spec = itr.next();
	      System.out.println(spec.getScanNum()+"\t"+spec.getPrecursorPeak().getMz() + "\t" + spec.getPrecursorPeak().getCharge()
	          +"\t"+spec.getActivationMethod());
	      if(spec.getActivationMethod() == ActivationMethod.ETD)
	      {
	    	  if(prevMethod != ActivationMethod.CID || 
	    			  spec.getScanNum() != prevScanNum+1 || spec.getPrecursorPeak().getMz() != prevPrecursorMz)
	    	  {
	    		  System.err.println("Wrong");
	    		  System.exit(-1);
	    	  }
	      }
	      else if(spec.getActivationMethod() == ActivationMethod.CID)
	      {
	    	  if(prevMethod != ActivationMethod.ETD)
	    	  {
	    		  System.err.println("Wrong");
	    		  System.exit(-1);
	    	  }
	      }
	      prevMethod = spec.getActivationMethod();
	      prevScanNum = spec.getScanNum();
	      prevPrecursorMz = spec.getPrecursorPeak().getMz();
	    }
	    System.out.println("NumSpecs: " + specNum);
  }
  
}
