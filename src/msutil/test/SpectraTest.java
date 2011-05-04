package msutil.test;

import msutil.Spectra;


/**
 * Test class for the object of this name
 * @author jung
 *
 */
public class SpectraTest {

  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    String mzXMLFile = userHome+"/Data/Cyclic/Javier/cyclicPeptides/MTVVI-20_1.medium.mzXML";
    new Spectra(mzXMLFile, 2);
  
  }
  
}
