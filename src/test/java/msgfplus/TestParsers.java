package msgfplus;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Iterator;

import org.junit.Assert;
import org.junit.Test;

import edu.ucsd.msjava.misc.ConvertToMgf;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;

public class TestParsers {

    @Test
    public void testReadingMgf() throws URISyntaxException {
        File mgfFile = new File(TestParsers.class.getClassLoader().getResource("test.mgf").toURI());
        SpectraAccessor specAcc = new SpectraAccessor(mgfFile);
        Iterator<Spectrum> itr = specAcc.getSpecItr();
        int numSpecs = 0;
        while(itr.hasNext()) {
            itr.next();
            numSpecs++;
        }
        Assert.assertTrue(numSpecs == 5760);
    }
    
    @Test
    public void testMzML() throws URISyntaxException, IOException {
        File mzMLFile = new File(TestParsers.class.getClassLoader().getResource("tiny.pwiz.mzML").toURI());
        File mgfFile = File.createTempFile("test", ".mgf");
        try {
            ConvertToMgf.convert(mzMLFile, mgfFile, false, null, null, -1, -1, -1, false); 
        } catch (Exception e) {
            e.printStackTrace();
        }
        mgfFile.deleteOnExit();
        
    }
    
    @Test
    public void testReadingMzXML() throws URISyntaxException {
        File mzXMLFile = new File(TestParsers.class.getClassLoader().getResource("testfile.mzXML").toURI());
        SpectraAccessor specAccessor = new SpectraAccessor(mzXMLFile);
        Iterator<Spectrum> specItr = specAccessor.getSpecItr();
        while(specItr.hasNext())
        {
            Spectrum spec = specItr.next();
            if(!spec.isCentroided())
            {
                System.out.println(spec.getScanNum() + " is not centroided.");
            }
        }
    }
    
    @Test
    public void testMzMLParser() throws URISyntaxException, IOException {
        File mzMLFile = new File(TestParsers.class.getClassLoader().getResource("tiny.pwiz.mzML").toURI());
        File mgfFile = File.createTempFile("test", ".mgf");
        try {
            ConvertToMgf.convert(mzMLFile, mgfFile, false, null, null, 43536, -1, -1, false);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
