package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.TSVResultParser;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

public class FindPSMIntersection {
    public static void main(String argv[]) throws Exception {
        find2();
    }

    public static void find() throws Exception {
        File ltqFile = new File("D:\\Research\\Data\\TrainingMSGFPlus\\AnnotatedSpectra\\HCD_HighRes_Tryp.mgf");
        File tofFile = new File("D:\\Research\\Data\\TrainingMSGFPlus\\AnnotatedSpectra\\CID_TOF_Tryp.mgf");

        HashSet<String> annotationSet = new HashSet<String>();

        SpectraIterator itr = new SpectraIterator(ltqFile.getPath(), new MgfSpectrumParser());
        while (itr.hasNext()) {
            Spectrum spec = itr.next();
            annotationSet.add(spec.getAnnotationStr());
        }

        itr = new SpectraIterator(tofFile.getPath(), new MgfSpectrumParser());
        while (itr.hasNext()) {
            String annotation = itr.next().getAnnotationStr();
            if (annotationSet.contains(annotation)) {
                System.out.println(annotation);
            }
        }
    }

    public static void find2() throws Exception {
        File trapFile = new File("D:\\Research\\Data\\TrapAndTOF\\TNT_biosensor_05_14Dec12_Phoenix_12-11-06_msgfplus.tsv");
        TSVResultParser parser = new TSVResultParser(trapFile);
        parser.parse(0.01f);

        Set<String> trapPepSet = parser.getPepSet();

        File tofFile = new File("D:\\Research\\Data\\TrapAndTOF\\CID_TOF_Tryp.mgf");
        File tofOutFile = new File("D:\\Research\\Data\\TrapAndTOF\\TOF.mgf");
        PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(tofOutFile)));

//		HashSet<String> annotationSet = new HashSet<String>();
        SpectraIterator itr = new SpectraIterator(tofFile.getPath(), new MgfSpectrumParser());
        while (itr.hasNext()) {
            Spectrum spec = itr.next();
            if (trapPepSet.contains(spec.getAnnotationStr())) {
                System.out.println(spec.getAnnotationStr());
                spec.outputMgf(out);
            }
//			annotationSet.add(spec.getAnnotationStr());
        }

        out.close();
        System.out.println("Done");
    }

}
