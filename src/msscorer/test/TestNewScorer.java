package msscorer.test;

import java.util.HashSet;

import parser.*;
import msgf.Tolerance;
import msscorer.IonProbability;
import msutil.*;

public class TestNewScorer {
	public static void main(String argv[]) throws Exception
	{
//		compareNewAndOldScorer();
//		rescoreMSGFResult(new File("/home/sangtaekim/Research/Data/ISBETD/MSGFDB0312"), new File("/home/sangtaekim/Research/Data/ISBETD/MSGFDBRescored"));
//		mgfReadTest();
		ionSelectionTest();
	}
	
	public static void ionSelectionTest() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/Zubarev/AnnotatedSpectra/Zubarev_HCD_Annotated.mgf";
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		HashSet<String> pepSet = new HashSet<String>();
		
		SpectraContainer container = new SpectraContainer();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			if(pepSet.contains(spec.getAnnotationStr()))
				continue;
			else
				pepSet.add(spec.getAnnotationStr());
			if(spec.getCharge() == 2 && spec.getParentMass() < 962.4145f)
				container.add(spec);
		}
		
		IonType[] ions = IonType.getAllKnownIonTypes(4, true).toArray(new IonType[0]);
		Tolerance tol = new Tolerance(0.5f);
		IonProbability probGen = new IonProbability(container.iterator(), ions, tol).filter(new WindowFilter(6,50));
		probGen.segment(0, 2);
		System.out.println(container.size());
		float[] ionProb = probGen.getIonProb();
		int index = 0;
		for(IonType ion : ions)
		{
			if(ionProb[index] > 0.15f)
				System.out.println(ion.getName()+"\t"+ion.getOffset()+"\t"+ion.getCharge()+"\t"+ionProb[index]);
			index++;
		}
	}
	
	public static void mgfReadTest() throws Exception
	{
		String specFileName = System.getProperty("user.home")+"/Research/Data/ISBETD/AnnotatedSpectra/annotatedISBETD_ipro09.mgf";
		SpectraIterator itr = new SpectraIterator(specFileName, new MgfSpectrumParser());
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			System.out.println(spec.getScanNum());
		}
	}

}
