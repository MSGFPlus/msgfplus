package misc;

import java.io.FileNotFoundException;
import java.util.HashSet;
import java.util.Random;

import parser.BufferedLineReader;
import parser.LineReader;
import parser.MgfSpectrumParser;

import msgf.Histogram;
import msgf.Profile;
import msscorer.NewScorerFactory;
import msutil.*;
import msgf.*;

public class SwedCAD {
	public static void main(String argv[]) throws Exception
	{
		makeAnnotatedMgf();
//		getBasicStats();
//		for(int length=7; length<=20; length++)
//			runMSGF(length);
		
		/*
		String specFileName = System.getProperty("user.home")+"/Research/Data/SwedCAD/swedCADAnnotated.mgf";
		String scoreParamFileName = System.getProperty("user.home")+"/Research/Data/RankScore/ionstat_FT_20ppm.txt";
		
		runMSProfile(specFileName, scoreParamFileName, new Tolerance(20, true), -1);
		*/
	}

	public static void simulateError() throws Exception
	{
		Tolerance[] fragTolerance = {
//				new Tolerance(0.5f, false),
//				new Tolerance(0.1f, false),
				new Tolerance(50, true),
				new Tolerance(30, true),
				new Tolerance(10, true),
		};
		String[] paramStr = {
//				"05",
//				"01",
				"50ppm",
				"30ppm",
				"10ppm",
		};
		for(int i=0; i<fragTolerance.length; i++)
		{
			Tolerance tol = fragTolerance[i];
			String specFileName = System.getProperty("user.home")
				+"/Research/Data/SwedCAD/swedCAD_"
				+tol.toString()
				+".mgf";
			String scoreParamFileName = System.getProperty("user.home")
				+"/Research/Data/SwedCAD/ionstat_FT_"
				+paramStr[i]
				+".txt";
			System.out.println("#"+tol);
//			runMSProfile(specFileName, scoreParamFileName, tol, -1);
		}		
	}
	
	
	public static void introduceError(Tolerance fragTolerance, Tolerance pmTolerance) throws Exception
	{
		String specFileName = System.getProperty("user.home")+"/Research/Data/SwedCAD/swedCADAnnotated.mgf";
		SpectraIterator itr = new SpectraIterator(specFileName, new MgfSpectrumParser());

		String outputFileName = System.getProperty("user.home")+"/Research/Data/SwedCAD/swedCAD_"+
			fragTolerance.toString();
		if(pmTolerance != null)
			outputFileName += "_"+pmTolerance.toString();
		outputFileName += ".mgf";
		
		SpectraContainer container = new SpectraContainer();

		Random rand = new Random();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			Peptide annotation = spec.getAnnotation();
			float theoPM = annotation.getParentMass();
			if(pmTolerance != null)
			{
				float newPM = getMassWithError(theoPM, pmTolerance, rand);
				float precursorMz = newPM/spec.getCharge()+(float)Composition.NEUTRON;
				spec.getPrecursorPeak().setMz(precursorMz);
			}
			
			for(Peak p : spec)
			{
				float mz = p.getMz();
				float newMz = getMassWithError(mz, fragTolerance, rand);
				p.setMz(newMz);
			}
			container.add(spec);
		}
		container.outputMgfFile(outputFileName);
		
		System.out.println(outputFileName);
	}
	
	private static float getMassWithError(float mass, Tolerance tol, Random rand)
	{
		double coeff;
		do {
			coeff = rand.nextGaussian();
		} while(coeff > 1 || coeff < -1);
		float newMass = mass+(float)coeff*tol.getToleranceAsDa(mass);
		return newMass;
		/*
		int sign = 1;
		if(rand.nextBoolean())
			sign = -1;
		float newMass = mass+sign*tol.getToleranceAsDa(mass)*rand.nextFloat();
		return newMass;
		*/
	}
	

	
	public static void getBasicStats()
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
		
		String specFileName = System.getProperty("user.home")+"/Research/Data/SwedCAD/swedCADAnnotated.mgf";
		SpectraIterator itr = null;
		try {
			itr = new SpectraIterator(specFileName, new MgfSpectrumParser());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		Histogram<Integer> hist = new Histogram<Integer>();
		Histogram<Integer> lengthHist = new Histogram<Integer>();
		HashSet<String> pepSet = new HashSet<String>();
		int numTryptic = 0;
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			Peptide annotation = spec.getAnnotation();
			if(annotation.hasTrypticCTerm())
				numTryptic++;
//			else
//				System.out.println(annotation);
			pepSet.add(annotation.toString());
			lengthHist.add(annotation.size());
			float errorPPM = (spec.getParentMass()-annotation.getParentMass())/annotation.getParentMass()*1e6f;
			int intErrorPPM = Math.round(errorPPM*10);
//			if(intErrorPPM == 3 || intErrorPPM == -3)
//				System.out.println(annotation+"\t"+spec.getParentMass()+"\t"+annotation.getParentMass()+"\t"+errorPPM);
			hist.add(intErrorPPM);
		}
		System.out.println("Error");
		hist.printSorted();
		
		System.out.println("\nLength");
		lengthHist.printSorted();
		
		System.out.println("\nNumPeptides: " + pepSet.size());
		System.out.println("\nNumTrypticPeptides: " + numTryptic);
	}
	
	public static void makeAnnotatedMgf()
	{
		String mgfOutputFileName = System.getProperty("user.home")+"/Research/Data/SwedCAD/swedCADAnnotated.mgf";
		mgfOutputFileName = System.getProperty("user.home")+"/Research/Data/SwedECD/swedECDAnnotated.mgf";
		
		SpectraContainer container = new SpectraContainer();
		
		String fileName = System.getProperty("user.home")+"/Research/Data/SwedCAD/all.txt";
		fileName = System.getProperty("user.home")+"/Research/Data/SwedECD/ECD.txt";
		LineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
		String s;
		Spectrum spec = null;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\\s+");
			if(token.length == 3)
			{
				String annotationStr = token[0];
				float parentMass = Float.parseFloat(token[1]);
				int charge = Integer.parseInt(token[2]);
				spec = new Spectrum(parentMass/charge+(float)Composition.NEUTRON, charge, 0);
				Peptide annotation = new Peptide(annotationStr, aaSet);
				spec.setAnnotation(annotation);
			}
			else if(token.length == 2)
			{
				assert(spec != null);
				spec.add(new Peak(Float.parseFloat(token[0]), Float.parseFloat(token[1]), 1));
			}
			else if(s.equalsIgnoreCase("<"))
			{
				assert(spec != null);
				container.add(spec);
			}
		}
		container.outputMgfFile(mgfOutputFileName);
		System.out.println(container.size()+" spectra are converted.");
	}
}
