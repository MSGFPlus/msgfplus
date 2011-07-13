package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

import msgf.Histogram;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.CompositionFactory;
import msutil.IonType;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;

import parser.BufferedLineReader;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraMap;

public class HCDCIDETD {
	
	public static void main(String argv[]) throws Exception
	{
//		generateAnnotatedSpectra();
//		analyzeFragmentErrors();
//		printParentMassErrorDist("CID");
//		compositionMSGF();
//		aminoAcidMSGF();
//		analyzeMSGFResults(false);
//		findGoodAAConst();
//		approxCompGraphTest();
	}

	static class ExtComposition 
	{
		int number;
//		public ExtComposition(int c, int h, int n, int o, int s) {
//			// MSB 6 (C) 13 (H) 5 (N) 5 (O) 3 (S)
//			number = c*0x04000000 + h*0x00002000 + n*0x00000100 + o*0x00000008 + s; 
//		}

		public ExtComposition(int c, int h, int n, int o, int s) 
		{
			int newC = c % 43;
			int newH = h + c/43*512 + n/19*264 + o/31*492 + s/7*222;
			int newN = n % 19;
			int newO = o % 31;
			int newS = s % 7;
			number = newC*0x04000000 + newH*0x00002000 + newN*0x00000100 + newO*0x00000008 + newS; 
		}
		
		public int getC() {	return (number & 0xFC000000) >>> 26; }
		public int getH() { return (number & 0x03FFE000) >> 13; }
		public int getN() { return (number & 0x00001F00) >> 8; } 
		public int getO() { return (number & 0x000000F8) >> 3; }
		public int getS() { return (number & 0x00000007); }
		public int getNumber()
		{
			return number;
		}
		
		@Override
		public int hashCode()
		{
			return number;
		}
		
		@Override
		public boolean equals(Object obj)
		{
			if(obj instanceof ExtComposition)
			{
				ExtComposition other = (ExtComposition)obj;
				if(number == other.number)
					return true;
			}
			return false;
		}
		
		@Override
		public String toString()
		{
			return new String(getC()+" "+getH()+" "+getN()+" "+getO()+" "+getS());
		}
		
		public float getMass()
		{
			return (float)(getC()*Composition.C+getH()*Composition.H+getN()*Composition.N+getO()*Composition.O+getS()*Composition.S);
		}
	}
	
	public static void approxCompGraphTest() throws Exception
	{
		int maxLength = 30;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		CompositionFactory allCompositions = new CompositionFactory(aaSet, null, maxLength);

		HashSet<ExtComposition> approxCompSet = new HashSet<ExtComposition>();
		System.out.println("Composition building done.");
		
		float maxErr = 0;
		Composition worstComp = Composition.NIL;
		for(int n : allCompositions.getData())
		{
			Composition c = new Composition(n);
			int numC = c.getC();
			int numH = c.getH();
			int numN = c.getN();
			int numO = c.getO();
			int numS = c.getS();
			
			ExtComposition newComp = new ExtComposition(numC, numH, numN, numO, numS);
			approxCompSet.add(newComp);
			float err = (newComp.getMass()-c.getMass());
			if(Math.abs(err) > maxErr)
			{
				maxErr = Math.abs(err);
				worstComp = c;
			}
//			System.out.println(c+"\t"+err);
		}
		System.out.println("MaxErr\t"+worstComp.toString()+"\t"+worstComp.getMass()+"\t"+maxErr);
		System.out.println("CompositionSize\t"+allCompositions.size());
		System.out.println("Size\t" + approxCompSet.size());
		
	}
	
	public static void findGoodAAConst() throws Exception
	{
//		System.out.println(Composition.C*Constants.INTEGER_MASS_SCALER-Math.round(Composition.C));
//		System.out.println(Composition.H*Constants.INTEGER_MASS_SCALER-Math.round(Composition.H));
//		System.out.println(Composition.N*Constants.INTEGER_MASS_SCALER-Math.round(Composition.N));
//		System.out.println(Composition.O*Constants.INTEGER_MASS_SCALER-Math.round(Composition.O));
//		System.out.println(Composition.S*Constants.INTEGER_MASS_SCALER-Math.round(Composition.S));
		
		double target = Composition.S;
		double base = Math.abs(target-Composition.H*Math.round(target/Composition.H));
//		double base = Math.abs(target*Constants.INTEGER_MASS_SCALER-Math.round(target));
		
		ArrayList<Integer> bestI = new ArrayList<Integer>();
		ArrayList<Integer> bestJ = new ArrayList<Integer>();
		ArrayList<Double> diffList = new ArrayList<Double>();
		for(int i=1; i<100; i++)
		{
			for(int j=1; j<1000; j++)
			{
				double diff = Math.abs(target*i - Composition.H*j);
				if(diff < base)
				{
					bestI.add(i);
					bestJ.add(j);
					diffList.add(target*i - Composition.H*j);
				}
			}
		}
		System.out.println("Base: " + base);
		for(int i=0; i<bestI.size(); i++)
			System.out.println(bestI.get(i)+"\t"+bestJ.get(i)+"\t"+diffList.get(i));
	}
	
	public static void analyzeMSGFResults(boolean target) throws Exception
	{
		String fileName = "/home/sangtaekim/Developments/MS_Java/bin/PNNLCompGraph.txt";
//		String fileName = "/home/sangtaekim/Developments/MS_Java/bin/PNNLAAGraph.txt";
		String s;
		BufferedLineReader in = new BufferedLineReader(fileName);
//		while(!(s=in.readLine()).endsWith("54220"))
//			continue;
		in.readLine();
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length != 20)
				continue;
			String annotation = token[10];
			String pepStr = annotation.substring(annotation.indexOf('.')+1,annotation.lastIndexOf('.'));
			if(pepStr.length() > 20)
				continue;
			String prot = token[8];
			if(target && prot.startsWith("Reversed"))
				continue;
			if(!target && !prot.startsWith("Reversed"))
				continue;
			if(token[19].startsWith("Parent"))
				continue;
			float specProb = Float.parseFloat(token[19]);
			float xcorr = Float.parseFloat(token[5]);
			System.out.println(specProb+"\t"+xcorr);
		}
	}
	
	private static AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
	public static void analyzeFragmentErrors() throws Exception
	{
		String[] methods = {"CID", "ETD", "HCD"};
		for(String method : methods)
			analyzeFragmentErrors(method);
		
	}
	
	public static void analyzeFragmentErrors(String method) throws Exception
	{
		IonType[] ionTypes;
		File dir = new File(System.getProperty("user.home")+"/Research/Data/CIDETDHCD");
		File specFile = new File(dir.getPath()+File.separator+"AnnotatedSpectra/Annotated_"+method+".mgf");
		
		double prefixMass = Composition.H;
		double suffixMass = Composition.H2O;
		
		SpectraIterator itr = new SpectraIterator(specFile.getPath(), new MgfSpectrumParser().aaSet(aaSet));
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			Peptide pep = spec.getAnnotation();
			for(int i=0; i<pep.size()-1; i++)
			{
				prefixMass += pep.get(i).getAccurateMass();
				suffixMass += pep.get(pep.size()-1-i).getAccurateMass();
			}
			
		}
	}
	
	public static void generateAnnotatedSpectra() throws Exception
	{
		String[] methods = {"CID", "ETD", "HCD"};
		for(String method : methods)
			generateAnnotatedSpectra(method);
	}
	
	public static void generateAnnotatedSpectra(String method) throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/CIDETDHCD");
		File resultFile = new File(dir.getPath()+File.separator+"AnalysisResults/"+method+"_FDR1pct.txt");
		File specFile = new File(dir.getPath()+File.separator+"HCD_CID_ETD_EBCP_LMW_1y_50ug_300rt.mzXML");
		File outputFile = new File(dir.getPath()+File.separator+"AnnotatedSpectra/Annotated_"+method+".mgf");

		MzXMLSpectraMap map = new MzXMLSpectraMap(specFile.getPath());
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

		Histogram<Integer> pmErrorHist = new Histogram<Integer>();
		BufferedLineReader in = new BufferedLineReader(resultFile.getPath());
		
		in.readLine();	// header
		String s;
		int numSpecs = 0;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length != 25)
				continue;
			
			String db = token[9];
			if(db.startsWith("Reverse"))
				continue;
			int scanNum = Integer.parseInt(token[2]);
			Spectrum spec = map.getSpectrumBySpecIndex(scanNum);
			int charge = Integer.parseInt(token[4]);
			spec.setCharge(charge);
			float expMass = spec.getParentMass();
			String annotation = token[11];
			Peptide pep = aaSet.getPeptide(annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.')));
			
			float theoMass = pep.getParentMass();
			float parentMassError = expMass - theoMass;
			float parentMassErrorPPM = parentMassError*1e6f/theoMass;
			if(parentMassErrorPPM > 50 || parentMassErrorPPM < -50)
			{
				System.out.println(annotation + " " + expMass + " " + theoMass + " " + parentMassErrorPPM);
			}
			else if(parentMassErrorPPM < 20 && parentMassErrorPPM > -20)
			{
				spec.setAnnotation(pep);
				spec.outputMgf(out);
				numSpecs++;
			}
				
			pmErrorHist.add(Math.round(parentMassErrorPPM));
		}
		out.close();
		System.out.println("Done: " + numSpecs + " spectra.");
		pmErrorHist.printSorted();
		
	}
	
	public static void printParentMassErrorDist(String method) throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/CIDETDHCD");
		File resultFile = new File(dir.getPath()+File.separator+"AnalysisResults/"+method+"_FDR1pct.txt");
		File specFile = new File(dir.getPath()+File.separator+"HCD_CID_ETD_EBCP_LMW_1y_50ug_300rt.mzXML");
		File outputFile = new File(dir.getPath()+File.separator+"AnnotatedSpectra/Annotated_"+method+".mgf");

		MzXMLSpectraMap map = new MzXMLSpectraMap(specFile.getPath());

		Histogram<Integer> pmErrorHistPPM = new Histogram<Integer>();
		Histogram<Integer> pmErrorHistDa = new Histogram<Integer>();
		BufferedLineReader in = new BufferedLineReader(resultFile.getPath());
		
		in.readLine();	// header
		String s;
		int numSpecs = 0;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length != 25)
				continue;
			
			String db = token[9];
			if(db.startsWith("Reverse"))
				continue;
			int scanNum = Integer.parseInt(token[2]);
			Spectrum spec = map.getSpectrumBySpecIndex(scanNum);
			int charge = Integer.parseInt(token[4]);
			spec.setCharge(charge);
			float expMass = spec.getParentMass();
			String annotation = token[11];
			Peptide pep = aaSet.getPeptide(annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.')));
			
			float theoMass = pep.getParentMass();
			float parentMassError = expMass - theoMass;
			float parentMassErrorPPM = parentMassError*1e6f/theoMass;
			if(parentMassErrorPPM > 50 || parentMassErrorPPM < -50)
			{
				System.out.println(annotation + " " + expMass + " " + theoMass + " " + parentMassErrorPPM);
			}
			numSpecs++;
				
			pmErrorHistPPM.add(Math.round(parentMassErrorPPM));
			pmErrorHistDa.add(Math.round(parentMassError*1000));
			
		}
		System.out.println("PPM");
		pmErrorHistPPM.printSorted();
		
		System.out.println("\tDa");
		pmErrorHistDa.printSorted();
	}		
		
	
}
