package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;


public class Clauser {
	public static void main(String argv[]) throws Exception
	{
		makeMGFWithActivationMethod();
	}
	
	public static void makeMGFWithActivationMethod() throws Exception
	{
		File mzXMLDir = new File("/home/sangtaekim/Test/cid_hcd_etd_triples/enc32");
		File mgfDir = new File("/home/sangtaekim/Test/cid_hcd_etd_triples/mgf");
		for(File f : mzXMLDir.listFiles())
		{
			String fileName = f.getName();
			if(!fileName.endsWith(".mzXML"))
				continue;
			
			System.out.println(fileName);
			String mgfFileName = fileName.substring(0, fileName.lastIndexOf('.'))+".mgf";
			File mgfFile = new File(mgfDir, mgfFileName);
			PrintStream out = null; 
			try {
				out = new PrintStream(new BufferedOutputStream(new FileOutputStream(mgfFile)));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			int numSpecs = 0;
			MzXMLSpectraIterator itr = new MzXMLSpectraIterator(f.getPath());
			while(itr.hasNext())
			{
				Spectrum spec = itr.next();
				spec.setTitle(f.getName()+":"+spec.getScanNum());
				if(spec.getActivationMethod() == null)
					spec.setActivationMethod(ActivationMethod.HCD);
				spec.outputMgf(out);
				numSpecs++;
			}
			out.close();
			System.out.println(numSpecs + " file converted.");
		}
	}
}
