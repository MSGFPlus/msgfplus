package edu.ucsd.msjava.msdbsearch;

import java.io.File;


import edu.ucsd.msjava.msutil.SpectraContainer;
import edu.ucsd.msjava.msutil.SpectraMap;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MgfSpectrumParser;

/**
 * Take CID mgf and ETD mgf and generate a single mgf with CID/ETD pairs
 * @author sangtaekim
 *
 */
public class MakePairedSpectra {

	public static void main(String argv[])
	{
		if(argv.length != 3)
			printErrorAndExit("Illegal parameters!");
		
		File cidFile = new File(argv[0]);
		if(!cidFile.isFile())
			printErrorAndExit(argv[0] + " is not a file.");
		String ext = cidFile.getName().substring(cidFile.getName().lastIndexOf('.')+1);
		if(!ext.equalsIgnoreCase("mgf"))
			printErrorAndExit(argv[0] + " must be a mgf file.");
			
		File etdFile = new File(argv[1]);
		if(!etdFile.isFile())
			printErrorAndExit(argv[1] + " is not a file.");
		ext = etdFile.getName().substring(etdFile.getName().lastIndexOf('.')+1);
		if(!ext.equalsIgnoreCase("mgf"))
			printErrorAndExit(argv[1] + " must be a mgf file.");
		
		File outFile = new File(argv[2]);
		
		merge(cidFile, etdFile, outFile);
	}
	
	public static void printErrorAndExit(String message)
	{
		System.err.println(message);
		System.out.println("usage: java -Xmx3500M MakePairedSpectra CIDSpectra (*.mgf) ETDSpectra (*.mgf) OutputSpectra (*.mgf)");
		System.exit(-1);
	}

	public static void merge(File cidFile, File etdFile, File outFile)
	{
		SpectraMap cidMap = new SpectraMap(cidFile.getPath(), new MgfSpectrumParser());
		SpectraMap etdMap = new SpectraMap(etdFile.getPath(), new MgfSpectrumParser());
		SpectraContainer mergedContainer = new SpectraContainer();
		
		for(Integer cidScanNum : cidMap.getSpecIndexList())
		{
			int etdScanNum = cidScanNum+1;
			Spectrum etdSpec = etdMap.getSpectrumBySpecIndex(etdScanNum);
			if(etdSpec != null)
			{
				Spectrum cidSpec = cidMap.getSpectrumBySpecIndex(cidScanNum);
				float cidPreMz = cidSpec.getPrecursorPeak().getMz();
				float etdPreMz = etdSpec.getPrecursorPeak().getMz();
				float error = Math.abs(cidPreMz-etdPreMz);
				if(error < 0.001f && cidSpec.getCharge() == etdSpec.getCharge())
				{
					mergedContainer.add(cidSpec);
					mergedContainer.add(etdSpec);
				}
			}
		}
		mergedContainer.outputMgfFile(outFile.getPath());
	}
}
