package ui;

import java.io.File;

public class ScoringParamGen {

	public static void main(String argv[])
	{
		
	}
	
	public static void printUsageAndExit(String message)
	{
		if(message != null)
			System.out.println("Error: " + message + "\n");
		System.out.println("MSGFDB v"+ MSGFDB.VERSION + " (" + MSGFDB.RELEASE_DATE + ")");
		System.out.print("Usage: java -Xmx2000M -cp MSGFDB.jar ScoringParamGen\n"
				+ "********** Disco"
				+ "\t-s SpectrumFile (*.mzXML, *.mzML, *.mgf, *.ms2, *.pkl or *_dta.txt)\n" //, *.mgf, *.pkl, *.ms2)\n"
				+ "\t-d Database (*.fasta or *.fa)\n"
				+ "\t-t ParentMassTolerance (e.g. 2.5Da, 30ppm or 0.5Da,2.5Da)\n"
				+ "\t   Use comma to set asymmetric values. E.g. \"-t 0.5Da,2.5Da\" will set 0.5Da to the left (expMass<theoMass) and 2.5Da to the right (expMass>theoMass).\n"
				+ "\t[-o outputFileName] (Default: stdout)\n"
				+ "\t[-thread NumOfThreads] (Number of concurrent threads to be executed, Default: Number of available cores)\n"
				+ "\t[-tda 0/1] (0: don't search decoy database (default), 1: search decoy database to compute FDR)\n"
				+ "\t[-m FragmentationMethodID] (0: as written in the spectrum or CID if no info (Default), 1: CID, 2: ETD, 3: HCD, 4: Merge spectra from the same precursor)\n"
				+ "\t[-inst InstrumentID] (0: Low-res LCQ/LTQ (Default for CID and ETD), 1: TOF , 2: High-res LTQ (Default for HCD))\n"
				+ "\t[-e EnzymeID] (0: No enzyme, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: endogenous peptides)\n"
				+ "\t[-c13 0/1/2] (Number of allowed C13, Default: 1)\n"
				+ "\t[-nnet 0/1/2] (Number of allowed non-enzymatic termini, Default: 1)\n"
				+ "\t[-mod ModificationFileName] (Modification file, Default: standard amino acids with fixed C+57)\n"
				+ "\t[-minLength MinPepLength] (Minimum peptide length to consider, Default: 6)\n"
				+ "\t[-maxLength MaxPepLength] (Maximum peptide length to consider, Default: 40)\n"
				+ "\t[-minCharge MinPrecursorCharge] (Minimum precursor charge to consider if charges are not specified in the spectrum file, Default: 2)\n"
				+ "\t[-maxCharge MaxPrecursorCharge] (Maximum precursor charge to consider if charges are not specified in the spectrum file, Default: 3)\n"
				+ "\t[-n NumMatchesPerSpec] (Number of matches per spectrum to be reported, Default: 1)\n"
				+ "\t[-uniformAAProb 0/1] (0: use amino acid probabilities computed from the input database (Default), 1: use probability 0.05 for all amino acids)\n"
				);
		System.out.println("Example (high-precision): java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 30ppm -c13 1 -nnet 0 -tda 1 -o testMSGFDB.tsv");
		System.out.println("Example (low-precision): java -Xmx2000M -jar MSGFDB.jar -s test.mzXML -d IPI_human_3.79.fasta -t 0.5Da,2.5Da -nnet 0 -tda 1 -o testMSGFDB.tsv");
		System.exit(-1);
	}	
}
