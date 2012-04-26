package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import msutil.DBFileFormat;
import msutil.SpecFileFormat;
import params.FileParameter;
import params.IntParameter;
import params.ParamManager;
import parser.BufferedLineReader;

public class ConvertFastaForNAcetyl {
	public static final int VERSION = 7611;
	public static final String DATE = "04/24/2012";
	
	public static void main(String argv[]) throws Exception
	{
		ParamManager paramManager = new ParamManager("ConvertFastaForNAcetyl", String.valueOf(VERSION), DATE,
			"java -Xmx2000M -cp MSGFDB.jar misc.ConvertFastaForNAcetyl");

		FileParameter inputDBParam = new FileParameter("i", "DatabaseFile", "Input fasta file (*.fa, *.fasta)");
		inputDBParam.addFileFormat(DBFileFormat.FASTA);
		inputDBParam.fileMustExist();
		inputDBParam.mustBeAFile();
		paramManager.addParameter(inputDBParam);

		FileParameter outputDBParam = new FileParameter("o", "DatabaseFile", "Input fasta file (*.fa, *.fasta)");
		outputDBParam.addFileFormat(DBFileFormat.FASTA);
		outputDBParam.fileMustNotExist();
		paramManager.addParameter(outputDBParam);
		
		IntParameter prefixLengthParam = new IntParameter("maxLength", "MaxPepLength", "Maximum peptide length to consider, Default: 40");
		prefixLengthParam.minValue(1);
		prefixLengthParam.defaultValue(40);
//		addParameter(prefixLengthParam);
		
		if(argv.length == 0)
		{
			paramManager.printUsageInfo();
			return;
		}
		
		// Parse parameters
		String errMessage = paramManager.parseParams(argv); 
		if(errMessage != null)
		{
			System.err.println("[Error] " + errMessage);
			System.out.println();
			paramManager.printUsageInfo();
			return;
		}
		
		// Run program
		long time = System.currentTimeMillis();
		
		paramManager.printToolInfo();
		String errorMessage = convert(paramManager);
		if(errorMessage != null)
		{
			System.err.println("[Error] " + errorMessage);
			System.out.println();
		}
		else
			System.out.format("AnnotatedMgfToMSGFInput complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis()-time)/(float)1000);
		
	}

	public static String convert(ParamManager paramManager) throws Exception
	{
		File specFile = paramManager.getFile("i");
		File outputFile = paramManager.getFile("o");
		
//		SpectraIterator itr = new SpectraIterator(specFile.getPath(), new MgfSpectrumParser());
		BufferedLineReader in = new BufferedLineReader(specFile.getPath());
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		out.println("#SpectrumFile\tSpecIndex\tAnnotation\tCharge");
		String s;
		int specIndex = 0;
		int charge = 0;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("BEGIN"))
				specIndex++;
			else if(s.startsWith("END"))
			{
				charge = 0;
			}
			else if(s.startsWith("CHARGE="))
			{
				charge = Integer.parseInt(s.substring(s.indexOf('=')+1, s.lastIndexOf('+')));
			}
			else if(s.startsWith("SEQ="))
			{
				String pepSeq = s.substring(4);
				out.println(specFile.getName()+"\t"+specIndex+"\t"+"."+pepSeq+"."+"\t"+charge);
			}
		}

		in.close();
		out.close();
		
		return null;
	}	
}
