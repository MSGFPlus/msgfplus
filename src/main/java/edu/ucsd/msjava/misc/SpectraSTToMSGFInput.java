package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.FileFormat;
import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.parser.BufferedLineReader;

public class SpectraSTToMSGFInput {
	public static final int VERSION = 7575;
	public static final String DATE = "04/13/2012";
	
	public static void main(String argv[]) throws Exception
	{
		ParamManager paramManager = new ParamManager("SpectraSTToMSGFInput", String.valueOf(VERSION), DATE,
			"java -Xmx2000M -cp MSGFDB.jar SpectraSTToMSGFInput");

		FileParameter stParam = new FileParameter("i", "SpectraSTResult", "SpectraST result file (*.txt)");
		stParam.addFileFormat(new FileFormat(".txt"));
		stParam.fileMustExist();
		stParam.mustBeAFile();
		paramManager.addParameter(stParam);
		
		FileParameter outputFileParam = new FileParameter("o", "OutputFile", "MS-GF input file");
		paramManager.addParameter(outputFileParam);

		paramManager.addSpecFileParam();
		
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
		File stFile = paramManager.getFile("i");
		File outputFile = paramManager.getFile("o");
		File specFile = paramManager.getFile("s");
		
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
		
		BufferedLineReader in = new BufferedLineReader(stFile.getPath());
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		out.println("#SpectrumFile\tSpecIndex\tAnnotation\tProtein\tCharge\tSpectraSTScore");
		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\\s+");
			if(token.length == 16)
			{
				String[] newToken = new String[17];
				for(int i=0; i<token.length; i++)
				{
					if(i < 2)
						newToken[i] = token[i];
					else if(i > 2)
						newToken[i+1] = token[i];
					else // i == 2
					{
						newToken[2] = token[2].substring(0, token[2].lastIndexOf('/')+2);
						newToken[3] = token[2].substring(token[2].lastIndexOf('/')+2);
					}
				}
				token = newToken;
			}
			if(token.length != 17)
				continue;
			String specIndex = token[0];
			
			String annotationStr = token[2].substring(0, token[2].lastIndexOf('/'));
			StringBuffer buf = new StringBuffer();
			int startIndex=0;
			if(annotationStr.startsWith("n[43]"))
			{
				buf.append("+42");
				startIndex = 5;
			}
			char prevAA = '\0';
			for(int i=startIndex; i<annotationStr.length(); i++)
			{
				char c = annotationStr.charAt(i);
				if(Character.isUpperCase(c))
					buf.append(c);
				else if(c == '[')
				{
					StringBuffer massBuf = new StringBuffer();
					while(annotationStr.charAt(++i) != ']')
						massBuf.append(annotationStr.charAt(i));
					int mass = Integer.parseInt(massBuf.toString());
					int residueMass = aaSet.getAminoAcid(prevAA).getNominalMass();
					int delMass = mass-residueMass;
					if(delMass > 0)
						buf.append("+");
					buf.append(delMass);
				}
				prevAA = c;
			}
			
			int charge = Integer.parseInt(token[2].substring(token[2].lastIndexOf('/')+1));
			float spectraSTScore = Float.parseFloat(token[11]);
			String protein = token[16];
			
			out.println(specFile.getName()+"\t"+specIndex+"\t."+buf.toString()+".\t"+protein+"\t"+charge+"\t"+spectraSTScore);
		}

		in.close();
		out.close();
		
		return null;
	}	
}
