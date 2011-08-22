package parser;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import msutil.ActivationMethod;
import msutil.Spectrum;

public class MzXMLToMgfConverter {
	public static void main(String argv[])
	{
		if(argv.length < 2 || argv.length % 2 != 0)
			printUsageAndExit("Illegal parameters");
		
		File mzXMLFile = null;
		File mgfFile = null;
		int minCharge = 0;
		int maxCharge = Integer.MAX_VALUE;
		ActivationMethod activationMethod = null;
		int minMSLevel = 2;
		int maxMSLevel = 2;
		
		for(int i=0; i<argv.length; i+=2)
		{
			if(argv[i].equalsIgnoreCase("-s"))
			{
				mzXMLFile = new File(argv[i+1]);
				if(!mzXMLFile.exists())
				{
					printUsageAndExit(argv[i+1] + " doesn't exist.");
				}
			}
			else if(argv[i].equalsIgnoreCase("-t"))
			{
				mgfFile = new File(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-c"))
			{
				String chargeStr = argv[i+1];
				String minChargeStr = "";
				String maxChargeStr = "";
				if(chargeStr.contains(":"))	// range
				{
					String[] chargeToken = chargeStr.split(":");
					if(chargeToken.length > 1)
						minChargeStr = chargeToken[0].trim();
					if(chargeToken.length == 2)
						maxChargeStr = chargeToken[1].trim();
				}
				else
					minChargeStr = maxChargeStr = chargeStr;
				try {
					if(minChargeStr.length() > 0)
						minCharge = Integer.parseInt(minChargeStr);
					if(maxChargeStr.length() > 0)
						maxCharge = Integer.parseInt(maxChargeStr);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal charge: " + argv[i+1]);
				} 
			}
			else if(argv[i].equalsIgnoreCase("-m"))
			{
				// 0: All, 1: CID, 2: ETD, 3: HCD
				if(argv[i+1].equalsIgnoreCase("1"))
					activationMethod = ActivationMethod.CID;
				else if(argv[i+1].equalsIgnoreCase("2"))
					activationMethod = ActivationMethod.ETD;
				else if(argv[i+1].equalsIgnoreCase("3"))
					activationMethod = ActivationMethod.HCD;
			}
			else if(argv[i].equalsIgnoreCase("-l"))
			{
				String msLevelStr = argv[i+1];
				String minMSLevelStr = "";
				String maxMSLevelStr = "";
				if(msLevelStr.contains(":"))	// range
				{
					String[] token = msLevelStr.split(":");
					if(token.length > 1)
						minMSLevelStr = token[0].trim();
					if(token.length == 2)
						maxMSLevelStr = token[1].trim();
				}
				else
					minMSLevelStr = maxMSLevelStr = msLevelStr;
				try {
					if(minMSLevelStr.length() > 0)
						minMSLevel = Integer.parseInt(minMSLevelStr);
					if(maxMSLevelStr.length() > 0)
						maxMSLevel = Integer.parseInt(maxMSLevelStr);
					if(minMSLevel <= 0 || minMSLevel > maxMSLevel)
					{
						printUsageAndExit("Illigal MS Level: " + argv[i+1]);
					}
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal MS Level: " + argv[i+1]);
				} 
			}
		}
		if(mzXMLFile == null)
			printUsageAndExit("MzXML file is missing!");
		if(mgfFile == null)
			printUsageAndExit("Mgf file is missing!");
		
		convert(mzXMLFile, mgfFile, minCharge, maxCharge, activationMethod, minMSLevel, maxMSLevel);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.err.println(message);
		System.out.print("usage: java MzXMLToMgfConverter\n"
				+ "\t-s mzXMLFileName\n"
				+ "\t-t mgfFileName\n"
				+ "\t[-c charge (Default: all, e.g. 2, 2:4, :3, 2:)]\n"
				+ "\t[-m 0/1/2 (0: All (default), 1: CID, 2: ETD, 3: HCD)]\n"
				+ "\t[-l msLevel (Default: 2, e.g. 2, 2:4, :3, 2:)]\n"
		);
		System.exit(-1);
	}
	
	public static void convert(File mzXMLFile, File mgfFile, int minCharge, int maxCharge, ActivationMethod actMethod, int minMSLevel, int maxMSLevel)
	{
		int numFiles = 0;
		PrintStream out = null;
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(mgfFile)));
		} catch (IOException e) {
			e.printStackTrace();
		}
		MzXMLSpectraIterator itr = new MzXMLSpectraIterator(mzXMLFile.getPath(), minMSLevel, maxMSLevel+1);
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int charge = spec.getCharge();
			if(charge < minCharge || charge > maxCharge)
				continue;
			if(spec.getActivationMethod() != actMethod)
				continue;
			spec.setTitle(mzXMLFile.getName()+":"+spec.getScanNum()+":"+spec.getActivationMethod());
			spec.outputMgf(out);
			numFiles++;
		}
		out.close();
		System.out.println(numFiles + " spectra are converted.");
	}
}
