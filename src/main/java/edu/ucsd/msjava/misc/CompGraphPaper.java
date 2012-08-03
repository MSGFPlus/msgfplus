package edu.ucsd.msjava.misc;

import edu.ucsd.msjava.msgf.*;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.MgfSpectrumParser;

public class CompGraphPaper {
	
	public static void main(String argv[]) throws Exception
	{
//		processErrorSimulResults("5ppm");
//		sizeCompGraph();
		sizeAllPeptides();
	}

	public static void sizeAllPeptides() throws Exception
	{
		float num = 0;
		for(int i=1; i<=20; i++)
		{
			num += Math.pow(20, i);
			System.out.println(i+"\t"+num);
		}
	}
	
	public static void sizeCompGraph() throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		for(int i=1; i<=20; i++)
		{
			CompositionFactory factory = new CompositionFactory(aaSet, null, i);
			System.out.println(i+"\t"+factory.size());
		}
	}
	
	public static void processErrorSimulResults(String tolStr) throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/SwedCAD/diffErr.txt";
		fileName = System.getProperty("user.home")+"/Research/Data/SwedCAD/swedFedMSGF.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		
		float[] specProbThreshold = {1e-8f, 1e-9f, 1e-10f, 1e-11f, 1e-12f, 1e-13f};
		int[] numIdentifiedSpec = new int[specProbThreshold.length];
		
		boolean process = false;
		int numProcessedSpecs = 0;
		int numCorrectDeNovo = 0;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
			{
				if(s.equalsIgnoreCase("#"+tolStr))
					process = true;
				else if(process == true)
					break;
			}
			if(process == false)
				continue;
			String[] token = s.split("\t");
			if(token.length != 8)
				continue;
			String peptide = token[1];
			if(peptide.length() < 7 || peptide.length() > 20)
				continue;
			numProcessedSpecs++;
			
			int msgfScore = Integer.parseInt(token[3]);
			int rawScore = Integer.parseInt(token[4]);
			if(msgfScore == rawScore)
				numCorrectDeNovo++;
			float specProb = Float.parseFloat(token[5]);
			for(int i=0; i<specProbThreshold.length; i++)
			{
				if(specProb <= specProbThreshold[i])
					numIdentifiedSpec[i]++;
			}
		}
		System.out.println("NumSpectra\t" + numProcessedSpecs);
		System.out.println("DeNovo\t" + numCorrectDeNovo + "\t" + numCorrectDeNovo/(float)numProcessedSpecs);
		System.out.println("NumId");
		for(int i=0; i<specProbThreshold.length; i++)
		{
			System.out.println(specProbThreshold[i]+"\t"+numIdentifiedSpec[i]);
		}
	}
}
