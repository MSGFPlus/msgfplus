package misc;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import parser.BufferedLineReader;
import parser.TSVParser;

import msdbsearch.CompactFastaSequence;
import msdbsearch.CompactSuffixArray;
import msscorer.NewRankScorer;
import msscorer.NewScorerFactory;
import msscorer.ScoringParameterGeneratorWithErrors;
import msscorer.NewScorerFactory.SpecDataType;
import msutil.ActivationMethod;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.InstrumentType;
import msutil.Protocol;
import msutil.SpectraContainer;

public class MSGFPlusPaper {
	public static void main(String argv[]) throws Exception
	{
//		nominalMassTable();
//		checkPeptidesWithNominalMassErrors();
		countTotalNumberOfPartitions();
//		aLPModel();
//		modInLib();
//		vennDiagram();
	}

	public static void vennDiagram() throws Exception
	{
		String dataset = "HCD_FT";
		File msgfdbFile = new File("/home/sangtaekim/Research/Data/Heck_DDDT/backup1229/"+dataset+"_7ppm.tsv");
		File mascotDATFile = new File("/home/sangtaekim/Research/Data/Heck_DDDT/Mascot/Mascot_"+dataset+".dat");
		File percolatorFile = new File("/home/sangtaekim/Research/Data/Heck_DDDT/Mascot/Percolator_"+dataset+".tsv");
		
		HashSet<Integer> msgfdbID = new HashSet<Integer>();
		
		BufferedLineReader in = new BufferedLineReader(msgfdbFile.getPath());
		
		int scanNumIndex = 2;
		int fdrIndex = 13;
		
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length < fdrIndex)
				continue;
			float fdr = Float.parseFloat(token[fdrIndex]);
			if(fdr <= 0.01f)
				msgfdbID.add(Integer.parseInt(token[scanNumIndex]));
		}		
		
		in.close();
		
		HashMap<Integer,Integer> queryNumToScanNum = new HashMap<Integer,Integer>();
		in = new BufferedLineReader(mascotDATFile.getPath());
		String keyWord = "Content-Type: application/x-Mascot; name=\"query";
		while((s=in.readLine()) != null)
		{
			if(s.startsWith(keyWord))
			{
				int queryNum = Integer.parseInt(s.substring(keyWord.length(), s.lastIndexOf('"')));
				in.readLine();
				in.readLine();
				s = in.readLine();
				assert(s.startsWith("scans="));
				int scanNum = Integer.parseInt(s.substring(s.lastIndexOf('=')+1));
				queryNumToScanNum.put(queryNum, scanNum);
			}
		}
		in.close();
		
		HashSet<Integer> percolatorID = new HashSet<Integer>();
		
		in = new BufferedLineReader(percolatorFile.getPath());
		
		int qValueIndex = 2;
		in.readLine();
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 5)
				continue;
			float qValue = Float.parseFloat(token[qValueIndex]);
			if(qValue <= 0.01f)
			{
				int queryNum = Integer.parseInt(token[0].substring(token[0].indexOf(':')+1, token[0].indexOf(';')));
				int scanNum = queryNumToScanNum.get(queryNum);
				percolatorID.add(scanNum);
			}
		}		
		
		int msgfdbOnly = 0;
		int percolatorOnly = 0;
		int both = 0;
		for(int id : msgfdbID)
		{
			if(percolatorID.contains(id))
				both++;
			else
				msgfdbOnly++;
		}
		for(int id : percolatorID)
		{
			if(!msgfdbID.contains(id))
				percolatorOnly++;
		}
		
		System.out.println("MS-GFDB all: " + msgfdbID.size());
		System.out.println("Percolator all: " + percolatorID.size());
		
		System.out.println("Both: " + both);
		System.out.println("MS-GFDB only: " + msgfdbOnly);
		System.out.println("Percolator only: " + percolatorOnly);
	}
	
	public static void modInLib() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/NISTLib/yeast_2011_05_24_it.pepidx";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		HashSet<String> modNames = new HashSet<String>();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#") || s.length() == 0)
				continue;

			String[] token = s.split("\\s+");
			String pepStr = token[0];
			String pepInfoStr = token[1];
			
			String[] tokenInfo = pepInfoStr.split("\\|");
			int charge = Integer.parseInt(tokenInfo[0]);
			
			String modInfo = tokenInfo[1];
			String[] tokenMod = modInfo.split("/");
			int numMods = Integer.parseInt(tokenMod[0]);
			for(int i=1; i<tokenMod.length; i++)
			{
				String[] mod = tokenMod[i].split(",");
				int location = Integer.parseInt(mod[0]);
				char residue = mod[1].charAt(0);
				String modName = mod[2];
				modNames.add(modName+":"+residue+":"+location);
			}			
		}
		
		for(String m : modNames)
			System.out.println(m);
	}
	
	public static void aLPModel() throws Exception
	{
		File specFile = new File("/home/sangtaekim/Research/Data/IonStat/SpectraForTraining/ETD_LowRes_aLP.mgf");
		SpecDataType dataType = new SpecDataType(ActivationMethod.ETD, InstrumentType.LOW_RESOLUTION_LTQ, Enzyme.ALP, Protocol.NOPROTOCOL);
		
//		File specFile = new File("/home/sangtaekim/Research/Data/IonStat/SpectraForTraining/ETD_LowRes_Tryp.mgf");
//		SpecDataType dataType = new SpecDataType(ActivationMethod.ETD, InstrumentType.LOW_RESOLUTION_LTQ, Enzyme.TRYPSIN, Protocol.NOPROTOCOL);
		
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		File outputDir = new File("/home/sangtaekim/Developments/MS_Java_Dev/bin");
		
		ScoringParameterGeneratorWithErrors.generateParameters(specFile, dataType, aaSet, outputDir, true, true, true);
	}
	
	public static void countTotalNumberOfPartitions() throws Exception
	{
		int numPartitions = 0;
		for(ActivationMethod method : ActivationMethod.getAllRegisteredActivationMethods())
		{
			if(method == ActivationMethod.FUSION || method == ActivationMethod.ASWRITTEN)
				continue;
			for(InstrumentType inst : InstrumentType.getAllRegisteredInstrumentTypes())
			{
				for(Enzyme enzyme : Enzyme.getAllRegisteredEnzymes())
				{
					for(Protocol protocol : Protocol.getAllRegisteredProtocols())
					{
						SpecDataType condition = new SpecDataType(method, inst, enzyme, protocol);
						InputStream is = ClassLoader.getSystemResourceAsStream("resources/ionstat/"+condition+".param");
						if(is != null)
						{
							System.out.println(condition);
							NewRankScorer scorer = new NewRankScorer(new BufferedInputStream(is));
							numPartitions += scorer.getParitionSet().size();
						}
					}
				}
			}
		}
		System.out.println("NumPartitions: " + numPartitions);
	}
	
	public static void nominalMassTable() throws Exception
	{
		String delimiter = " & ";
		System.out.println("Residue" + delimiter + "NominalMass" + delimiter + 
				"RealMass" + delimiter + "RescaledMass" + delimiter + 
				"Error(RealMass)" + delimiter + "Error(RescaledMass)"
				);
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		aaSet = AminoAcidSet.getStandardAminoAcidSet();
		float sumReal = 0;
		float sumRes = 0;
		float sumPPM = 0;
		float absSumReal = 0;
		float absSumRes = 0;
		for(AminoAcid aa : aaSet)
		{
			float mass = aa.getMass();
			int nominalMass = aa.getNominalMass();
			float rescaledMass = mass*0.9995f;
			float error = mass-nominalMass;
			float rescaledError = rescaledMass-nominalMass;
			float errorPPM = (rescaledMass-nominalMass)/mass*1e6f;
//			System.out.println(aa.getResidue()+delimiter+mass+delimiter+nominalMass+delimiter+rescaledMass+delimiter+error+delimiter+errorPPM+"\\\\");
			System.out.format("%c%s%d%s%.3f%s%.3f%s%.3f%s%.3f\\\\\n",
					aa.getResidue(),delimiter,
					nominalMass,delimiter,
					mass,delimiter,
					rescaledMass,delimiter,
					error,delimiter,
					rescaledError,delimiter
					);
			sumReal += error;
			sumRes += rescaledError;
			absSumReal += Math.abs(error);
			absSumRes += Math.abs(rescaledError);
			
			sumPPM += errorPPM;
		}
		System.out.println("AverageRealError\t"+sumReal/20);
		System.out.println("AverageRescaledError\t"+sumRes/20);
		System.out.println("AbsAverageRealError\t"+absSumReal/20);
		System.out.println("AbsAverageRescaledError\t"+absSumRes/20);
		System.out.println("AverageErrorPPM\t"+sumPPM/20);
	}
	
	public static void checkPeptidesWithNominalMassErrors() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/IPI/IPI_human_3.87.fasta";
		int maxPeptideLength = 100;
		CompactFastaSequence fastaSequence = new CompactFastaSequence(fileName);
		CompactSuffixArray sa = new CompactSuffixArray(fastaSequence, maxPeptideLength);
		sa.measureNominalMassError(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
	}
}
