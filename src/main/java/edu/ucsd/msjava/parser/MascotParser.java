package edu.ucsd.msjava.parser;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedList;

import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Modification;
import edu.ucsd.msjava.msutil.Peptide;


public class MascotParser {
	
	private MascotParser() {}

	
	public static PSMList<PSM> parseFromDat(String fileName, boolean isDecoy)
	{
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		
		String s;
		int mode = 0;	// 0: parameters, 1: header 2: peptides, 3: queries, 6: decoy_peptides
		String specFileName = null;
//		String enzymeName = null;
//		String[] fixMod = null;
//		String[] variableMods = null; 
		int numQueries = 0;
		class SimplePSMList extends PSMList<PSM>
		{
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;
		}

		PSMList<PSM> psmList = new PSMList<PSM>();
		SimplePSMList[] psmArr = null;
		
		class SimpleMod {
			public SimpleMod(String name, float mass, String residue) {
				this.name = name;
				this.mass = mass;
				this.residue = residue;
			}
			String name;
			float mass;
			String residue;
		}
		
		SimpleMod[] mods = null;
		
		Hashtable<Integer, String> queryTitleMap = new Hashtable<Integer, String>();
		Hashtable<Integer, Integer> queryChargeMap = new Hashtable<Integer, Integer>();
		try {
			while((s = in.readLine()) != null)
			{
				if(s.startsWith("Content-Type:"))
				{
					if(s.equalsIgnoreCase("Content-Type: application/x-Mascot; name=\"parameters\""))
					{
						mode = 0;
						while((s=in.readLine()) != null) {
							if(s.startsWith("FILE="))
							{
								String[] token = s.split("[/\\\\]");
								if(specFileName == null)
									specFileName = token[token.length-1];
							}
							else if(s.startsWith("CLE="))
							{
								
							}
							else if(s.startsWith("MODS="))
							{
								
							}
							else if(s.startsWith("IT_MODS"))
							{
								int numMods = s.split(",").length;
								mods = new SimpleMod[numMods];
							}
							else if(s.equalsIgnoreCase("--gc0p4Jq0M2Yt08jU534c0p"))
								break;
						} 
					}
					else if(s.equalsIgnoreCase("Content-Type: application/x-Mascot; name=\"masses\""))
					{
						mode = 5;
					}
					else if(s.equalsIgnoreCase("Content-Type: application/x-Mascot; name=\"header\""))
					{
						mode = 1;
						// ex: queries=11905
						do {
							s = in.readLine();
						} while(!s.startsWith("queries"));
						numQueries = Integer.parseInt(s.substring(s.indexOf('=')+1));
					}
					else if(s.equalsIgnoreCase("Content-Type: application/x-Mascot; name=\"peptides\""))
					{
						if(!isDecoy)
						{
							mode = 2;
							psmArr = new SimplePSMList[numQueries];
						}
						else
							mode = -1;
					}
					else if(s.equalsIgnoreCase("Content-Type: application/x-Mascot; name=\"decoy_peptides\""))
					{
						if(isDecoy)
						{
							mode = 2;
							psmArr = new SimplePSMList[numQueries];
						}
						else
							mode = -1;
					}
					else if(s.startsWith("Content-Type: application/x-Mascot; name=\"query"))
					{
						mode = 3;
						int queryNum = Integer.parseInt(s.substring(s.lastIndexOf("query")+5, s.lastIndexOf('"')));
						do {
							s = in.readLine();
						} while(s != null && !s.startsWith("title="));
						String titleStr = s.substring(s.lastIndexOf('=')+1);
						String title = URLDecoder.decode(titleStr, "UTF-8");
						queryTitleMap.put(queryNum, title);
						do {
							s = in.readLine();
						} while(s != null && !s.startsWith("charge="));
						// charge=+n
						int charge;
						if(s.endsWith("Mr") || s.contains(","))
							charge = 0;
						else
						{
							if(s.contains("=+"))	// charge=+n
								charge = Integer.parseInt(s.substring(s.lastIndexOf('+')+1));
							else	// charge=n+
								charge = Integer.parseInt(s.substring(s.indexOf('=')+1, s.lastIndexOf('+')));
						}
						queryChargeMap.put(queryNum, charge);
					}
					else
						mode = -1;
				}
				if(mode == -1)
					continue;
				else if(mode == 5)	// masses
				{
					String keyword = "delta";
					if(s.length() > 0)
					{
						if(s.startsWith(keyword))
						{
							int modIndex = Integer.parseInt(s.substring(keyword.length(), s.indexOf('='))) - 1;
							String ptmStr = s.substring(s.indexOf('=')+1);
							String[] token = ptmStr.split(",");
							float mass = Float.parseFloat(token[0]);
							String[] token2 = token[1].split("\\s+");
							String ptmName = token2[0];
							StringBuffer residueStrBuf = new StringBuffer();
							for(int i=1; i<token2.length; i++)
								residueStrBuf.append(token2[i]);
							String residueStr = residueStrBuf.toString();
							String ptmResidues = residueStr.substring(residueStr.indexOf('(')+1, residueStr.lastIndexOf(')'));
							mods[modIndex] = new SimpleMod(ptmName, mass, ptmResidues);
						}
					}
				}
				else if(mode == 2)	// peptides or decoy_peptides
				{
					if(s.length() > 0)
					{
						String[] token = s.split("[=,;]");
						if(token.length > 3)	
						{
							// ex: q1_p1=1,798.496307,-0.001229,3,LKDLLR,6,20000000,11.61,0000002000000000000,0,0;"gi|60682641":0:323:328:3
							// miscleavages, peptide Mr, delta, #ions matched, peptide, peaks used from ions1,
							// variable modifications string, ions score, ion series found,
							// peaks used from Ions2, peaks used from Ions3,
							// accession string:data for second protein,frame number,start,end,multiplicity
							if(s.matches("q\\d+_p\\d+=.*"))
							{
//								float peptideMr = Float.parseFloat(token[2]);
								int queryNum = Integer.parseInt(token[0].substring(token[0].indexOf('q')+1, token[0].indexOf('_')));
//								int pepNum = Integer.parseInt(token[0].substring(token[0].lastIndexOf('p')+1));
								String peptide = token[5];
								String modVector = token[7];
								float score = Float.parseFloat(token[8]);
								String protein = token[12];
								do {
									s = in.readLine();
								} while (s != null && !s.matches("q\\d+_p\\d+_terms=.*"));
								// ex: q7_p7_terms=K,A
								if(s != null)
								{
									String[] token2 = s.split("[=,]");
									char nTermResidue = token2[1].charAt(0);
									char cTermResidue = token2[2].charAt(0);
									StringBuffer pepStr = new StringBuffer();
									
									// N-term modification
									if(modVector.charAt(0) != '0')
									{
										char v = modVector.charAt(0);
										int modIndex = v - '0' - 1;	// 0 based index
										float modMass = mods[modIndex].mass;
										if(modMass > 0)
											pepStr.append("+");
										pepStr.append(modMass);
									}
									for(int i=1; i<modVector.length()-1; i++)
									{
										char v = modVector.charAt(i);
										pepStr.append(peptide.charAt(i-1));
										if(v != '0')	// modified residue
										{
											int modIndex = v - '0' - 1;	// 0 based index
											float modMass = mods[modIndex].mass;
											if(modMass > 0)
												pepStr.append("+");
											pepStr.append(modMass);
										}
									}
									// C-term modification
									if(modVector.charAt(modVector.length()-1) != '0')
									{
										char v = modVector.charAt(modVector.length()-1);
										int modIndex = v - '0' - 1;	// 0 based index
										float modMass = mods[modIndex].mass;
										if(modMass > 0)
											pepStr.append("+");
										pepStr.append(modMass);
									}
									
									Peptide pep = new Peptide(pepStr.toString());
									if(!pep.isInvalid())
									{
										PSM psm = new PSM().peptide(pep).precedingResidue(nTermResidue).succeedingResidue(cTermResidue);
										psm.rawScore(score);
										psm.ptm(modVector);
										if(isDecoy)
											protein = "XXX" + protein;
										psm.protein(protein);
										if(psmArr[queryNum-1] == null)
											psmArr[queryNum-1] = new SimplePSMList();
										psmArr[queryNum-1].add(psm);
//										psm.score("QueryNum", queryNum);
									}
								}
							}
						}
					}
				}
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		}

		for(int queryNum=1; queryNum<=numQueries; queryNum++)
		{
			String title = queryTitleMap.get(queryNum);
			Integer charge = queryChargeMap.get(queryNum);
			if(title == null || charge == null || psmArr[queryNum-1] == null)
				continue;
			for(PSM psm : psmArr[queryNum-1])
			{
				psm.title(title);
				psm.charge(charge);
				psmList.add(psm);
			}
		}
		
		return psmList;
	}	
	
	public static HashMap<Integer,Integer> getQueryNumChargeMap(String fileName)
	{
		HashMap<Integer,Integer> map = new HashMap<Integer,Integer>();
		
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		String s;
		boolean isSummary = false;
		while((s=in.readLine()) != null)
		{
			if(!isSummary)
			{
				if(s.startsWith("--gc0p4Jq0M2Yt08jU534c0p"))
					isSummary = true;
				continue;
			}
			if(s.startsWith("qexp"))
			{
				int qNum = Integer.parseInt(s.substring(4, s.lastIndexOf('=')));
				int charge = Integer.parseInt(s.substring(s.lastIndexOf(',')+1, s.lastIndexOf('+')));
				map.put(qNum, charge);
			}
		}
		try {
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return map;
	}
}