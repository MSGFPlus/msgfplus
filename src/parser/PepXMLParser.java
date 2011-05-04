package parser;

import java.io.IOException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Peptide;

import org.w3c.dom.*;
import org.xml.sax.SAXException;

public class PepXMLParser {
	public static PSMList<PSM> parse(String fileName)
	{
		PSMList<PSM> psmList = new PSMList<PSM>();
		
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = null;
		try {
			db = dbf.newDocumentBuilder();
		} catch (ParserConfigurationException e) {
			e.printStackTrace();
		}
		Document doc = null;
		try {
			doc = db.parse(fileName);
		} catch (SAXException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		Element docEle = doc.getDocumentElement();
		NodeList runSummaryList = docEle.getElementsByTagName("msms_run_summary");
		if(runSummaryList == null)
			return null;
		
		for(int i=0; i<runSummaryList.getLength(); i++)
		{
			Element rs = (Element)runSummaryList.item(i);
			String specExt = rs.getAttribute("raw_data"); 
			
			NodeList searchSummaryList = rs.getElementsByTagName("search_summary");
			assert(searchSummaryList.getLength() == 1);
			Element searchSummary = (Element)searchSummaryList.item(0);
			NodeList aaModList = searchSummary.getElementsByTagName("aminoacid_modification");
			AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
			if(aaModList != null && aaModList.getLength() > 0)
			{
				for(int j=0; j<aaModList.getLength(); j++)
				{
					Element aaMod = (Element)aaModList.item(j);
					char residue = aaMod.getAttribute("aminoacid").charAt(0);
//					float massDiff = Float.parseFloat(aaMod.getAttribute("massdiff"));
					boolean isFixed = aaMod.getAttribute("variable").equalsIgnoreCase("N");
					String name = aaMod.getAttribute("description");
					if(residue == 'C' && isFixed && name.contains("Carbamidomethyl"))
						aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
				}
			}
			//TODO: it is not rigorous
			
			NodeList specQueryList = rs.getElementsByTagName("spectrum_query");
			if(specQueryList == null)
				continue;
			for(int j=0; j<specQueryList.getLength(); j++)
			{
				Element sq = (Element)specQueryList.item(j);
				String spectrum = sq.getAttribute("spectrum");
				String specFileName = spectrum.substring(0, spectrum.indexOf('.')) + specExt;
				int startScan = Integer.parseInt(sq.getAttribute("start_scan"));
//				int endScan = Integer.parseInt(sq.getAttribute("end_scan"));
				float precursorNeutralMass = (float)Double.parseDouble(sq.getAttribute("precursor_neutral_mass"));
				int charge = Integer.parseInt(sq.getAttribute("assumed_charge"));
				float precursorMz = (precursorNeutralMass+charge*(float)Composition.PROTON)/charge;
				NodeList searchResultList = sq.getElementsByTagName("search_result");
				if(searchResultList == null)
					continue;
				for(int k=0; k<searchResultList.getLength(); k++)
				{
					Element sr = (Element)searchResultList.item(k);
					NodeList searchHitList = sr.getElementsByTagName("search_hit");
					if(searchHitList == null)
						continue;
					for(int l=0; l<searchHitList.getLength(); l++)
					{
						Element sh = (Element)searchHitList.item(l);
						
						String peptideStr = sh.getAttribute("peptide");
						String prevAA = sh.getAttribute("peptide_prev_aa");
						String nextAA = sh.getAttribute("peptide_next_aa");
						String protein = sh.getAttribute("protein");
						
						Peptide peptide = new Peptide(peptideStr);
						if(peptide.contains(null))
							continue;
						PSM psm = new PSM().peptide(new Peptide(peptideStr)).precedingResidue(prevAA.charAt(0)).succeedingResidue(nextAA.charAt(0))
							.charge(charge).specFileName(specFileName).scanNum(startScan).precursorMz(precursorMz).protein(protein);
						
						NodeList modInfoList = sh.getElementsByTagName("modification_info");
						if(modInfoList != null && modInfoList.getLength() > 0)
						{
							assert(modInfoList.getLength() == 1);
							Element modInfo = (Element)modInfoList.item(0);
							String modPeptide = modInfo.getAttribute("modified_peptide");
							StringBuffer peptideBuf = new StringBuffer();
							for(int index=0; index<modPeptide.length(); index++)
							{
								char c = modPeptide.charAt(index);
								if(c == '[')
								{
									char modResidue = modPeptide.charAt(index-1); 
									int modMass = 0;
									int sign = 1;
									assert(index+1<modPeptide.length());
									int nextChar = modPeptide.charAt(index+1); 
									if(nextChar == '+')
									{
										index++;
									}
									else if(nextChar == '-')
									{
										sign = -1;
										index++;
									}
									while((c=modPeptide.charAt(++index)) != ']')
									{
										assert(Character.isDigit(c));
										modMass = 10*modMass + (c-'0');
									}
									modMass *= sign;
									modMass -= aaSet.getAminoAcid(modResidue).getNominalMass();
									if(modMass >= 0)
										peptideBuf.append("+");
									peptideBuf.append(modMass);
								}
								else
									peptideBuf.append(c);
							}
							psm.peptide(new Peptide(peptideBuf.toString()));
							
							NodeList modList = modInfo.getElementsByTagName("mod_aminoacid_mass");
							if(modList != null)
							{
								for(int m=0; m<modList.getLength(); m++)
								{
									Element mod = (Element)modList.item(m);
									int position = Integer.parseInt(mod.getAttribute("position"))-1;	// starting from zero
									float mass = Float.parseFloat(mod.getAttribute("mass"));
									psm.ptm(position+":"+mass);
								}
							}
						}

						NodeList searchScoreList = sh.getElementsByTagName("search_score");
						if(searchScoreList != null)
						{
							for(int m=0; m<searchScoreList.getLength(); m++)
							{
								Element ss = (Element)searchScoreList.item(m);
								psm.score(ss.getAttribute("name"), Float.parseFloat(ss.getAttribute("value")));
							}
						}
						
						NodeList analysisResultList = sh.getElementsByTagName("analysis_result");
						if(analysisResultList != null)
						{
							for(int m=0; m<analysisResultList.getLength(); m++)
							{
								Element ar = (Element)analysisResultList.item(m);
								String analysisName = ar.getAttribute("analysis");
								NodeList analysisList = ar.getElementsByTagName(analysisName+"_result");
								assert(analysisList != null && analysisList.getLength() == 1);
								Element pr = (Element)analysisList.item(0);
								float prob = Float.parseFloat(pr.getAttribute("probability"));
								psm.score(analysisName, prob);
								if(analysisName.equalsIgnoreCase("interprophet"))
									psm.probScore(1-prob);
							}
						}
						psmList.add(psm);
					}
				}
			}
		}
		return psmList;
	}
	
	public static void main(String argv[]) throws Exception
	{
		String xmlFileName = "/home/sangtaekim/Research/Data/ISBETD/BioRep2TechRep1_ETD/OMScp_YeastCombNR_20070207_ForwDecoy/interact-ipro.pep.xml";
		PSMList<PSM> psmList = parse(xmlFileName);
		System.out.println(psmList.size());
	}
}
