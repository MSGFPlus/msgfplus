package trex.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import msutil.AminoAcid;
import msutil.AminoAcidSet;

import parser.BufferedLineReader;

public class ProteinSearch {
	static Hashtable<String, String> proteinTable = null;
	static Hashtable<String, String> contaminantTable = null;
	static Hashtable<String, String> peptideSet = null; // peptide w/o hydroxy proline, no Q no I.
	static Hashtable<Integer, HashSet<String>> specSet = null;
	static Hashtable<Integer, ArrayList<String>> specWithUniquePep = null;
	
	public static int numFoundTh = 0, numSpec=0, numPep=0, numPro1=0, numPro2=0;
	public static boolean writeSpec = false;
	public static boolean writeScanNum = false;
	public static boolean writePrecursor_Successor = true;
	public static boolean writeProtein = false;
	public static boolean writePeptide = false;
	public static boolean onlyCOL = false;
	public static boolean onlyNONCOL = false;
	public static boolean allowHydroxyProlineforNonCOL = true;
	
	//static int num_p_th = 100;
	static int num_KR_th = 100;
	
	final public static String[] dbName = {
		"/home/kwj/workspace/inputs/DataBases/ipi.CHICK.fasta",
		"/home/kwj/workspace/inputs/DataBases/DOG.fasta",
		"/home/kwj/workspace/inputs/DataBases/ipi.HUMAN.fasta",
		"/home/kwj/workspace/inputs/DataBases/ipi.BOVIN.fasta",
		"/home/kwj/workspace/inputs/DataBases/ipi.MOUSE.fasta",
		"/home/kwj/workspace/inputs/DataBases/ipi.RAT.fasta",
		
		"/home/kwj/workspace/inputs/uniprot_sprot.fasta",
		
		
	//	"/home/kwj/workspace/inputs/DataBases/panda.fasta",
	//	"/home/kwj/workspace/inputs/DataBases/reverse_uniprot_sprot.fasta",
	};
	
	public static String DataSet;
	public static String[] dirName;
	
	public static void setDataSet(String d){
		DataSet = d;
		String[] dn={
				"/home/kwj/workspace/outputs/"+ DataSet + "/dbresults/chick",
				"/home/kwj/workspace/outputs/"+ DataSet + "/dbresults/dog",
				"/home/kwj/workspace/outputs/"+ DataSet + "/dbresults/human",
				"/home/kwj/workspace/outputs/"+ DataSet + "/dbresults/bovine",
				"/home/kwj/workspace/outputs/"+ DataSet + "/dbresults/mouse",
				"/home/kwj/workspace/outputs/"+ DataSet + "/dbresults/rat",
				
				"/home/kwj/workspace/outputs/"+ DataSet + "/dbresults/uniprot",
		//		"/home/kwj/workspace/outputs/"+ DataSet + "/dbresults/panda",
				
		//		"/home/kwj/workspace/outputs/"+ DataSet + "/dbresults/reverse_uniprot",
			};
		
		dirName = dn;
		
	}
	
	
	
	
	
	/*
	 * final public static String[] dirName ={
		"/home/kwj/workspace/outputs/cavebear/dbresults/chick",
		"/home/kwj/workspace/outputs/cavebear/dbresults/dog",
		"/home/kwj/workspace/outputs/cavebear/dbresults/human",
		"/home/kwj/workspace/outputs/cavebear/dbresults/bovin",
		"/home/kwj/workspace/outputs/cavebear/dbresults/mouse",
		"/home/kwj/workspace/outputs/cavebear/dbresults/rat",
		"/home/kwj/workspace/outputs/cavebear/dbresults/uniprot"
	};*/
	
	static class outFilter implements FilenameFilter{
		String[] ext = { ".out"};
		
		public boolean accept(File dir, String name) {
			name = name.toLowerCase();
			return name.endsWith(ext[0]);// || name.endsWith(ext[1]);
		} 
	}

	public static void initDB() {proteinTable = null;}
	
	final static String[] keratinKeyWord = {
	 "ENSCAFT00000025346",
	 "ENSCAFT00000025319",
	 "ENSCAFT00000038760",
	 "ENSCAFT00000025381",
	 "ENSCAFT00000025254",
	 "ENSCAFT00000025317",
	 "ENSCAFT00000025250",
	 "ENSCAFT00000025270",
	 "ENSCAFT00000025351",
	 "ENSCAFT00000025410",
	 "ENSCAFT00000011645",
	 "ENSCAFT00000025255",
	 "ENSCAFT00000011630",
	 "ENSCAFT00000025265",
	 "ENSCAFT00000025272",
	 "ENSCAFT00000011636",
	 "ENSCAFT00000011623",
	 "ENSCAFT00000025391",
	 "ENSCAFT00000011567",
	 "KERATIN",
	 "=KRT"
	};
	
	public static void loadDB(String filename) throws IOException{
		BufferedLineReader in = null;
		in = new BufferedLineReader(filename);
		if(proteinTable == null) proteinTable =  new Hashtable<String, String>();
		if(contaminantTable == null) contaminantTable =  new Hashtable<String, String>();
		String s;
		String key = new String(), value = new String();
		
		while((s=in.readLine()) != null){
			
			if(s.startsWith(">")){
				if(!key.isEmpty()){
					boolean keratinfound = false;
					for(String keratin: keratinKeyWord){
						if(key.toUpperCase().contains(keratin)){
							keratinfound = true;
							break;
						}
					}
					if(!keratinfound) 
						proteinTable.put(key, value);
					else contaminantTable.put(key, value);
				}
				key = s;
				value = new String();
				if(onlyNONCOL && (key.toLowerCase().contains("=col") || key.toLowerCase().contains("collagen") || key.toLowerCase().contains("loc363458"))) key = ""; //
				if(onlyCOL && (!key.toLowerCase().contains("=col") && !key.toLowerCase().contains("collagen") && !key.toLowerCase().contains("loc363458"))) key = "";
			}
			else{
				value += s;
			}	
		}
		if(!key.isEmpty()) proteinTable.put(key, value);
		in.close();
	}
	
	
	public static void loadContaminantDB(String filename) throws IOException{
		BufferedLineReader in = null;
		in = new BufferedLineReader(filename);
		if(contaminantTable == null) contaminantTable =  new Hashtable<String, String>();
		String s;
		String key = new String(), value = new String();
		
		while((s=in.readLine()) != null){
			if(s.startsWith(">")){
				if(!key.isEmpty()) contaminantTable.put(key, value);
				key = s;
				value = new String();
			}
			else{
				value += s;//.replace('Q', 'K').replace('I', 'L');
			}	
		}
		if(!key.isEmpty()) contaminantTable.put(key, value);
		in.close();
	}
	
	public static void initSpecAndPeptides(){
		peptideSet = new Hashtable<String, String>();
		specSet = new Hashtable<Integer, HashSet<String>>();
	};
	
	public static void loadSpecs(String filename) throws IOException{ // should be fixed later!!!! Now each file should have different spec nums
		File f = new File(filename);
		File[] files = null;
		
		if(f.isDirectory())
			files = f.listFiles(new outFilter());
		else{
			files = new File[1];
			files[0] = f;
		}
		// peptide, spectrum with highest match score
		String s;
		//numSpec = 0;
		
		BufferedLineReader in = null;
		
		if(specWithUniquePep == null)
			specWithUniquePep = new Hashtable<Integer, ArrayList<String>>();
		
		for(File file : files){
			//System.out.println(file);
			in = new BufferedLineReader(file.toString());

			int prevScanNum = -1;
			
			while((s=in.readLine()) != null){
				if(s.startsWith("#") || s.startsWith("Spec Title")) continue;
				String t[] = s.split("\t");
				if(t.length < 5) continue;
				
				int scanNum = Integer.parseInt(t[1]);
				if(prevScanNum == scanNum){
					int prevscore = Integer.parseInt(specWithUniquePep.get(prevScanNum).get(0).split("\t")[5]);
					int score = Integer.parseInt(t[5]);
					
					ArrayList<String> sa = null;
					if(prevscore < score){ 
						sa = new ArrayList<String>();	
					}else if(prevscore == score){
						sa = specWithUniquePep.get(scanNum);
					}
					
					if(prevscore <= score){
						sa.add(s);
						specWithUniquePep.put(scanNum, sa);
					}
					// Integer.parseInt(t[5]);
				}else{
					ArrayList<String> sa = new ArrayList<String>();	 sa.add(s);
					specWithUniquePep.put(scanNum, sa);
				}
				
				
				prevScanNum = scanNum;
			}
			
			in.close();
		}
	}
	
	public static void loadPeptides() throws IOException{
	
		
		Hashtable<String, ArrayList<String>> peptideSetTmp = new Hashtable<String, ArrayList<String>>();
		// peptide, spectrum with highest match score
		String s;
			
		for(int n : specWithUniquePep.keySet()){
			ArrayList<String> sa = specWithUniquePep.get(n); 
			for(int i=0; i<sa.size(); i++){
				s = sa.get(i);
		//	}
		//	while((s=in.readLine()) != null){
				if(s.startsWith("#") || s.startsWith("Spec Title")) continue;
				
				String t[] = s.split("\t");
				if(t.length < 5) continue;
				
				int scanNum = Integer.parseInt(t[1]);
				
				HashSet<String> sp = null;
				if(specSet.containsKey(scanNum)){
					sp = specSet.get(scanNum);
				}else{
					sp = new HashSet<String>();
				}
				sp.add(t[4].substring(t[4].indexOf('.')+1, t[4].lastIndexOf('.')));
				specSet.put(scanNum, sp);
				
				String p = t[4].substring(t[4].indexOf('.')+1, t[4].lastIndexOf('.'));
			
				if(onlyNONCOL && !allowHydroxyProlineforNonCOL && p.contains("p")) continue;
				
				p = p.replace("p", "P");
			
				if(!peptideSetTmp.containsKey(p)){
					ArrayList<String> u = new ArrayList<String>(); u.add(s);
					peptideSetTmp.put(p, u);
				}else{
					ArrayList<String> u = peptideSetTmp.get(p); u.add(s);
					peptideSetTmp.put(p, u);
				}
			}
		}
		
		for(String k : peptideSetTmp.keySet()){
			ArrayList<String> v = peptideSetTmp.get(k);
			int maxscore = Integer.MIN_VALUE;
			String specWithHighestScore = new String();
			for(String spec : v){
				String t[] = spec.split("\t");
				int score =  Integer.parseInt(t[5]);
				maxscore = maxscore > score ? maxscore : score;
				if(score == maxscore) specWithHighestScore = spec;
			}
			peptideSet.put(k, specWithHighestScore);
		
			
		}
		
		//System.out.println(filename + "\n# spec: " + numSpec + "\n# pep with contaminant: " + peptideSet.size());
	}
	
	
				
	static void eraseContaminant(){
		if(contaminantTable != null){
			for(String pro : contaminantTable.keySet()){
				String proSeq = contaminantTable.get(pro);
				String proSeq_QI_replaced = proSeq.replace('Q', 'K').replace('I', 'L');
				//String s = new String();
	
				ArrayList<String> toErase = new ArrayList<String>();
				ArrayList<Integer> toEraseScanNum = new ArrayList<Integer> ();
				for(int scanNum : specSet.keySet()){
					for(String sp : specSet.get(scanNum)){
						if(proSeq_QI_replaced.contains(sp.toUpperCase())){
							toEraseScanNum.add(scanNum);
							break;
						}
					}
				}
				
				for(String pep : peptideSet.keySet()){
					if(proSeq_QI_replaced.contains(pep))	toErase.add(pep);
					if(toEraseScanNum.contains(Integer.parseInt(peptideSet.get(pep).split("\t")[1]))) toErase.add(pep);
				}
				
				for(int scanNum: toEraseScanNum)
					specSet.remove(scanNum);
				
				for(String k:toErase)
					peptideSet.remove(k);
			}
		}
		
	}
	
	
	
	public static void writeOutput() throws IOException{
		writeOutput("");
	}
	
	public static void writeOutput(String filename) throws IOException{
		boolean write = true;
		if(filename.isEmpty()) write = false;
		
		BufferedWriter out = null;
		if(write) out =  new BufferedWriter(new FileWriter(filename));
		Hashtable<Integer, ArrayList<String>> outString = new Hashtable<Integer, ArrayList<String>>();
		int maxNumFound = -1;
		
		//contaminantTable
		eraseContaminant();
		
		
		
		//System.out.println("# pep before matching: " + peptideSet.size());
		
		HashSet<String> matchedPep = new HashSet<String>();
		
		for(String pro : proteinTable.keySet()){
		//	HashSet<String> matchedPep_to_this_protein = new HashSet<String>();
			String proSeq = proteinTable.get(pro);
			String proSeq_QI_replaced = proSeq.replace('Q', 'K').replace('I', 'L');
			char[] proSeqCovered = new char[proSeq_QI_replaced.length()];
			ArrayList<Integer> coveredIndices = new ArrayList<Integer>();
			String s = new String();
			int numFound = 0;
			
			for(String pep : peptideSet.keySet()){
				boolean isFound = false;
				if(proSeq_QI_replaced.contains(pep)){
					
					ArrayList<String> pep_from_DBs = new ArrayList<String>();
					for(int i=0; i<proSeq_QI_replaced.length(); i++){ // not used replace to avoid overlapping peptides...
						if(proSeq_QI_replaced.startsWith(pep, i)){
							if(i>0 && !(proSeq_QI_replaced.charAt(i-1) == 'K' || proSeq_QI_replaced.charAt(i-1) == 'R') ) continue;
							String pep_from_DB = new String();
							
							int num_p =0, num_KR = 0;
							
							for(int j=i; j<i+pep.length(); j++){
								coveredIndices.add(j);
								pep_from_DB += proSeq.charAt(j);
								if(proSeq.charAt(j) == 'K' || proSeq.charAt(j) == 'R') num_KR ++;
								
							}
							//if(i-1>0);
							//if(i+pep.length() + 1 < proSeq.length())
							//pep_from_DB = proSeq.charAt(i-1) + "."+ pep_from_DB + "." + proSeq.charAt(i+pep.length());
							if(!(pep_from_DB.endsWith("K") || pep_from_DB.endsWith("R"))) continue;
							if(num_KR > num_KR_th) continue;
							isFound = true;
							
							if(writePeptide)
								s += "\t" + pep_from_DB + "\n";
							if(!pep_from_DBs.contains(pep_from_DB)) pep_from_DBs.add(pep_from_DB);
						}
						
					}
					
					// num_P, num_cleavage, ends_with_KR check!!
					//isCol
					
					
					if(isFound){ 
						
						if(writeSpec){
							for(String pep_from_DB : pep_from_DBs){
								//System.out.println(pep_from_DB);
								String t = peptideSet.get(pep);
								
								String u = t.split("\t")[4];
								char[] k = u.substring(u.indexOf('.')+1, u.lastIndexOf('.')).toCharArray();
								//System.out.println(new String(k));

								for(int i=0; i<k.length; i++){
									if(k[i] == 'p'){
										
									}else{
										k[i] = pep_from_DB.charAt(i);
									}
								}	
								if(!writeScanNum){
									if(writePrecursor_Successor){
										s += "\t" +
											u.substring(0, u.indexOf('.')+1) + new String(k) + u.substring(u.lastIndexOf('.'))
											+ "\t" + t.split("\t")[5] + "\t" + t.split("\t")[6] + "\n";
									}else
										s += "\t" +
										 new String(k) 
										+ "\t" + t.split("\t")[5] + "\t" + t.split("\t")[6] + "\n";
								}else{
									if(writePrecursor_Successor){
										s += "\t" + t.split("\t")[1] + "\t" +
											u.substring(0, u.indexOf('.')+1) + new String(k) + u.substring(u.lastIndexOf('.'))
											+ "\t" + t.split("\t")[5] + "\t" + t.split("\t")[6] + "\n";
									}else
										s += "\t" +  t.split("\t")[1] + "\t" +
										 new String(k) 
										+ "\t" + t.split("\t")[5] + "\t" + t.split("\t")[6] + "\n";
								}
							//	s += "\t" + t + "\n";
							}
						}
						
						
						numFound ++;
						matchedPep.add(pep);
					}
					
					
				}
				
			}
			
			
			if(numFound > numFoundTh){
				
				float numCovered = 0;
				for(int i=0; i<proSeq.length(); i++) {
					if(coveredIndices.contains(i)){
						proSeqCovered[i] = '*';
						numCovered++;
					}else
						proSeqCovered[i] = proSeq.charAt(i);
				}
				//new String(proSeqCovered);	
				
				float coverage = numCovered/proSeqCovered.length;				
				
				maxNumFound = maxNumFound >= numFound ? maxNumFound : numFound;
				

				if(writeProtein) s = pro + "\t# distinct matching peptides: " + numFound + "\tCoverage: " + coverage*100 +  "%\n" +
					"\t" + proSeq + "\n" + "\t"  + (new String(proSeqCovered)) + "\n" + s;
				else s = pro + "\t# distinct matching peptides: " + numFound + "\tCoverage" + coverage*100 +  "%\n" + s;
				
				if(!outString.containsKey(numFound)){
					ArrayList<String> t = new ArrayList<String>(); t.add(s);
					outString.put(numFound, t);
				}else{
					ArrayList<String> t = outString.get(numFound); t.add(s);
					outString.put(numFound, t);
				}
			}
		}
		
		//System.out.println("# pep : " + matchedPep.size());
		
		//if(numPep > 0 )System.out.println(matchedPep);
		
		
		
		
		
		ArrayList<String> toErase = new ArrayList<String>();
		
		for(String p: peptideSet.keySet()){
			if(!matchedPep.contains(p)) toErase.add(p);
		}
		
		Hashtable<String, String> peptideSetTmp = new Hashtable<String, String>(peptideSet);
		
		for(String k:toErase)
			peptideSetTmp.remove(k);
	//	out.write("# pep : " + matchedPep.size() + "\n");
		
		numPep = peptideSetTmp.size();
		
		HashSet<Integer> usedScanNum = new HashSet<Integer>();
		for(String k: peptideSetTmp.keySet()){
			for(int scanNum : specSet.keySet()){
				for(String spec : specSet.get(scanNum)){
					if(spec.toUpperCase().equals(k)){
						usedScanNum.add(scanNum);break;
					}
				}
			}
		}
		
		numSpec = usedScanNum.size();
		
		numPro1 = 0;
		numPro2 = 0;
		for(int i=maxNumFound ; i>numFoundTh; i--){
			if(outString.containsKey(i)){
				for(String s : outString.get(i)){
					if(write)out.write(s);
					numPro1++;
					if(i > 1) numPro2++;
				}
			}
		}
		
		//System.out.println("# protein: " + numPro +"\n# protein matched by 2 or more pep: " + numPro2);
		if(write)out.close();
	}
	
	/*
	static void rewritePeptides(String filename) throws IOException{
		File f = new File(filename);
		File[] files = null;
		
		BufferedWriter out = null;
		
		
		if(f.isDirectory())
			files = f.listFiles(new outFilter());
		else{
			files = new File[1];
			files[0] = f;
		}
		
		String s;

		BufferedLineReader in = null;
		for(File file : files){
			//System.out.println(file);
			in = new BufferedLineReader(file.toString());
			out =  new BufferedWriter(new FileWriter(file.toString() + "_revised"));
			
			while((s=in.readLine()) != null){
				if(s.startsWith("#") || s.startsWith("Spec Title")) continue;
				
				boolean specUsed = false;
				for(String pep : peptideSet.keySet()){
					if(s.contains(pep)){
						specUsed = true;
						break;
					}
				}
				
				if(specUsed) out.write(s + "\n");	//
			}
			
			in.close();
			out.close();
		}
	}
	
	*/
	static void applyMutation(String pro_file , String mut_file) throws IOException{
		
		BufferedLineReader in = new BufferedLineReader(pro_file);
		String s;
		String original_protein = null;
		String cov_protein = null;
		
		while((s=in.readLine()) != null){
			if(s.startsWith(">") || s.isEmpty()) continue;
			if(s.contains("*")) cov_protein = s.trim();
			else original_protein = s.trim();
			
		}
		
		
		
		System.out.println(cov_protein);
		
		//System.out.println(original_protein.length());
		float num_cov = 0;
		for(int i=0; i<cov_protein.length(); i++)
			if(cov_protein.charAt(i) == '*') num_cov++;
		System.out.println(num_cov / original_protein.length());
		
		in = new BufferedLineReader(mut_file);
		
		
		//Hashtable<Integer, Integer> mutated_pep = new Hashtable<Integer, Integer>();// mutated pep beginning position, mutated position
		int[] mutated_position = new int[original_protein.length()];
		int[] highest_score = new int[original_protein.length()];
		int[] end_position = new int[original_protein.length()];
		
		for(int i=0; i<mutated_position.length; i++){// init
			mutated_position[i] = -1;
			highest_score[i] = -1;
			end_position[i] = -1;
		}
		
		int current_position = -1;
		int current_end_position = -1;
		int current_mutated_position = -1;
		
		while((s=in.readLine()) != null){
			if(s.startsWith(">")){
				current_position = Integer.parseInt(s.substring(1).split(",")[0]);
				current_end_position = Integer.parseInt(s.substring(1).split(",")[1]);
				current_mutated_position = Integer.parseInt(s.substring(1).split(",")[2]);
			}else{
				int score = Integer.parseInt(s.split("\t")[2]);
				highest_score[current_position] = highest_score[current_position] > score? highest_score[current_position] : score;
				if(highest_score[current_position] == score){
					mutated_position[current_position] = current_mutated_position;
					end_position[current_position] = current_end_position;
				}
			}
		}
		
		char[] resulting_pro = cov_protein.toCharArray();
		
		for(int i=0; i<resulting_pro.length; i++){
			//if(resulting_pro[i] == '*') continue;
			if(end_position[i] > 0){
				for(int j=i; j<=end_position[i]; j++){
					resulting_pro[j] = '#';
				}
				resulting_pro[i+mutated_position[i]] = '@';
				i = end_position[i];
			}
		}
		//System.out.println(new String(resulting_pro));
		
		String protein = new String(resulting_pro);
		num_cov = 0;
		float num_mutated = 0, num_cov_by_mu = 0;
		char c = '\0';
		
		for(int i=0; i<protein.length(); i++){
		
			if(c == '*' && protein.charAt(i) != '*' || c == '#' && protein.charAt(i) != '#' || c == '@' && protein.charAt(i) != '@'){
				System.out.print("}");
			}
			
		/*	if(i==12||i==904||i==985||i==1058||i==1301||i==1453)
				System.out.print("\\textcolor{yellow}{");
			if(i==13||i==905||i==986||i==1059||i==1302||i==1454)
				System.out.print("}");
			*/	
			if(c != '*' &&  protein.charAt(i) == '*'){
				System.out.print("\\textcolor{red}{");
			}
			
			if(c != '#' &&  protein.charAt(i) == '#'){
				System.out.print("\\textcolor{blue}{");
			}
			
			if(c != '@' &&  protein.charAt(i) == '@'){
				System.out.print("\\textcolor{green}{");
			}
			if(i%80 == 0)System.out.print("\\\\");
			System.out.print(original_protein.charAt(i));
			
			c = protein.charAt(i);

			if(c == '*') num_cov++;
			if(c == '@') num_mutated++;
			if(c == '#' || c == '@') num_cov_by_mu++;
		}
		System.out.println();
		System.out.println(num_cov / original_protein.length());
		System.out.println(num_cov_by_mu / original_protein.length());
		System.out.println((num_cov + num_cov_by_mu) / original_protein.length());
		System.out.println(num_mutated + "\t" + original_protein.length());
		in.close();
	//System.out.println(new String(pro).length());
	}
	
	static void generateMutatedProteinDB(String infile, String outfile) throws IOException{
		ArrayList<String> proteins = new ArrayList<String>();
		BufferedLineReader in = new BufferedLineReader(infile);
		String s;
		
		while((s=in.readLine()) != null){
			if(s.startsWith(">")) continue;
			if(s.contains("*"))
				proteins.add(s.trim());
		}
		
		BufferedWriter out =new BufferedWriter(new FileWriter(outfile));
		
		for(String protein: proteins){
			ArrayList<String> subProList = new ArrayList<String>();
			ArrayList<String> subProIndexList = new ArrayList<String>();
			
			if(!protein.endsWith("*")) protein += "*";
			
			int beginIndex = 0;
			String subPro = new String();
			for(int i=0; i<protein.length(); i++){
				char aa = protein.charAt(i);
				//char aa_before = protein.charAt(i-1);
				if(aa == '*'){
					subPro = new String();
					beginIndex = i+1;
				}else{
					if(aa == 'K' || aa == 'R'){
						subPro += aa;
						subProList.add(subPro);
						subProIndexList.add(beginIndex+","+i);
						subPro = new String();
						beginIndex = i+1;
						//}
					}else{
						subPro += aa;
					}
				}
			}
			
			for(int i=0; i<subProList.size();i++){
				String t = subProList.get(i);
				for(int j=0; j<t.length();j++){
					char c[] = t.toCharArray();
					char original = c[j];
					
					for(AminoAcid aa : AminoAcidSet.getStandardAminoAcidSet()){
						c[j] = aa.getResidueStr();
						if(original != c[j]){
							out.write(">" + subProIndexList.get(i)+","+ j + "," +aa +"\n");
							out.write(new String(c) + "\n");
						}
					}
				}
			}
		}
		in.close();
		out.close();
	//	for(String s : subProIndexList)
		//	System.out.println(s);
		
	}
	
	static void printPep(String filename) throws IOException{
		BufferedLineReader in = new BufferedLineReader(filename);
		String s;
		Hashtable<String, HashSet<String>> pepTable = new Hashtable<String,  HashSet<String>>();
		String value = new String();
		while((s=in.readLine()) != null){
			//if(s.startsWith("&")) break;
			if(s.startsWith("#")){
				value = s.substring(1);
			}
			else{
				String[] keys = s.split(",");
				//System.out.println(keys.length);
				for(String key : keys){
					if(pepTable.containsKey(key)){
						HashSet<String> t = pepTable.get(key);
						t.add(value);
						pepTable.put(key, t);
					}
					else{
						HashSet<String> t = new HashSet<String>();
						t.add(value);
						pepTable.put(key, t);
					}
				}
			}
		}
		
		int h=0,d=0,hd=0;
		for(int i = 6; i>=1;i--){
			for(String key : pepTable.keySet()){
				HashSet<String> v = pepTable.get(key);
				
					if(v.size() != i) continue;
					System.out.println(key + "\t" + pepTable.get(key));
					if(v.contains("O")) h++;
					if(v.contains("D")) d++;
					
					if(v.contains("O") && v.contains("D")) hd++;
				}
			}
		System.out.println(h+"\t"+d +"\t" + hd);
		in.close();
	}
	
	public static void main(String[] argv) throws IOException{
		
		ProteinSearch.setDataSet("Had");
		//ProteinSearch.setDataSet("cavebear");
		ProteinSearch.loadContaminantDB("/home/kwj/workspace/inputs/CommonContaminants.fasta");
		
		ProteinSearch.writeSpec = true;
		ProteinSearch.initDB();ProteinSearch.initSpecAndPeptides();
		ProteinSearch.loadDB("/home/kwj/workspace/inputs/DataBases/ipi.HUMAN.fasta");
		ProteinSearch.loadSpecs("/home/kwj/workspace/outputs/"+ProteinSearch.DataSet+"/dbresults/human/");
		ProteinSearch.loadPeptides();
		ProteinSearch.writeOutput("/home/kwj/workspace/outputs/"+ProteinSearch.DataSet+"/dbresults/human/spec");
		
		
		//generateMutatedProteinDB("/home/kwj/workspace/outputs/cavebear/dbresults/tmp", "/home/kwj/workspace/inputs/pandacol.fasta");
		
		
		/*
		 initDB();initPeptides();
		loadDB("/home/kwj/workspace/inputs/elephant_mutated1.fasta");
		//loadDB("/home/kwj/workspace/inputs/DataBases/ipi.HUMAN.fasta");
		loadPeptides(dir + "elephant_mutated1.out");
		writeOutput(dir  + "elephant_mutated1_spec");
		
		applyMutation(
				"SPPCLIESSEQGPTGPPGRDGEDGIPGPPGPPGPPGPPGLGGNFAAQYDAKGIGLGPGPMGLMGPRGPPGATGPPGSPGFQGPPGEPGEPGQTGPAGSRGPAGPPGKAGEDGHPGKPGRPGERGVVGPQGARGFPGTPGLPGKGNQGHNGLDGLKGQPGAPGVKGEPGAPGENGTPGQIGARGLPGERGRVGGPGPAGARGSDGSVGPVGPAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGNPGANGLAGAKGAAGLPGVAGAPGLPGPRGIPGPVGAAGATGARGLVGEPGPAGSKGESGSKGEPGSAGPQGPPGPSGEEGKRGSSGEAGSAGPAGPPGLRGGPGSRGLPGADGRAGVMGPPGSRGASGPAGVRGPSGDSGRPGEPGVMGPRGLPGSPGNVGPAGKEGPAGLPGIDGRPGPIGPAGARGEPGNIGFPGPKGPAXXXXKNGDKGHAGLAGPRGAPGPDGNNGAQGPPGLQGVQGGKGEQGPAGPPGFQGLPGPSGTAGEAGKPGERGIPGEFGLPGPAGPRGERGPPGQSGAAGPTGPIGSRGPSGPPGPDGNKGEPGVVGAPGTAGPSGPGGLPGERGAAGIPGGKGEKGETGLRGDTGNTGRDGARGAPGAVGAPGPAGATGDRGEAGPAGSAGPAGPRGSPGERGEVGPAGPNGFAGPAGAAGQAGAKGERGTKGPKGENGPVGPTGPVGAAGPAGPNGPPGPAGSRGDGGPPGATGFPGAAGRTGPPGPAGITGPPGPPGAAGKEGLRGPRGDQGPVGRTGETGASGPPGFAGEKGSSGEPGTAGPPGTPGPQGILGPPGILGLPGSRGERGLPGVAGAVGEPGPLGIAGPPGARGPPGAVGSPGVNGAPGEAGRDGNPGSDGPPGRDGLPGHKGERGYPGNAGPVGTAGAPGPQGPLGPAGKHGNRGEPGPAGSVGPVGAVGPRGPSXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGQHGDQGSPGSVGPAGPRGPAGPSGPVGKDGRPGHAGAVGPAGVRGSQGSQGPSGPPGPPGPPGPPGPSGGGYDFGYDGDFYRADQPRSPPSLRPKDYEVDATLKSLNNQIETLLTPEGSKKNPARTCRDLRLSHPEWSSGYYWIDPNQGCTMDAIKVYCDFSTGETCIRAQPENIPAKNWYRSSKAKKHIWFGETINGGTQIEYNEEGVTTKDMATQLAFMRLLANHASQNITYHCKNSIAYMDEETGNLKKAVILQGSNDVELVAEGNSRFTYTVLVDGC",
				"SPPCLIESSEQGPTGPPGRDGEDGIPGPPGPPGPPGPPGLGGNFAAQYDAKGIGLGPGPMGLMGPRGPPGATGPPGSPGFQGPPGEPGEPGQTGPAGSRGPAGPPGKAGEDGHPGKPGRPGERGVVGPQGARGFPGTPGLPGKGNQGHNGLDGLKGQPGAPGVKGEPGAPGENGTPGQIGARGLPGERGRVGGPGPAGARGSDGSVGPVGPAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGNPGANGLAGAKGAAGLPGVAGAPGLPGPRGIPGPVGAAGATGARGLVGEPGPAGSKGESGSKGEPGSAGPQGPPGPSGEEGKRGSSGEAGSAGPAGPPGLRGGPGSRGLPGADGRAGVMGPPGSRGASGPAGVRGPSGDSGRPGEPGVMGPRGLPGSPGNVGPAGKEGPAGLPGIDGR**********GEPGNIGFPGPKGPAXXXXKNGDK*********GAPGPDGNNGAQGPPGLQGVQGGKGEQGPAGPPGFQGLPGPSGTAGEAGKPGERGIPGEFGLPGPAGPRGERGPPGQSGAAGPTGPIGSRGPSGPPGPDGNKGEPGVVGAPGTAGPSGPGGLPGERGAAGIPGGKGEKGETGLRGDTGNTGRDGARGAPGAVGAPGPAGATGDRGEAGPAGSAGPAGPRGSPGERGEVGPAGPNGFAGPAGAAGQAGAKGERGTKGPKGENGPVGPTGPVGAAGPAGPNGPPGPAGSRGDGGPPGATGFPGAAGRTGPPGPAGITGPPGPPGAAGKEGLRGPRGDQGPVGRTGETGASGPPGFAGEKGSSGEPGTAGPPGTPGPQGILGPPGILGLPGSRGERGLPGVAGAVGEPGPLGIAGPPGARGPPGAVGSPGVNGAPGEAGRDGNPGSDGPPGRDGLPGHKGERGYPGNAGPVGTAGAPGPQGPLGPAGKHGNR******************GPSXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGQHGDQGSPGSVGPAGPR***********DGRPGHAGAVGPAGVRGSQGSQGPSGPPGPPGPPGPPGPSGGGYDFGYDGDFYRADQPRSPPSLRPKDYEVDATLKSLNNQIETLLTPEGSKKNPARTCRDLRLSHPEWSSGYYWIDPNQGCTMDAIKVYCDFSTGETCIRAQPENIPAKNWYRSSKAKKHIWFGETINGGTQIEYNEEGVTTKDMATQLAFMRLLANHASQNITYHCKNSIAYMDEETGNLKKAVILQGSNDVELVAEGNSRFTYTVLVDGC",
			"/home/kwj/workspace/outputs/Mas/dbresults/elephant_mutated1_spec");
		*/
		
		/*
		 initDB();initPeptides();
			loadDB("/home/kwj/workspace/inputs/elephant_mutated2.fasta");
			//loadDB("/home/kwj/workspace/inputs/DataBases/ipi.HUMAN.fasta");
			loadPeptides(dir + "elephant_mutated2.out");
			writeOutput(dir  + "elephant_mutated2_spec");
			
			applyMutation(
					"GPKGDTGPRGPRGPAGPPGRDGIPGQPGLPGPPGPPGPPGPPGLGGNFAPQLSYGYDEKSAGGISVPGPMPSGLRGLLAPGAPXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGPTGPAGPPGFPGAVGAKGEAGPQGARGSEGPQGVRGEPGPPGPAGAAGPAGNPGADGQPGAKGANGAPGIAGAPGFPGARGPAGPQGPSGAPGPKGNSGEPGAPGSKGDAGAKGEPGPVGIQGPPGPAGEEGKRGARGEPGPTGLPGPPGERGGPGSRGFPGADGVAGPKGPAGERGSPGPAGPKGSPGEAGRPGEAGLPGAKGLTGSPGSPGPDGKTGPPVPAGQDGRPGPPVPGARGQAGVMGFPGPKGAAGEPGKAGERGVLGSSGAVXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGVPGDLGAPGPSGARGERGFPGERGVQGPPGPAGPRGSNGAPGNDGAKGDAGAPGAPGSQGAPGLQGMPGERGAAGLPGPKGDRGDAGPKGADGSPGKDGPRGLTGPIGPPGPAGAPGDKGEAGPSGPAGPTGARGAPGDRGEPGPPGPAGFAGPPGADGQPGAKGEPGDAGAKGDAGPPGPAGPTGAPGPIGNVGAPGPKGARGSAGPPGATGFPGAAGRVGPPGPSGNAGPPGPPGPAGKEGGKGPRGETGPAGRPGEVGPPGPPGPAGEKGSPGADGPAGAPGTPGPQGIGGQRGVVGLPGQRGERGFPGLPGPSGEPGKQGPSGSSGERGPPGPAGPPGLAGPPXXXXXXSAPDAESSPRRDGSPGPKSDRGETGPSEPPGAPGAPGAPDPVDPAGKSEDRGETGPAGPAGPAGPAGVRGPAGPQGPRGDKGETGEQGDRGLKGHRGFSGLQGPPGPPGSPGEQGPSGASGPAGPRGPPGSAGAPGKDGLNGLPGPIGPPGPRGRTGDAGPVGPPGPPGPPGPPGPPSGAFDFSFLPQPPQEKAHDGGRYYRADDANVVRDRDLEVDTTLKSLSQQIENIRSPEGSRKNPARTCRDLKMCHSDWKSGEYWIDPNQGCNLDAIKVFCNMETGETCVYPTQPSVVQNWYISNPKKRHVWYGESMTDEQFEYGGEGSDPADVAIQLTFLRLMSTEASQNITYHCKNSVAYMDQQTGNLKKALLLQGSNEIEIRAEGNSRFTYSVTEDGCTSHTGTWGKTIIEYKTTKTSRLPIIDVAPLDVGAPDQEFGFDIGPVCFL",
					"GPKGDTGPRGPRGPAGPPGRDGIPGQPGLPGPPGPPGPPGPPGLGGNFAPQLSYGYDEKSAGGISVPGPMPSGLRGLLAPGAPXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGPTGPAGPPGFPGAVGAKGEAGPQGARGSEGPQGVRGEPGPPGPAGAAGPAGNPGADGQPGAKGANGAPGIAGAPGFPGARGPAGPQGPSGAPGPKGNSGEPGAPGSKGDAGAKGEPGPVGIQGPPGPAGEEGKRGARGEPGPTGLPGPPGERGGPGSRGFPGADGVAGPKGPAGERGSPGPAGPKGSPGEAGRPGEAGLPGAKGLTGSPGSPGPDGKTGPPVPAGQDGRPGPPVPGARGQAGVMGFPGPKGAAGEPGKAGERGVLGSSGAVXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGVPGDLGAPGPSGARGERGFPGERGVQGPPGPAGPRGSNGAPGNDGAKGDAGAPGAPGSQGAPGLQGMPGERGAAGLPGPKGDRGDAGPKGADGSPGKDGPRGLTGPIGPPGPAGAPGDKGEAGPSGPAGPTGARGAPGDRGEPGPPGPAGFAGPPGADGQPGAKGEPGDAGAKGDAGPPGPAGPTGAPGPIGNVGAPGPKGARGSAGPPGATGFPGAAGRVGPPGPSGNAGPPGPPGPAGKEGGKGPRGETGPAGRPGEVGPPGPPGPAGEKGSPGADGPAGAPGTPGPQGIGGQRGVVGLPGQRGER**************QGPSGSSGERGPPGPAGPPGLAGPPXXXXXXSAPDAESSPRRDGSPGPKSDRGETGPSEPPGAPGAPGAPDPVDPAGKSEDRGETGPAGPAGPAGPAGVRGPAGPQGPRGDKGETGEQGDRGLKGHRGFSGLQGPPGPPGSPGEQGPSGASGPAGPRGPPGSAGAPGKDGLNGLPGPIGPPGPRGRTGDAGPVGPPGPPGPPGPPGPPSGAFDFSFLPQPPQEKAHDGGRYYRADDANVVRDRDLEVDTTLKSLSQQIENIRSPEGSRKNPARTCRDLKMCHSDWKSGEYWIDPNQGCNLDAIKVFCNMETGETCVYPTQPSVVQNWYISNPKKRHVWYGESMTDEQFEYGGEGSDPADVAIQLTFLRLMSTEASQNITYHCKNSVAYMDQQTGNLKKALLLQGSNEIEIRAEGNSRFTYSVTEDGCTSHTGTWGKTIIEYKTTKTSRLPIIDVAPLDVGAPDQEFGFDIGPVCFL",
					"/home/kwj/workspace/outputs/Mas/dbresults/elephant_mutated2_spec");
					*/
/*
		initDB();
		loadDB("/home/kwj/workspace/inputs/DataBases/Panda6Frame/Panda_found.fasta");
		initPeptides();
		loadPeptides(dir+"/Panda6Frame/");
		writeOutput(dir   + "/Panda6Frame/"+ outfileName);
		System.out.println("&"+ProteinSearch.numSpec+"&"+ProteinSearch.numPep+"&"+ProteinSearch.numPro1+"&"+ProteinSearch.numPro2+"\\\\");
	*/	
		/*
		initDB();initPeptides();
		loadDB("/home/kwj/workspace/inputs/DataBases/elephant_col.fasta");
		loadPeptides(dir + "elephant_col/");
		writeOutput(dir + "elephant_col/" + outfileName);
		*/
		
		/*
		initDB();initPeptides();
		loadDB("/home/kwj/workspace/inputs/DataBases/panda.fasta");
		loadPeptides(dir + "panda/");
		writeOutput(dir  + "panda/" + outfileName);
		*/
		/*
		initDB();initPeptides();
		loadDB("/home/kwj/workspace/inputs/uniprot_sprot.fasta");
		loadPeptides(dir + "uniprot/");
		writeOutput(dir  + "uniprot/" + outfileName);
		
		initDB();initPeptides();
		loadDB("/home/kwj/workspace/inputs/DataBases/ipi.MOUSE.fasta");
		loadPeptides(dir + "mouse/");
		writeOutput(dir + "mouse/" + outfileName);
		
		
		initDB();initPeptides();
		loadDB("/home/kwj/workspace/inputs/DataBases/ipi.HUMAN.fasta");
		loadPeptides(dir + "human/");
		writeOutput(dir + "human/" + outfileName);
		
		initDB();initPeptides();
		loadDB("/home/kwj/workspace/inputs/DataBases/DOG.fasta");
		loadPeptides(dir + "dog/");
		writeOutput(dir + "dog/" + outfileName);
		
		initDB();initPeptides();
		loadDB("/home/kwj/workspace/inputs/DataBases/ipi.RAT.fasta");
		loadPeptides(dir + "rat/");
		writeOutput(dir + "rat/" + outfileName);
		
		initDB();initPeptides();
		loadDB("/home/kwj/workspace/inputs/DataBases/ipi.CHICK.fasta");
		loadPeptides(dir + "chick/");
		writeOutput(dir + "chick/" + outfileName);
		
		*/
		
	//	initDB();initPeptides();
	//	loadDB("/home/kwj/workspace/inputs/uniprot_sprot.fasta");
	//	loadPeptides(dir + "uniprot/");
	//	loadDB("/home/kwj/workspace/inputs/DataBases/ipi.MOUSE.fasta");
	//	loadPeptides(dir + "mouse/");
	//	loadDB("/home/kwj/workspace/inputs/DataBases/ipi.HUMAN.fasta");
	//	loadPeptides(dir + "human/");
	//	loadDB("/home/kwj/workspace/inputs/DataBases/DOG.fasta");
	//	loadPeptides(dir + "dog/");
	//	loadDB("/home/kwj/workspace/inputs/DataBases/ipi.RAT.fasta");
	//	loadPeptides(dir + "rat/");
		//loadDB("/home/kwj/workspace/inputs/DataBases/ipi.CHICK.fasta");
	//	loadPeptides(dir + "chick/");
	//	writeOutput(dir  + outfileName);
		
	}
}
