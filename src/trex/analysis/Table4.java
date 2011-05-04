package trex.analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.Hashtable;

import parser.BufferedLineReader;

public class Table4 {
	static int numDB = ProteinSearch.dbName.length;
	
	static String[] name = {
		"Chicken",
		"Dog", 
		"Human", 
		"Bovin",
		"Mouse", 
		"Rat", 
		"SwissProt",
		"Panda",
		};
	
	public static void main(String[] argv) throws IOException{
		ProteinSearch.setDataSet("Had");
	    //ProteinSearch.setDataSet("Mas");
		//ProteinSearch.setDataSet("TRex");
	//	ProteinSearch.setDataSet("cavebear");
		System.out.println("\\begin{table}");
		System.out.println("\\center\\small");
		System.out.println("\\begin{tabular}{|c|c|c|c|c|}");
		System.out.println("\\hline");
		
		ProteinSearch.loadContaminantDB("/home/kwj/workspace/inputs/CommonContaminants.fasta");
		ProteinSearch.writeSpec = true;
		ProteinSearch.writePrecursor_Successor = false;
		ProteinSearch.writeScanNum = true;
		ProteinSearch.initDB();ProteinSearch.initSpecAndPeptides();
		ProteinSearch.loadDB("/home/kwj/workspace/inputs/DataBases/DOG.fasta");
		ProteinSearch.loadSpecs("/home/kwj/workspace/outputs/"+ProteinSearch.DataSet+"/dbresults/dog/");
		ProteinSearch.loadPeptides();
		ProteinSearch.writeOutput("/home/kwj/workspace/outputs/"+ProteinSearch.DataSet+"/dbresults/dog/spec_for_col");
		ProteinSearch.writeSpec = false;
		
		
		BufferedWriter out = null;
		out =  new BufferedWriter(new FileWriter("/home/kwj/workspace/outputs/"+ ProteinSearch.DataSet + "/dbresults/col/collagen"));
		
		System.out.println("Scan No. & Sequence & Match score & Spec. Prob &DBs that peptides are found in\\\\\\hline\\hline");
		//out.write("Sequence& Sequence & Match score & Spec. Prob &DBs that peptides are found in\\\\\\hline\\hline\n");
		
		
		ProteinSearch.writeSpec = true;
		ProteinSearch.writePrecursor_Successor = false;
		//ProteinSearch.writePeptide = true;
		ProteinSearch.onlyCOL = true;
		
		
		//ProteinSearch.loadContaminantDB("/home/kwj/workspace/inputs/CommonContaminants.fasta");
		
		ProteinSearch.writeScanNum = true;
		ProteinSearch.initSpecAndPeptides();
		/*
		for(int i=0; i<numDB;i++){
			ProteinSearch.loadSpecs(ProteinSearch.dirName[i]);
			//System.out.println(ProteinSearch.specSet.size()+" " +ProteinSearch.peptideSet.size());
		}
		
		ProteinSearch.loadPeptides();
		
		for(int i=0; i<numDB;i++){ // dog -> by hand
			if(!ProteinSearch.dbName[i].contains("DOG")) ProteinSearch.onlyCOL = true;
			else ProteinSearch.onlyCOL = false;
			ProteinSearch.initDB();
			ProteinSearch.loadDB(ProteinSearch.dbName[i]);
		//	ProteinSearch.loadPeptides(ProteinSearch.dirName[i]);
			//System.out.println(ProteinSearch.specSet.size()+" " +ProteinSearch.peptideSet.size());
			ProteinSearch.writeOutput("/home/kwj/workspace/outputs/"+ ProteinSearch.DataSet + "/dbresults/col/"+name[i] +"col.txt");
			//System.out.println(ProteinSearch.specSet.size()+" " +ProteinSearch.peptideSet.size());
		}
		
		
		*/
		Hashtable<String, BitSet> pep_DB = new Hashtable<String, BitSet>();
		
		for(int i=0; i<numDB;i++){
			//if(ProteinSearch.dbName[i].contains("DOG")) continue;
			//System.out.println(numDB);
			String filename = "/home/kwj/workspace/outputs/"+ ProteinSearch.DataSet + "/dbresults/col/"+name[i] +"col.txt";
			BufferedLineReader in = null;
			in = new BufferedLineReader(filename);
			String s;
			boolean otherSpecies = true;
			while((s = in.readLine()) != null){
				if(s.startsWith(">")){
					if( i<numDB && ProteinSearch.dbName[i].contains("uniprot_sprot")){
						if(s.contains("Canis familiaris ") || s.contains("Bos taurus")|| s.contains("Homo sapiens")
								||s.contains("Rattus norvegicus")||s.contains("Mus musculus ") || s.contains("Gallus gallus"))
							otherSpecies = false;
						else otherSpecies = true;
					}else otherSpecies = true;
					continue;
				}
				
				BitSet v;
				if(otherSpecies){
					if(!pep_DB.containsKey(s)){
						v = new BitSet();
					}else{
						v = pep_DB.get(s);
					}
					//System.out.println(v.toString());
					v.set(i);
					
					pep_DB.put(s, v);
				}
				
			}
			in.close();
		}
		
		for(String s:pep_DB.keySet()){
			String p1 = s.substring(1).split("\t")[1].toLowerCase();
			BitSet v1 = pep_DB.get(s);
			for(String t:pep_DB.keySet()){
				if(t.equals(s)) continue;
				String p2 = t.substring(1).split("\t")[1].toLowerCase();
				BitSet v2 = pep_DB.get(t);
				
				if(p1.contains(p2)){
					v2.or(v1);
					pep_DB.put(t, v2);
				}
				
				if(p2.contains(p1)){
					v1.or(v2);
					pep_DB.put(s, v1);
					
				}
			}
		}
		
		char[] dbInitial = {'C','D','H','B','M','R','O','P'};
		
		ArrayList<String> usedkey = new ArrayList<String>();
		ArrayList<String> multiplePepScanNum = new ArrayList<String>();
		HashSet<String> tmpScanNum = new HashSet<String>();
		

		for(String s: pep_DB.keySet()){
			String[] p = s.substring(1).split("\t");
			tmpScanNum.add(p[0]);
			multiplePepScanNum.add(p[0]);
		}
		
		for(String toremove : tmpScanNum)
			multiplePepScanNum.remove(toremove);
		
		int minHashCode = Integer.MAX_VALUE, maxHashCode = Integer.MIN_VALUE;
		for(String s: pep_DB.keySet()){
			int hashCode = pep_DB.get(s).hashCode();
			minHashCode = minHashCode < hashCode? minHashCode:hashCode;
			maxHashCode = maxHashCode > hashCode? maxHashCode:hashCode;
		}
		
		for(int l=name.length;l>0;l--){
			for(int i=minHashCode;i<=maxHashCode;i++){
	
				for(String s: pep_DB.keySet()){
					//System.out.println(s);
				//	pep_DB.get(s).hashCode()
					if(pep_DB.get(s).cardinality() != l) continue;
					if(pep_DB.get(s).hashCode() != i) continue;
					if(usedkey.contains(s)) continue;
					
					//System.out.println(pep_DB.size() + " " + cntr);
					
					
					//System.out.println(usedkey.size());
					String[] p = s.substring(1).split("\t");
					float specProb = Float.parseFloat(p[3]);
					
					if(multiplePepScanNum.contains(p[0]))
						System.out.print("${}^{" + p[0] + "}$");
					
				//	System.out.print(p[0] + "&"+ p[1] + "&"+ p[2]+"&");
					System.out.print( p[1] + "&"+ p[2]+"&");
					System.out.printf("%1.2e&", specProb);
	
					out.write(s.substring(1).replace('\t', '&')+"&");
					String dbs = new String();
					for(int j=0; j<name.length;j++){
						if(pep_DB.get(s).get(j)) dbs+=dbInitial[j];
					}
					System.out.println(dbs+"\\\\");
					out.write(dbs+"\\\\\n");
					usedkey.add(s);
	
				
				}
			}
		}
		System.out.println("\\hline");
		//out.write("\\hline\n");
		out.close();
		//System.out.println(pep_DB);
		
		/*updateNonCol("/home/kwj/workspace/outputs/"+ ProteinSearch.DataSet + "/dbresults/non_col/non_col.txt");
		
		for(String key : nonCol.keySet()){
			float specProb = Float.parseFloat(key.substring(key.lastIndexOf('&')+1));
			System.out.print(nonCol.get(key) + key.substring(0, key.lastIndexOf('&')+1));
			System.out.printf("%1.2e \\\\\n", specProb);
		}
		
		System.out.println("\\hline");
		*/
		
		System.out.println("\\end{tabular}");
		System.out.println("\\vspace{1cm}");
		System.out.println("\\caption{}\\label{}");
		System.out.println("\\end{table}");
	}
}
