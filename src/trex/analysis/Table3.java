package trex.analysis;

import java.io.IOException;
import java.util.Hashtable;

import parser.BufferedLineReader;

public class Table3 {
	static int numDB = ProteinSearch.dbName.length;
	
	static Hashtable<String, String> nonCol = new Hashtable<String, String>();
	
	static void updateNonCol(String filename) throws IOException{
		BufferedLineReader in = null;
		in = new BufferedLineReader(filename);
		String s;
		String key = new String(), value = new String();
		
		while((s=in.readLine()) != null){
			if(s.startsWith(">")){
				if(!key.isEmpty() && !key.contains("p")&& !nonCol.containsKey(key)){
					nonCol.put(key, value);
				}
				value = s.substring(1, s.indexOf('#'));
				key = new String();
			}else{
				if(key.isEmpty())
					key = s.replace('\t', '&');
				else
					key += "\\\\\n" + s.replace('\t', '&');
				//System.out.println(key);
			}	
		}
		
		if(!key.isEmpty() &&!key.contains("p")&& !nonCol.containsKey(key)){
			nonCol.put(key, value);
		}
	
	}
	
	public static void main(String[] argv) throws IOException{
		ProteinSearch.setDataSet("Had");
	//	ProteinSearch.setDataSet("TRex");
	//	ProteinSearch.setDataSet("cavebear");
		System.out.println("\\begin{table}[t]");
		System.out.println("\\center\\small");
		System.out.println("\\begin{tabular}{|c|c|c|c|}");
		System.out.println("\\hline");
		System.out.println("Gene Symbol / Protein name & Sequence & Match score & Spec. Prob \\\\\\hline\\hline");
		
		
		/*
		ProteinSearch.writeSpec = true;
		ProteinSearch.onlyNONCOL = true;
		ProteinSearch.allowHydroxyProlineforNonCOL = false;
		ProteinSearch.initDB(); ProteinSearch.initSpecAndPeptides();
		
		ProteinSearch.loadContaminantDB("/home/kwj/workspace/inputs/CommonContaminants.fasta");
		
		for(int i=0; i<numDB;i++){ // dog -> by hand
			if(!ProteinSearch.dbName[i].contains("DOG")){
				ProteinSearch.loadDB(ProteinSearch.dbName[i]);
				ProteinSearch.loadSpecs(ProteinSearch.dirName[i]);
				ProteinSearch.loadPeptides();
			}
		}
		
		ProteinSearch.writeOutput("/home/kwj/workspace/outputs/"+ ProteinSearch.DataSet + "/dbresults/non_col/non_col.txt");
		
		updateNonCol("/home/kwj/workspace/outputs/"+ ProteinSearch.DataSet + "/dbresults/non_col/non_col.txt");
		*/
		
		updateNonCol("/home/kwj/workspace/outputs/"+ ProteinSearch.DataSet + "/dbresults/non_col/non_col_with_dog.txt");
		
		for(String key : nonCol.keySet()){
			float specProb = Float.parseFloat(key.substring(key.lastIndexOf('&')+1));
			System.out.print(nonCol.get(key) + key.substring(0, key.lastIndexOf('&')+1));
			System.out.printf("%1.2e \\\\\n", specProb);
		}
		
		System.out.println("\\hline");
		System.out.println("\\end{tabular}");
		System.out.println("\\vspace{1cm}");
		System.out.println("\\caption{}\\label{}");
		System.out.println("\\end{table}");
	}
}
