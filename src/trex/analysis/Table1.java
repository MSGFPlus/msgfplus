package trex.analysis;

import java.io.IOException;

public class Table1 {
	
	static String[] columnName = {
		"Chicken",
		"Dog", 
		"Human", 
		"Bovin",
		"Mouse", 
		"Rat", 
		
		"SwissProt",
		"Panda",
	};
	
	
	static int numDB = ProteinSearch.dbName.length;
	
	public static void main(String[] argv) throws IOException{
		ProteinSearch.setDataSet("cavebear");
		//ProteinSearch.setDataSet("Mas");
		//ProteinSearch.setDataSet("TRex");
		//ProteinSearch.setDataSet("Had");
		System.out.println("\\begin{table}");
		System.out.println("\\center\\small");
		System.out.println("\\begin{tabular}{|c|c|c|c|c|}");
		System.out.println("\\hline");
		              
		System.out.println("Database&\\# matching&\\# matched&\\# matched proteins&\\# matched proteins\\\\");
		System.out.println("&spectra&peptides&(matched by at least 1 peptide)&(matched by at least 2 peptides)\\\\\\hline\\hline");
		
		
		
		ProteinSearch.loadContaminantDB("/home/kwj/workspace/inputs/CommonContaminants.fasta");
		
		for(int i=0; i<numDB;i++){
			ProteinSearch.initDB(); ProteinSearch.initSpecAndPeptides();
			ProteinSearch.loadDB(ProteinSearch.dbName[i]);
			ProteinSearch.loadSpecs(ProteinSearch.dirName[i]);
			ProteinSearch.loadPeptides();
			ProteinSearch.writeOutput();
			System.out.println(columnName[i] + "&"+ProteinSearch.numSpec+"&"+ProteinSearch.numPep+"&"+ProteinSearch.numPro1+"&"+ProteinSearch.numPro2+"\\\\");
		}
		System.out.println("\\hline");
		
		System.out.println("\\end{tabular}");
		System.out.println("\\vspace{1cm}");
		System.out.println("\\caption{}\\label{}");
		System.out.println("\\end{table}");
	}
}
