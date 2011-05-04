package trex.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;

public class Table2 {
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
		ProteinSearch.setDataSet("TRex");
	//	ProteinSearch.setDataSet("Mas");
		ProteinSearch.setDataSet("Had");
	//	ProteinSearch.setDataSet("cavebear");
		ProteinSearch.loadContaminantDB("/home/kwj/workspace/inputs/CommonContaminants.fasta");
		
		ProteinSearch.writeSpec = true;
		ProteinSearch.initDB();ProteinSearch.initSpecAndPeptides();
		ProteinSearch.loadDB("/home/kwj/workspace/inputs/DataBases/DOG.fasta");
		ProteinSearch.loadSpecs("/home/kwj/workspace/outputs/"+ProteinSearch.DataSet+"/dbresults/dog/");
		ProteinSearch.loadPeptides();
	
		ProteinSearch.writeOutput("/home/kwj/workspace/outputs/"+ProteinSearch.DataSet+"/dbresults/dog/spec");
		ProteinSearch.writeSpec = false;
		
		System.out.println("\\begin{table}");
		System.out.println("\\center\\small");
		System.out.println("\\begin{tabular}{|c|c|c|c|}");

		System.out.println("\\hline");
		              
		System.out.println("Database&Collagen/non-collagen&\\# proteins&\\# proteins\\\\");
		System.out.println("&&(matched by at least 1 peptide)&(matched by at least 2 peptides)\\\\\\hline\\hline");
		
		
		
		ProteinSearch.onlyCOL = true;
		
		
		
		int[][] numPro = new int[4][numDB];
		
		for(int i=0; i<numDB;i++){
			ProteinSearch.initDB(); ProteinSearch.initSpecAndPeptides();
			ProteinSearch.loadDB(ProteinSearch.dbName[i]);
			ProteinSearch.loadSpecs(ProteinSearch.dirName[i]);
			ProteinSearch.loadPeptides();

			ProteinSearch.writeOutput();
			
			numPro[0][i] = ProteinSearch.numPro1;
			numPro[1][i]= ProteinSearch.numPro2;
		}
		
		ProteinSearch.onlyCOL = false;
		ProteinSearch.onlyNONCOL = true;
	//	ProteinSearch.allowHydroxyProlineforNonCOL = false;
		
		for(int i=0; i<numDB;i++){
			ProteinSearch.initDB(); ProteinSearch.initSpecAndPeptides();
			ProteinSearch.loadDB(ProteinSearch.dbName[i]);
			ProteinSearch.loadSpecs(ProteinSearch.dirName[i]);
			ProteinSearch.loadPeptides();
			ProteinSearch.writeOutput();
			
			numPro[2][i] = ProteinSearch.numPro1;
			numPro[3][i]= ProteinSearch.numPro2;
		}
		
		for(int i=0; i<numDB;i++){
			System.out.println(columnName[i] + "&"+"collagen&" + numPro[0][i]+"&"+ numPro[1][i]+"\\\\");
			System.out.println("&non-collagen&" + numPro[2][i]+"&"+ numPro[3][i]+"\\\\\\hline");
		}
		
		System.out.println("\\end{tabular}");
		System.out.println("\\vspace{1cm}");
		System.out.println("\\caption{}\\label{}");
		System.out.println("\\end{table}");
		//System.out.println("\\hline");
	}
}
