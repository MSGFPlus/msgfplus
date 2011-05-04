package trex.analysis;

import java.io.IOException;

public class Figure1 {
	public static void main(String[] argv) throws IOException{
		String species = "dog";
		boolean isShort = true;
		
		ProteinSearch.setDataSet("cavebear");
		
		
		String outfileName = "pro_for_figure1";
		ProteinSearch.writeProtein = true;
		//ProteinSearch.writePeptide = true;
		ProteinSearch.loadContaminantDB("/home/kwj/workspace/inputs/CommonContaminants.fasta");
		
		ProteinSearch.initDB();
		ProteinSearch.initSpecAndPeptides();
		// xxColSeq file should contain whole seq
		ProteinSearch.loadDB("/home/kwj/workspace/outputs/"+ProteinSearch.DataSet + "/dbresults/" + species +"ColSeq");
		if(isShort)
			ProteinSearch.loadSpecs("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/dogpanda_short/");
		else
			ProteinSearch.loadSpecs("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/" + species + "/");
		
		
		ProteinSearch.loadPeptides();
		
		if(isShort){
			ProteinSearch.writeOutput("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated_short/" + outfileName);
			ProteinSearch.writePeptide = true;
			ProteinSearch.writeOutput("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated_short/" + outfileName+"_spec");
		}
		else{
			ProteinSearch.writeOutput("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated/" + outfileName);
			ProteinSearch.writePeptide = true;
			ProteinSearch.writeOutput("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated/" + outfileName+"_spec");
		}
		ProteinSearch.writePeptide = false;
		ProteinSearch.writeProtein = false;
		
		if(isShort)
			ProteinSearch.generateMutatedProteinDB("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated_short/" + outfileName,
					"/home/kwj/workspace/inputs/" + species + "col.fasta");
		else
			ProteinSearch.generateMutatedProteinDB("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated/" + outfileName,
					"/home/kwj/workspace/inputs/" + species + "col.fasta");

		// do MSgap.jar to generate .out
		
		String outfileName2 = "spec_for_figure1";
		ProteinSearch.writeSpec = true;
		
		ProteinSearch.loadContaminantDB("/home/kwj/workspace/inputs/CommonContaminants.fasta");
		
		ProteinSearch.initDB();
		ProteinSearch.initSpecAndPeptides();
		
		ProteinSearch.loadDB("/home/kwj/workspace/inputs/" + species + "col.fasta");
		
		if(isShort){
			ProteinSearch.loadSpecs("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/" + species + "col_mutated_short/");
			ProteinSearch.loadPeptides();
			ProteinSearch.writeOutput("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated_short/" + outfileName2);
			
			ProteinSearch.applyMutation(
					"/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated_short/" + outfileName
					, "/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated_short/" + outfileName2);

		}else{
			ProteinSearch.loadSpecs("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/" + species + "col_mutated/");
			ProteinSearch.loadPeptides();
			ProteinSearch.writeOutput("/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated/" + outfileName2);
			
			ProteinSearch.applyMutation(
					"/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated/" + outfileName
					, "/home/kwj/workspace/outputs/" + ProteinSearch.DataSet + "/dbresults/"  + species + "col_mutated/" + outfileName2);
		}
	}
}
