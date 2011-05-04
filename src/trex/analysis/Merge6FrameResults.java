package trex.analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import parser.BufferedLineReader;

public class Merge6FrameResults {
	public static void generateMatchedProteinsDB(String outDir, String outPrefix , String fastaDir, String fastaPrefix) throws IOException{
		
		//String outPrefix = "panda", fastaPrefix = "panda_pro";
		int numDBs = 39;
		BufferedWriter mergedFasta = new BufferedWriter(new FileWriter(fastaDir + fastaPrefix + ".fasta"));
		BufferedWriter mergedOut = new BufferedWriter(new FileWriter(outDir + outPrefix + ".out"));
	//	Hashtable<String, String> proteinTable = new Hashtable<String, String>();
		
		for(int i=1; i<=numDBs; i++){
			//if(i==23) continue;
			System.out.println("processing " + i + "th file...");
			BufferedLineReader in = new BufferedLineReader(outDir + outPrefix + i + ".out");
			BufferedLineReader fasta = new BufferedLineReader(fastaDir + fastaPrefix + i + ".fasta");
			String s;
			
			ArrayList<String> peptides = new ArrayList<String>();
			
			while((s=in.readLine()) != null){
				mergedOut.write(s+"\n");
				if(s.startsWith("#") || s.startsWith("Spec")) continue;
				String[] t = s.split("\t");
				String toadd = t[4].toUpperCase().substring(t[4].indexOf('.')+1, t[4].lastIndexOf('.'));
				if(!peptides.contains(toadd))
					peptides.add(toadd);
			}
			
			String key = new String(), value = new String();
			
			while((s=fasta.readLine()) != null){
				if(s.startsWith(">")){
					String value_QI_replaced = value.replace('Q', 'K').replace('I', 'L');
					
					for(String pep : peptides){
						if(value_QI_replaced.contains(pep)){
							mergedFasta.write(key + "\n" +  value + "\n");
						//	proteinTable.put(key, value);//
							break;
						}
					}
					
					key = s;
					value = new String();
				}
				else{
					value += s;
				}	
			}
			
			String value_QI_replaced = value.replace('Q', 'K').replace('I', 'L');
			
			for(String pep : peptides){
				if(value_QI_replaced.contains(pep)){
					mergedFasta.write(key + "\n" +  value + "\n");
				//	proteinTable.put(key, value);//
					break;
				}
			}
			
			in.close(); fasta.close();
		}
		mergedFasta.close(); mergedOut.close();
	}
	
	public static void main(String[] argv) throws IOException{
		//generateMatchedProteinsDB("/home/kwj/workspace/outputs/cavebear/dbresults/Panda6Frame/" 
		//		, "panda" , "/home/kwj/workspace/inputs/DataBases/Panda6Frame/", "panda_pro");
		
		ProteinSearch.loadContaminantDB("/home/kwj/workspace/inputs/CommonContaminants.fasta");
		ProteinSearch.writeSpec = true;
		ProteinSearch.writePrecursor_Successor = false;
		ProteinSearch.initDB();ProteinSearch.initSpecAndPeptides();
		ProteinSearch.loadDB("/home/kwj/workspace/inputs/DataBases/Panda6Frame/panda_pro.fasta");
		ProteinSearch.loadSpecs("/home/kwj/workspace/outputs/cavebear/dbresults/Panda6Frame/panda.out");
		ProteinSearch.loadPeptides();
		ProteinSearch.writeOutput("/home/kwj/workspace/outputs/cavebear/dbresults/Panda6Frame/spec");
		
		
	}
	
}
