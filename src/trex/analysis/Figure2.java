package trex.analysis;

import java.io.IOException;

import parser.BufferedLineReader;

public class Figure2 {
	public static void main(String[] argv) throws IOException{
		ProteinSearch.setDataSet("cavebear");
		BufferedLineReader in = null;
		in = new BufferedLineReader("/home/kwj/workspace/outputs/"+ ProteinSearch.DataSet + "/dbresults/col/collagen");
		String s;
		int numD = 0, numH=0, numB=0, numM=0, numC=0, numR=0, numO=0, numP=0;
		int numDH=0, numDB=0, numDM=0, numDC=0, numDR=0, numDO = 0, numDP = 0;
		int numPH=0, numPB=0, numPM=0, numPC=0, numPR=0,numPO=0, numPD = 0;
		while((s=in.readLine()) != null){
			String t = s.split("&")[4];
			float specProb = Float.parseFloat(s.split("&")[3]);
			//System.out.println(specProb);
		//	if(specProb > 1e-11) continue;
			if(t.contains("D")){
				numD++;
				if(t.contains("H")) numDH++;
				if(t.contains("B")) numDB++;
				if(t.contains("C")) numDC++;
				if(t.contains("R")) numDR++;
				if(t.contains("M")) numDM++;
				if(t.contains("O")) numDO++;
				if(t.contains("P")) numDP++;
			}
			if(t.contains("H")) numH++;
			if(t.contains("B")) numB++;
			if(t.contains("C")) numC++;
			if(t.contains("R")) numR++;
			if(t.contains("M")) numM++;
			if(t.contains("O")) numO++;
			if(t.contains("P")){
				numP++;
				if(t.contains("D")) numPD++;
				if(t.contains("H")) numPH++;
				if(t.contains("B")) numPB++;
				if(t.contains("C")) numPC++;
				if(t.contains("R")) numPR++;
				if(t.contains("M")) numPM++;
				if(t.contains("O")) numPO++;
			}
		}
		System.out.println("DH\t"+(numD-numDH)+"\t"+numDH+"\t"+(numH-numDH));
		System.out.println("DM\t"+(numD-numDM)+"\t"+numDM+"\t"+(numM-numDM));
		System.out.println("DR\t"+(numD-numDR)+"\t"+numDR+"\t"+(numR-numDR));
		System.out.println("DC\t"+(numD-numDC)+"\t"+numDC+"\t"+(numC-numDC));
		System.out.println("DB\t"+(numD-numDB)+"\t"+numDB+"\t"+(numB-numDB));
		System.out.println("DO\t"+(numD-numDO)+"\t"+numDO+"\t"+(numO-numDO));
		System.out.println("DP\t"+(numD-numDP)+"\t"+numDP+"\t"+(numP-numDP));
		System.out.println();
		System.out.println("PH\t"+(numP-numPH)+"\t"+numPH+"\t"+(numH-numPH));
		System.out.println("PM\t"+(numP-numPM)+"\t"+numPM+"\t"+(numM-numPM));
		System.out.println("PR\t"+(numP-numPR)+"\t"+numPR+"\t"+(numR-numPR));
		System.out.println("PC\t"+(numP-numPC)+"\t"+numPC+"\t"+(numC-numPC));
		System.out.println("PB\t"+(numP-numPB)+"\t"+numPB+"\t"+(numB-numPB));
		System.out.println("PO\t"+(numP-numPO)+"\t"+numPO+"\t"+(numO-numPO));
		System.out.println("PD\t"+(numP-numPD)+"\t"+numPD+"\t"+(numD-numPD));
		
		
		in.close();
	}
}
