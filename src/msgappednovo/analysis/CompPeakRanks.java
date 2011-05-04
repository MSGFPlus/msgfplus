package msgappednovo.analysis;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;

import parser.MgfSpectrumParser;

import msgf.Tolerance;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;

public class CompPeakRanks {
	public static void write(String spec1, String spec2, String outfile, int num) throws IOException{
		int charge = 2;
		int maxRank = 200;
		
		Iterator<Spectrum> iterator;
		Iterator<Spectrum> iterator2;
		PrintStream out = new PrintStream(outfile);
		
		iterator = new SpectraIterator(spec1, new MgfSpectrumParser());
		//	idealiterator = new SpectraIterator(idealmgf, new MgfSpectrumParser());
		iterator2 = new SpectraIterator(spec2, new MgfSpectrumParser());
		
		int n=0;
		out.println("ranks = [");
		while(iterator.hasNext()){
			Spectrum so = iterator.next();
			if(so.getAnnotation().isModified()) continue;
		//	Spectrum si = idealiterator.next();
			
			if(so.getCharge() != charge) continue;
		
			if(n ++ > num) break;
			
			Spectrum sr = iterator2.next();
			System.out.println(sr.getScanNum()+ " " +so.getAnnotation() + "\t" + sr.getAnnotation());
			assert(so.getAnnotation().equals(sr.getAnnotation()));
			
			so.setRanksOfPeaks(); sr.setRanksOfPeaks();
			
			for(Peak p : so){
				int[] rank = new int[2];
				
				rank[0] = p.getRank();
				if(rank[0] > maxRank) continue;
				rank[1] = sr.getPeakByMass(p.getMz(), new Tolerance(0, false)).getRank();
				if(rank[1] > maxRank) continue;
				out.println(rank[0] + "\t" + rank[1]);
			}
			
			
		}
		out.println("];");
		out.close();
	}
	
	static public void main(String[] args) throws IOException{
			
			int charge = 2; 

			String inputmgf = "/home/kwj/workspace/inputs/Training/CID_Tryp_Confident.mgf";
			String outmgf = inputmgf.substring(0, inputmgf.lastIndexOf(".")) + "_ranked.mgf";
			String outm = "/home/kwj/workspace/inputs/Training/prank.m";
				
			write(inputmgf, outmgf, outm, 200);
			
	}
}
