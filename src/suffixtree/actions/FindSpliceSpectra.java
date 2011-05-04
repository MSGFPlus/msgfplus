package suffixtree.actions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.SortedSet;
import java.util.TreeSet;

import msgap.results.GappedPeptideResults;
import msgap.results.MSGDResultFileParser;

public class FindSpliceSpectra {
	
	private static class MatchInfo {
	  private int scanNum;
	  private String filename;
	  private float prob;
	  private String id;
	  private int start, mid, end;
	  private int intronSize;
	  
	  public MatchInfo(int scanNum, String filename, float prob, String id, int start, int end, int mid) {
	    this.scanNum = scanNum;
	    this.filename = filename;
	    this.prob = prob;
	    this.id = id;
	    this.start = start;
	    this.mid = mid;
	    this.end = end;
	  }
	 
	  @Override
	  public String toString() {
	    int perc = getLocationPerc();
	    return String.format("%s\t%d\t%.2e\t%s\t%d\t%d", this.filename, this.scanNum, this.prob, this.id, perc, intronSize);
	  }
	  
	  public int getLocationPerc() {
	    return (this.mid - this.start) * 100 / (this.end - this.start);
	  }
	  
	}
	
	
	 public static void processHumanResults(String resultFile, String grcFile, String outGrcFile, String solutionFile) {
	    
	    ArrayList<MatchInfo> matches = new ArrayList<MatchInfo>();
	    
	    try {
	      BufferedReader in = new BufferedReader(new FileReader(resultFile));
	      
	      String line;
	      while ((line = in.readLine()) != null) {
	        //System.out.println(line);
	        
	        String[] tokens = line.split("\t");
	        String annotation = tokens[6];
	        
	        if (tokens.length < 13) continue;
	        
	        // parse the coordinates, coordinates are inclusive
	        String[] annoTokens = annotation.split(";");
	        String[] coors = annoTokens[0].split(":");
	        String[] intronStr = annoTokens[1].split(":");
	        
	        if (coors.length>1) { // the last coor is bogus because it is the end
	          TreeSet<Integer> breaks = new TreeSet<Integer>();
	          HashMap<Integer,Integer> introns = new HashMap<Integer,Integer>();
	          for (int index=0; index<coors.length-1; index++) {
              int pos = Integer.parseInt(coors[index]);
              if (pos>=0) {
	              breaks.add(pos);
	              introns.put(pos, Integer.parseInt(intronStr[index]));
              }
	          }
	          
	          int start = Integer.parseInt(tokens[9]);
	          int end = Integer.parseInt(tokens[10]);
	          
	          //System.err.println(start + " " + end + " " + line);
	          SortedSet<Integer> mids = breaks.subSet(start, false, end, false);
	          if (mids.size()==1) {
              MatchInfo mi = new MatchInfo(Integer.parseInt(tokens[1]), tokens[0], Float.parseFloat(tokens[8]), tokens[5], start, end, mids.first());
              int percCut = mi.getLocationPerc();
              if (40 < percCut && percCut < 60 && (mi.end-mi.start) >= 25) {
                mi.intronSize = introns.get(mids.first());
                matches.add(mi);
              }
	          }
	          
	        }
	      }
	      
	      PrintWriter solutionFD = new PrintWriter(solutionFile);
	       
	      GappedPeptideResults gpr = new MSGDResultFileParser(grcFile, Integer.MAX_VALUE).iterator().next();
	      HashSet<Integer> retain = new HashSet<Integer>();
	      for (MatchInfo mi : matches) {
	        solutionFD.println(mi.toString());
	        int specId = gpr.getSpecId(mi.filename, mi.scanNum);
	        if (specId < 0) {
	          System.out.println(mi.toString());
	        }
	        else {
	          retain.add(specId);
	        }
	      }
	      solutionFD.close();
	      
	      GappedPeptideResults selectedGpr = gpr.retain(retain);
	      System.out.println("Total selected spectra " + retain.size());
	      
	      PrintWriter newGrcFile = new PrintWriter(outGrcFile);
	      selectedGpr.toFile(newGrcFile);
	      newGrcFile.close();
	      
	      
	      
	    }
	    catch (IOException ioe) {
	      System.err.println(ioe);
	      System.exit(-1);
	    }
	    
	  }
	 
	 
	 
	 
	public static void processResults(String resultFile, String grcFile, String outGrcFile, String solutionFile) {
		
	  ArrayList<MatchInfo> matches = new ArrayList<MatchInfo>();
	  
		try {
			BufferedReader in = new BufferedReader(new FileReader(resultFile));
			
			String line;
			while ((line = in.readLine()) != null) {
				//System.out.println(line);
				
				String[] tokens = line.split("\t");
				String annotation = tokens[6];
				
				// parse the coordinates, coordinates are inclusive
				String[] annoTokens = annotation.split(", ");
				String[] coors = annoTokens[1].split("from ")[1].split(",");
				
				if (coors.length>1) {
					ArrayList<Integer> breaks = new ArrayList<Integer>();
					int cumCount = 0;
					for (String coorPair : coors) {
					  String[] pair = coorPair.split("-");
	  				cumCount += Math.abs(Integer.parseInt(pair[1])-Integer.parseInt(pair[0]))+1;
      			breaks.add(cumCount/3);
					}
					breaks.remove(breaks.size()-1); // the last break is not a splice break
					
					int start = Integer.parseInt(tokens[9]);
	        int end = Integer.parseInt(tokens[10]);
	        
	        for (int b : breaks) {
	          if (start < b && b < end) {
	            MatchInfo mi = new MatchInfo(Integer.parseInt(tokens[1]), tokens[0], Float.parseFloat(tokens[8]), tokens[5], start, end, b);
	            int percCut = mi.getLocationPerc();
	            if (40 < percCut && percCut < 60) {
	              matches.add(mi);
	            }
	          }
	        }
	        
				}
			}
			
			System.out.println("Total selected spectra " + matches.size());
			
			PrintWriter solutionFD = new PrintWriter(solutionFile);
		   
		  GappedPeptideResults gpr = new MSGDResultFileParser(grcFile, Integer.MAX_VALUE).iterator().next();
		  HashSet<Integer> retain = new HashSet<Integer>();
		  for (MatchInfo mi : matches) {
		    solutionFD.println(mi.toString());
		    int specId = gpr.getSpecId(mi.filename, mi.scanNum);
		    if (specId < 0) {
		      System.out.println(mi.toString());
		    }
		    else {
		      retain.add(specId);
		    }
		  }
		  solutionFD.close();
		  
		  GappedPeptideResults selectedGpr = gpr.retain(retain);
		  
		  PrintWriter newGrcFile = new PrintWriter(outGrcFile);
		  selectedGpr.toFile(newGrcFile);
		  newGrcFile.close();
		  
		  
			
		}
		catch (IOException ioe) {
			System.err.println(ioe);
			System.exit(-1);
		}
		
	}

	public static void runYeast() {
  	String userHome = System.getProperty("user.home");
    String resultFile = userHome + "/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/results6ORF.txt";
    String grcFile = userHome + "/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/results.grc";
    String outGrcFile = userHome + "/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/resultsSplicing.grc";
    String solutionFile = userHome + "/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/solution.txt";
    processResults(resultFile, grcFile, outGrcFile, solutionFile);
	}
	
	public static void runHuman() {
    String userHome = System.getProperty("user.home");
    String resultFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/lys/results6ORF.txt";
    String grcFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/lys/results6.grc";
    String outGrcFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/lys/resultsSplicing.grc";
    String solutionFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/lys/solution.txt";
    processHumanResults(resultFile, grcFile, outGrcFile, solutionFile);
  }
	
	public static void main(String[] args) {
	  //runYeast();
	  runHuman();
	}
	
	
}
