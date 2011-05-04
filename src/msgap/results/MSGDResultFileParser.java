package msgap.results;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;

import suffixtree.Constants;

/**
 * This class understand the format of the .grc files of the MSGappedDictionary
 * output. 
 * @author jung
 *
 */
public class MSGDResultFileParser implements Iterable<GappedPeptideResults> {

  private class MyIterator implements Iterator<GappedPeptideResults> {

    private String nextLine;
    private BufferedReader fid;
    
    private MyIterator() {
      try {
        fid = new BufferedReader(new FileReader(filePath));
      }
      catch (IOException ioe) {
        System.err.println(ioe);
        System.exit(-9);
      }
    }
    
    @Override
    public boolean hasNext() {
      
      // we have more to parse
      if (nextLine!=null) return true;
      
      try {
        nextLine = fid.readLine();
        if (nextLine!=null) return true;
      }
      catch (IOException ioe) {
        System.err.println(ioe);
        System.exit(-9);
      }
      
      return false;
    }

    @Override
    public GappedPeptideResults next() {
      
      // we have initialized nextLine by calling the check
      if (!hasNext()) return null;
      
      GappedPeptideResults gpr = new GappedPeptideResults();
      
      // parse the file now
      String line = this.nextLine, ident = null, filename = null, actMethod = null;
      boolean added = false;
      int currentCharge = 0;
      int currentSpecId = -1;
      int currentScanNumber = -1;
      int queriesCount = 0, discarded = 0;
      float currentPm = 0.0f;
      
      try {
        
        do {
          if (line.charAt(0)=='#') {
            
            if (queriesCount>=chunkSize) {
              // start parsing from here next time the iterator is called
              nextLine = line;
              return gpr;
            }
            
            added = false; 
            
            // parse the identification line
            StringTokenizer tk = new StringTokenizer(line);
            
            // token 1 is the scan number, after discarding the first char
            currentSpecId = Integer.parseInt(tk.nextToken().substring(1));
            
            // scanNumber
            currentScanNumber = Integer.parseInt(tk.nextToken());
            
            // filename
            filename = tk.nextToken();
            
            // parent mass
            currentPm = Float.parseFloat(tk.nextToken());
            
            // charge
            currentCharge = Integer.parseInt(tk.nextToken());
            
            // activation method. For backward compatibility we check for next token
            if (tk.hasMoreTokens())    actMethod = tk.nextToken();
            
            // identification. For backward compatibility we check for next token
            if (tk.hasMoreTokens())    ident = tk.nextToken();
          }
          else {
            
            // parse the array
            String trimmed = line.trim();
            StringTokenizer tk = new StringTokenizer(trimmed.substring(1, trimmed.length()-1), ", ");
            ArrayList<Integer> sequence = new ArrayList<Integer>();
            while (tk.hasMoreTokens()) {
              int mass = Integer.parseInt(tk.nextToken());
              sequence.add(mass);
            }
            
            if (sequence.size() >= minLength) {
              // add the Spectrum if not added already
              if (!added) {
                gpr.addSpectrum(currentSpecId, currentScanNumber, filename, currentPm, currentCharge, actMethod, ident);
                added = true;  // signal we have added this already
              }
              
              gpr.addSequence(currentSpecId, sequence);
              queriesCount++;
            }
            else {
              discarded++;
            }
            
          }
        } while ((line=fid.readLine())!=null);
      }
      catch (IOException ioe) {
        System.err.println(ioe);
        System.exit(-9);
      }
      
      // we are done parsing the file
      nextLine = null;
      if (discarded > 0) {
        System.out.printf("Discarded %d gapped peptides due to insufficient length\n", discarded);
      }
      return gpr;
    }

    @Override
    public void remove() {
      throw new RuntimeException("Cannot remove entry from the MSGDResultFileParser iterator. This is not supported");
    }
    
  }
  
  private String filePath;
  private int chunkSize;
  private int minLength;
  
  /**
   * Constructor taking the path of the file to parse and the maximum number of
   * queries parsed in each GappedPeptideResult object.
   * @param filePath the path the .grc file
   * @param chunkSize the number of gapped queries in each GappedPeptideResult object
   *                  when the iterator is called.
   */
  public MSGDResultFileParser(String filePath, int chunkSize) {
    this(filePath, chunkSize, -1);
  }
  
  /**
   * Constructor taking the path of the file to parse and the maximum number of
   * queries parsed in each GappedPeptideResult object.
   * @param filePath the path the .grc file
   * @param chunkSize the number of gapped queries in each GappedPeptideResult object
   *                  when the iterator is called.
   * @param minLength the minimum length of the queries to parse
   */
  public MSGDResultFileParser(String filePath, int chunkSize, int minLength) {
    this.filePath = filePath;
    this.chunkSize = chunkSize;
    this.minLength = minLength;
  }
  
  /**
   * Default constructor taking the path to the .grc file and using the predetermined
   * bundling count for queries per each GappedPeptideResult object.
   * @param filePath
   */
  public MSGDResultFileParser(String filePath) {
    this(filePath, Constants.MAX_QUERY_BUNDLING_COUNT);
  }
  
  @Override
  public Iterator<GappedPeptideResults> iterator() {
    return new MyIterator();
  }
  
}
