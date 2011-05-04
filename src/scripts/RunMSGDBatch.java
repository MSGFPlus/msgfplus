package scripts;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import msgap.NewMSGappedDictionary;
import msgap.Parameters;


/**
 * This is a easy way to MSGD multiple times in the computer, generating the .grc and
 * .spr files only. No searching is done.
 * @author jung
 *
 */
public class RunMSGDBatch {
  
  // default output name 
  private final static String outName = "output6"; 
  //private final static int maxGapMass = 500;
  private final static int minGapLength = 6;
  
  private List<Parameters> work;

  private RunMSGDBatch(String[] inputs) {
    //int minGapLength = 6;
    
    List<Parameters> work = new ArrayList<Parameters>();
    for (String input : inputs) {
      if (new File(input).isDirectory()) {
        // running an entire directory
        String[] specs = {String.format("InputFile %s", input), 
                          String.format("OutputFile %s/%s", input, outName),
                          //String.format("MaxGapMass %d", maxGapMass)};
                          //String.format("SpecProb %e", Constants.PROB_CUTOFF),
                          String.format("Delta %d", minGapLength)};
        work.add(new Parameters(specs));
      }
      else {
        // running a specific file
        System.out.println(input);
        String prefix = input.substring(0, input.lastIndexOf('.'));
        String[] specs = {String.format("InputFile %s", input), 
                          String.format("OutputFile %s", prefix),
                          //String.format("MaxGapMass %d", maxGapMass)};
                          //String.format("SpecProb %e", Constants.PROB_CUTOFF),
                          String.format("Delta %d", minGapLength)};
        work.add(new Parameters(specs));
      }
    }
    this.work = Collections.synchronizedList(work);
  }
  
  private void launch(int threads) {
    for (int i=0; i<threads; i++) {
      new Thread(new MSGDRunner()).start();
    }
  }
  
  private class MSGDRunner implements Runnable {
    @Override
    public void run() {
      while (!work.isEmpty()) {
        long time = System.currentTimeMillis();
        Parameters p = work.remove(0);
        System.out.printf("Running %s\n", p.getOutFileName());
        NewMSGappedDictionary.run(p, false); 
        int rt = (int)((System.currentTimeMillis()-time)/1000.0);
        System.out.printf("Done running MSGD for %s in %d seconds\n", p.getOutFileName(), rt);
      }
    }
  }
  
  
  
  public static void runBatch() {
    // number of threads to launch
    final int threads = 1;
    
    String userHome = System.getProperty("user.home");
    
    // give the accurate parent mass priority
    /*
    String[] inputs = {String.format("%s/Data/Spectra/Asp/ORG014_LTQ_FT_0", userHome),
                       String.format("%s/Data/Spectra/Csp/ORG033_LTQ_Orb_0", userHome),
                       String.format("%s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_0", userHome),
                       String.format("%s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_1", userHome),
                       String.format("%s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_2", userHome),
                       String.format("%s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_3", userHome),
                       String.format("%s/Data/Spectra/Scerv/ORG105_LTQ_FT_0", userHome),
                       String.format("%s/Data/Spectra/Scerv/ORG105_LTQ_Orb_0", userHome),
                       String.format("%s/Data/Spectra/Asp/ORG014_LTQ_0", userHome),
                       String.format("%s/Data/Spectra/Asp/ORG014_LTQ_1", userHome),
                       String.format("%s/Data/Spectra/Asp/ORG014_LTQ_2", userHome),
                       String.format("%s/Data/Spectra/Asp/ORG014_LTQ_3", userHome),
                       String.format("%s/Data/Spectra/Csp/ORG033_LTQ_0", userHome),
                       String.format("%s/Data/Spectra/Avar/ORG013_LTQ_0", userHome),
                       String.format("%s/Data/Spectra/Avar/ORG013_LTQ_1", userHome),
                       String.format("%s/Data/Spectra/Ecoli/ORG045_LTQ_0", userHome),
                       String.format("%s/Data/Spectra/Scerv/ORG105_LTQ_0", userHome),
                       String.format("%s/Data/Spectra/Scerv/ORG105_LCQ_0", userHome)
                      };
                      */
    String[] inputs = {String.format("%s/Data/Spectra/Sone/LTQFT0", userHome),
                       String.format("%s/Data/Spectra/Sone/LTQFT1", userHome),
                       String.format("%s/Data/Spectra/Sone/LTQFT2", userHome),
                       String.format("%s/Data/Spectra/Sone/LTQFT3", userHome),
                       String.format("%s/Data/Spectra/Sone/LTQFT4", userHome),
                       String.format("%s/Data/Spectra/Sone/LTQFT5", userHome)
                      };
    new RunMSGDBatch(inputs).launch(threads);
  }
  
  public static void runSingle() {
    // number of threads to launch
    final int threads = 3;
    
    String userHome = System.getProperty("user.home");
    
    /*
    String[] inputs = {String.format("%s/Data/Spectra/Hsapiens/Heck/mzXML/lys", userHome),
                       String.format("%s/Data/Spectra/Hsapiens/Heck/mzXML/tryp", userHome)};
    */
    
    String[] inputs = {String.format("%s/Data/Spectra/Maize/SQS02/embryo", userHome),
                       String.format("%s/Data/Spectra/Maize/SQS02/endosperm-10", userHome),
                       String.format("%s/Data/Spectra/Maize/SQS02/endosperm-12", userHome),
                       String.format("%s/Data/Spectra/Maize/SQS02/endosperm-8", userHome),
                       String.format("%s/Data/Spectra/Maize/SQS02/pericarp", userHome)};
    
    new RunMSGDBatch(inputs).launch(threads);
  }
  
  
  public static void main(String[] args) {
    //runBatch();
    runSingle();
  }

}
