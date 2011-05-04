package msgap;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import msgf.Tolerance;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.Spectrum;
import parser.BufferedLineReader;

public class Parameters {
  
  private static final String DATE_FORMAT_NOW = "yyyy-MM-dd-HH-mm";
  public static String getNowTimeStr() {
    return new SimpleDateFormat(DATE_FORMAT_NOW).format(Calendar.getInstance().getTime());
  }
  
	static class SpecFileFilter implements FilenameFilter{
		String[] ext = { ".mgf", ".mzxml", "ms2"};
		
		public boolean accept(File dir, String name) {
			name = name.toLowerCase();
			for(String extension : ext){
				if(name.endsWith(extension)) return true;
			}
			return false;
		} 
	}
	
	
	
	// this is the root directory used to store the temp files
	private String outputRootDir;
	private String sprPath;
	private String grcPath;
	
	   private int numHubs = 20; 
	   private int dictionarySize = 60;
	   private int delta = 5;
	   private int maxGapMass = Integer.MAX_VALUE;
	   
	   private boolean allowNonEnzymaticCleavage = true;
	   private float specProb = 1e-9f;
	   //   Tolerance tolerance = new Tolerance(0.5f, false);
	   private String[] specFiles = null;
	   private String outFileName  = null;
	   private String dbFileName = null;
	 //  private String decoydbFileName = null;
	   private String aaFileName = null;
	   private int minSpecCharge = 0;	
	   private int maxSpecCharge = 100;
	   private int msgfScoreThreshold = 0;
	   private int matchScoreThreshold = 0;
	   private Tolerance ionTolerance = new Tolerance(0.5f, false);
	   private Tolerance pmTolerance = new Tolerance(2.0f, false);
	   private boolean isETD = false;
	   private Enzyme enzyme = Enzyme.TRYPSIN;
	   private int maxModNum = 0;
	   private int maxModNumPerOneModification = 0;
	   private boolean useGeneratedGrcFile = false;
	   private AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
	   boolean correctPM = false;
	   private String msgfParaFile = null;
	   
	   String parseParameters(String s){
		   String ret = "";
		   
		   String[] token = s.split(" ");
		   String para = token[0];
		   String paraValue = token[token.length-1];
		                
		   if(para.equals("InputFile") || para.equals("-s")){
			   File a = new File(paraValue);
			   File[] specFiles;
			   
			   if(!a.exists()) return s + " : input spectrum file does not exist";
			   
			   if(a.isDirectory()){
				   specFiles = a.listFiles(new SpecFileFilter());
			   }else{
				   specFiles = new File[1];
				   specFiles[0] = a;
			   }
			  // if(specFiles == null) ret =  s + " [input spectra file] needs to be specified)";
			   
			   this.specFiles = new String[specFiles.length];
			   for(int i=0; i< specFiles.length; i++)
				   this.specFiles[i] = specFiles[i].getAbsolutePath();
		   }
		   else if(para.equals("OutputFile") || para.equals("-o")){
			   outFileName = paraValue;
			 //  if(outFileName == null) ret =  s + " [output file prefix] needs to be specified";
		   }
		   else if(para.equals("DBFile") || para.equals("-d")){
			   dbFileName = paraValue;
			   File dbfile =  new File(dbFileName);
			   if(! dbfile.exists()) return s + " : data base file does not exist";
			   if(dbfile.isDirectory());
			   else if(!dbFileName.toLowerCase().endsWith("fasta") && !dbFileName.toLowerCase().endsWith("fa")) 
				   return s + " : MS-GappedDictionary accepts only *.fasta file";
		   }
	//	   else if(para.equals("DecoyDBFile")){
		//	   decoydbFileName = paraValue;
	//	   }   
		   else if(para.equals("SpecCharge") || para.equals("-c")){
			   String[] t = paraValue.split(":");
			   
			   if(t[0].isEmpty())	minSpecCharge = 1;
			   else
				   maxSpecCharge = minSpecCharge = Integer.parseInt(t[0]);
			   if(t.length == 2){
				   if(!t[1].isEmpty()) 
					   maxSpecCharge = Integer.parseInt(t[1]);
				   else 
					   maxSpecCharge = 100;
			   }
		   }
		   else if(para.equals("SpecProb") || para.equals("-p")){
			   specProb = Float.parseFloat(paraValue);
		   }
		   else if(para.equals("NumHubs") || para.equals("-h")){ // hidden
			   numHubs = Integer.parseInt(paraValue);
		   }
		   else if(para.equals("DictSize") || para.equals("-ds")){ // hidden
			   dictionarySize = Integer.parseInt(paraValue);
		   }   
		   else if(para.equals("Delta") || para.equals("-l")){
			   delta = Integer.parseInt(paraValue);
			   if(delta <=5) dictionarySize = 60;
			   else if(delta ==6) dictionarySize = 80;
			   else dictionarySize = 100;
		   }
		   else if(para.equals("MaxGapMass")){ // hidden
			   maxGapMass = Integer.parseInt(paraValue);
		   }
		   else if(para.equals("AllowMissedCleavage") || para.equals("-amc")){//hidden
			   allowNonEnzymaticCleavage = Boolean.parseBoolean(paraValue);
		   }
		  // else if(para.equals("IsSpecProbExclusive")){
			//   isSpecProbExclusive = Boolean.parseBoolean(paraValue);
		  // }
		   else if(para.equals("MSGFScore") || para.equals("-filter")){
			   msgfScoreThreshold = Integer.parseInt(paraValue);
		   }
		   else if(para.equals("AminoAcidFile") || para.equals("-aa")){ // hidden for now
			   aaFileName = paraValue;

			   aaSet = AminoAcidSet.getAminoAcidSet(aaFileName);
		   }
		   else if(para.equals("MatchScore") || para.equals("-ps")){
			   matchScoreThreshold = Integer.parseInt(paraValue);
		   }
		   else if(para.equals("IonTolerance") || para.equals("-it")){//hidden and no effect
			   if(paraValue.toUpperCase().endsWith("PPM"))
				   ionTolerance = new Tolerance(Float.parseFloat(paraValue.substring(0, paraValue.toUpperCase().lastIndexOf("PPM"))), true);
			   else if(paraValue.toUpperCase().endsWith("DA"))
				   ionTolerance = new Tolerance(Float.parseFloat(paraValue.substring(0, paraValue.toUpperCase().lastIndexOf("DA"))), false);
			   else
				   ionTolerance = new Tolerance(Float.parseFloat(paraValue), false);
		   }
		   else if(para.equals("PMTolerance") || para.equals("-t")){
			   if(paraValue.toUpperCase().endsWith("PPM"))
				   pmTolerance = new Tolerance(Float.parseFloat(paraValue.substring(0, paraValue.toUpperCase().lastIndexOf("PPM"))), true);
			   else if(paraValue.toUpperCase().endsWith("DA"))
				   pmTolerance = new Tolerance(Float.parseFloat(paraValue.substring(0, paraValue.toUpperCase().lastIndexOf("DA"))), false);
			   else
				   pmTolerance = new Tolerance(Float.parseFloat(paraValue), false);
		   }
		   else if(para.equals("Enzyme") || para.equals("-e")){
			   int id = Integer.parseInt(paraValue);
			   if(id == 0)
				   enzyme = null;
			   else if(id == 1)
				   enzyme = Enzyme.TRYPSIN;
			   else if(id == 2)
				   enzyme = Enzyme.CHYMOTRYPSIN;
			   else if(id == 3)
				   enzyme = Enzyme.LysC;
			   else if(id == 4)
				   enzyme = Enzyme.LysN;
			   else if(id == 5)
				   enzyme = Enzyme.GluC;
			   else if(id == 6)
				   enzyme = Enzyme.ArgC;
			   else if(id == 7)
				   enzyme = Enzyme.AspN;
			   else return s + " : Wrong enzyme ID";
				   
			   /*
			   if(paraValue.toUpperCase().equals("TRYPSIN"))
				   enzyme = Enzyme.TRYPSIN;
			   else if(paraValue.toUpperCase().equals("CHYMOTRYPSIN"))
				   enzyme = Enzyme.CHYMOTRYPSIN;
			   else if(paraValue.toUpperCase().equals("LYSC"))
				   enzyme = Enzyme.LysC;
			   else if(paraValue.toUpperCase().equals("LYSN"))
				   enzyme = Enzyme.LysN;
			   else if(paraValue.toUpperCase().equals("GLUC"))
				   enzyme = Enzyme.GluC;
			   else if(paraValue.toUpperCase().equals("ARGC"))
				   enzyme = Enzyme.ArgC;
			   else if(paraValue.toUpperCase().equals("ASPN"))
				   enzyme = Enzyme.AspN;
				   */
		   }
		   else if(para.equals("isETD") || para.equals("-m")){
			   int id = Integer.parseInt(paraValue);
			   if(id == 1)
				   isETD = false;
			   else if(id == 2)
				   isETD = true;
			   
			   else ret = s + " : Wrong fragmentation method ID";
			   		
			  // isETD = Boolean.parseBoolean(paraValue);
		   }
		   else if(para.equals("NumMod") || para.equals("-nm")){ // hidden
			   maxModNum = Integer.parseInt(paraValue);
			   if(maxModNumPerOneModification == 0) maxModNumPerOneModification = maxModNum;
		   }
		   else if(para.equals("NumMod2") || para.equals("-nm2")){// hidden
			   maxModNumPerOneModification = Integer.parseInt(paraValue);
		   }
		   else if(para.equals("UsePrevGRCFile")|| para.equals("-u")){
			   useGeneratedGrcFile = true;
		   }else if(para.equals("-cpm")){ // hidden
			   correctPM = true;
		   }else if(para.equals("-fixMod")){
			   int id = Integer.parseInt(paraValue);
			   if(id == 0)
				   aaSet = AminoAcidSet.getStandardAminoAcidSet();
			   else if(id == 1);
			   else if(id == 2)
				   aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarboxymethylatedCys();
			   else ret = s + " : Wrong fixMod ID";
		   }else if(para.equals("-para")){ // hidden
			   msgfParaFile = paraValue;
		   }
		   else{
			   ret = para;
		   }
		   
		   return ret;
	   }
	   
   
  	/**
  	 * Alternative way to build a parameter file by passing in the list of
  	 * parameters as strings instead of reading from a file
  	 * @param specs the specifications of the parameters, one per line
  	 */
    public Parameters(String[] specs) {
      for (String s : specs) { 
        s = s.replace('"', ' ').trim();
        if(s.startsWith("#") || s.isEmpty()) continue;
        if(s.contains("#")) s=s.substring(0, s.indexOf("#"));
        String ret;
        ret = parseParameters(s);
        if(!ret.isEmpty()){
          System.out.println("Invalid parameter: " + ret);
          System.exit(-1);
        }
      }
      
      if(specFiles == null){
  		System.out.println("-s [spectrum file] option missing");
      	System.exit(-1);
  	}else if(outFileName == null){
  		System.out.println("-o [output file prefix] option missing");
  	}
         
      // initialize the outputRootDirectory based on the output file
      File outfile = new File(getOutFileName());
      this.outputRootDir = outfile.getAbsoluteFile().getParent();
         
      this.grcPath = new File(this.outputRootDir, new File(getOutFileName()).getName() + ".grc").getAbsolutePath();
      this.sprPath = new File(this.outputRootDir, new File(getOutFileName()).getName() + ".spr").getAbsolutePath();
    }

    // parameters from command line. int numspecs is dummy
    public Parameters(String[] specs, int numspecs) {
    	String s="";
    	for(String spec : specs){
    		if(spec.startsWith("-")){
    			if(s.startsWith("-")){
    				String ret = parseParameters(s);
    		        if(!ret.isEmpty()){
    		        	System.out.println("Invalid parameter: " + ret);
    		        	System.exit(-1);
    		        }
    			}
    			s = spec;
    		}
    		else{
    			s+=" " + spec;
    			String ret = parseParameters(s);
		        if(!ret.isEmpty()){
		        	System.out.println("Invalid parameter: " + ret);
		        	System.exit(-1);
		        }
		        s="";
    		}
    	}
    	
    	if(s.startsWith("-")){
			String ret = parseParameters(s);
	        if(!ret.isEmpty()){
	        	System.out.println("Invalid parameter: " + ret);
	        	System.exit(-1);
	        }
		}
    	
    	if(specFiles == null){
    		System.out.println("-i [spectrum file] option missing");
        	System.exit(-1);
    	}else if(outFileName == null){
    		System.out.println("-o [output file prefix] option missing");
    	}
    	
    	
    	
    	  // initialize the outputRootDirectory based on the output file
        File outfile = new File(getOutFileName());
        this.outputRootDir = outfile.getAbsoluteFile().getParent();
           
        this.grcPath = new File(this.outputRootDir, new File(getOutFileName()).getName() + ".grc").getAbsolutePath();
        this.sprPath = new File(this.outputRootDir, new File(getOutFileName()).getName() + ".spr").getAbsolutePath();
   
    	
    }
  
  
  public Parameters(String confFileName, boolean verbose) {
    try {
      BufferedLineReader in =  new BufferedLineReader(confFileName);
      String s; 
      while((s=in.readLine()) != null){
        s = s.replace('"', ' ').trim();
        if(s.startsWith("#") || s.isEmpty()) continue;
        if(s.contains("#")) s=s.substring(0, s.indexOf("#"));
        String ret = parseParameters(s);
        if(!ret.isEmpty()){
          System.out.println("Invalid parameter: " + ret);
          System.exit(-1);
        }
        if(verbose) System.out.println(s);
      }
      
      if(specFiles == null){
  		System.out.println("-s [spectrum file] option missing");
      	System.exit(-1);
  	}else if(outFileName == null){
  		System.out.println("-o [output file prefix] option missing");
  	}
      
      // initialize the outputRootDirectory based on the output file
      File outfile = new File(getOutFileName());
      this.outputRootDir = outfile.getAbsoluteFile().getParent();
      
      //this.grcPath = new File(this.outputRootDir, new File(outFileName()).getName() + "_" + getNowTimeStr()+".grc").getAbsolutePath();
      //this.sprPath = new File(this.outputRootDir, new File(outFileName()).getName() +"_"+ getNowTimeStr()+".spr").getAbsolutePath();
      
      this.grcPath = new File(this.outputRootDir, new File(getOutFileName()).getName() + ".grc").getAbsolutePath();
      this.sprPath = new File(this.outputRootDir, new File(getOutFileName()).getName() + ".spr").getAbsolutePath();
      
    }
    catch (IOException e) {
      System.err.println("Problem reading the parameter file " + confFileName);
      System.err.println(e.getStackTrace());
      System.exit(-1);
    }
  }
	   
	  
	   public String getMSGFParaFile() {return msgfParaFile;}
	   public int numHubs() {return numHubs;}
	   public int dictionarySize() {return dictionarySize;}
	   public int delta() {return delta;}
	   public boolean allowNonEnzymaticCleavage() {return allowNonEnzymaticCleavage;}
	   public float getSpecProb() {return specProb;}
	   public String[] specFiles() {return specFiles;}
	   
	 //  public String rankScorerFileName() {return rankScorerFileName;}
	   public String dbFileName() {return dbFileName;}
	   public int minSpecCharge() {return minSpecCharge;}
	   public int maxSpecCharge() {return maxSpecCharge;}
	   //public boolean isSpecProbExclusive() {return isSpecProbExclusive;}
	   public int msgfScore() {return msgfScoreThreshold;}
	   public int matchScore() {return matchScoreThreshold;}
	   public AminoAcidSet aaSet() {
		   return aaSet;
	   }
	   public Tolerance ionTolerance() {return ionTolerance;}
	   public Tolerance pmTolerance() {return pmTolerance;}
	   public Enzyme enzyme() {return enzyme;}
	   //public boolean isETD() {return isETD;}
	   public int maxModNum() {return maxModNum;}
	   public int maxModNumPerOneModification() {return maxModNumPerOneModification;}
	   public boolean useGeneratedGrcFile() {return useGeneratedGrcFile;}
	 //  public String decoydbFileName() {return decoydbFileName;}
	   public int maxGapMass() {return maxGapMass;}
	   public float specProbThresholdBeforeConsideringFlankingAAs() {return specProb * aaSet.size()/2;}
	   
	   public ActivationMethod getActivationMethod(Spectrum s){
		  
			ActivationMethod method =  isETD? ActivationMethod.ETD : ActivationMethod.CID;
			
			  if(s.getActivationMethod() != null)
			      method = s.getActivationMethod();
			  
			return method;	
		}
	   
	public String getSPRPath()      { return sprPath; }
	public String getGRCPath()      { return grcPath; }
	public String getDBPath()       { return this.dbFileName; }
	public String getOutFileName()  { return outFileName; }
	
}
