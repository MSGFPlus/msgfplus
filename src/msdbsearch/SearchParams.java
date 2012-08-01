package msdbsearch;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import params.FileParameter;
import params.IntRangeParameter;
import params.ParamManager;
import params.ToleranceParameter;

import msgf.Tolerance;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.InstrumentType;
import msutil.Protocol;
import msutil.DBSearchIOFiles;
import msutil.SpecFileFormat;

public class SearchParams {
	private List<DBSearchIOFiles> dbSearchIOList;
	private File databaseFile;
	private Tolerance leftParentMassTolerance; 
	private Tolerance rightParentMassTolerance; 
	private int min13C;
	private int max13C;
	private Enzyme enzyme;
	private int numTolerableTermini;
	private ActivationMethod activationMethod;
	private InstrumentType instType;
	private Protocol protocol;
	private AminoAcidSet aaSet;
	private int numMatchesPerSpec;
	private int startSpecIndex;
	private int endSpecIndex;
	private boolean useTDA;
	private boolean showFDR;
	private boolean showDecoy;
	private int minPeptideLength;
	private int maxPeptideLength;
	private int minCharge;
	private int maxCharge;
	private int numThreads;
	private boolean replicateMergedResults;
	private boolean doNotDseEdgeScore;
	private File dbIndexDir;
	
	public SearchParams()	{}
	
	public List<DBSearchIOFiles> getDBSearchIOList() {
		return dbSearchIOList;
	}

	public File getDatabaseFile() {
		return databaseFile;
	}

	public Tolerance getLeftParentMassTolerance() {
		return leftParentMassTolerance;
	}

	public Tolerance getRightParentMassTolerance() {
		return rightParentMassTolerance;
	}

	public int getMin13C() {
		return min13C;
	}

	public int getMax13C() {
		return max13C;
	}
	
	public Enzyme getEnzyme() {
		return enzyme;
	}

	public int getNumTolerableTermini() {
		return numTolerableTermini;
	}

	public ActivationMethod getActivationMethod() {
		return activationMethod;
	}

	public InstrumentType getInstType() {
		return instType;
	}

	public Protocol getProtocol() {
		return protocol;
	}

	public AminoAcidSet getAASet() {
		return aaSet;
	}

	public int getNumMatchesPerSpec() {
		return numMatchesPerSpec;
	}

	public int getStartSpecIndex() {
		return startSpecIndex;
	}

	public int getEndSpecIndex() {
		return endSpecIndex;
	}

	public boolean useTDA() {
		return useTDA;
	}

	public boolean showFDR() {
		return showFDR;
	}

	public boolean showDecoy() {
		return showDecoy;
	}

	public int getMinPeptideLength() {
		return minPeptideLength;
	}

	public int getMaxPeptideLength() {
		return maxPeptideLength;
	}

	public int getMinCharge() {
		return minCharge;
	}

	public int getMaxCharge() {
		return maxCharge;
	}

	public int getNumThreads() {
		return numThreads;
	}

	public boolean replicateMergedResults() {
		return replicateMergedResults;
	}

	public boolean doNotDseEdgeScore() {
		return doNotDseEdgeScore;
	}

	public File getDBIndexDir() {
		return dbIndexDir;
	}
	
	public String parse(ParamManager paramManager)
	{
		// Spectrum file
		FileParameter specParam = paramManager.getSpecFileParam();
		File specPath = specParam.getFile();
		
		dbSearchIOList = new ArrayList<DBSearchIOFiles>();
		
		if(!specPath.isDirectory())
		{
			// Spectrum format
			SpecFileFormat specFormat = (SpecFileFormat)specParam.getFileFormat();
			// Output file
			File outputFile = paramManager.getOutputFileParam().getFile();
			
			dbSearchIOList = new ArrayList<DBSearchIOFiles>();
			dbSearchIOList.add(new DBSearchIOFiles(specPath, specFormat, outputFile));
		}
		else	// spectrum directory
		{
			dbSearchIOList = new ArrayList<DBSearchIOFiles>();
			for(File f : specPath.listFiles())
			{
				SpecFileFormat specFormat = SpecFileFormat.getSpecFileFormat(f.getName());
				if(specParam.isSupported(specFormat))
				{
					String outputFileName = f.getName().substring(0, f.getName().lastIndexOf('.'))+".tsv";
					File outputFile = new File(outputFileName);
					if(outputFile.exists())
						return outputFile.getPath() + " already exists!";
					dbSearchIOList.add(new DBSearchIOFiles(f, specFormat, outputFile));
				}
			}
			return null;
		}
		
		// DB file
		databaseFile = paramManager.getDBFileParam().getFile();
		
		// PM tolerance
		ToleranceParameter tol = ((ToleranceParameter)paramManager.getParameter("t"));
		leftParentMassTolerance = tol.getLeftTolerance();
		rightParentMassTolerance = tol.getRightTolerance();
		
		int toleranceUnit = paramManager.getIntValue("u");
		if(toleranceUnit != 2)
		{
			boolean isTolerancePPM;
			if(toleranceUnit == 0)
				isTolerancePPM = false;
			else 
				isTolerancePPM = true;
			leftParentMassTolerance = new Tolerance(leftParentMassTolerance.getValue(), isTolerancePPM);
			rightParentMassTolerance = new Tolerance(rightParentMassTolerance.getValue(), isTolerancePPM);
		}
		
		IntRangeParameter c13Param = (IntRangeParameter)paramManager.getParameter("c13");
		this.min13C = c13Param.getMin();
		this.max13C = c13Param.getMax();
		
		if(rightParentMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f || leftParentMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f)
			min13C = max13C = 0;
		
		enzyme = paramManager.getEnzyme();
		numTolerableTermini = paramManager.getIntValue("ntt");
		activationMethod = paramManager.getActivationMethod();
		instType = paramManager.getInstType();
		if(activationMethod == ActivationMethod.HCD)
			instType = InstrumentType.HIGH_RESOLUTION_LTQ;
		
		protocol = paramManager.getProtocol();
		
		aaSet = null;
		File modFile = paramManager.getModFileParam().getFile();
		if(modFile == null)
			aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		else
		{
			String modFileName = modFile.getName();
			String ext = modFileName.substring(modFileName.lastIndexOf('.')+1);
			if(ext.equalsIgnoreCase("xml"))
				aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(modFile.getPath());
			else
				aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
			if(aaSet.containsPhosphorylation())
			{
				protocol = Protocol.PHOSPHORYLATION;
			}
		}
		
		numMatchesPerSpec = paramManager.getIntValue("n");
		
		startSpecIndex = ((IntRangeParameter)paramManager.getParameter("index")).getMin();
		endSpecIndex = ((IntRangeParameter)paramManager.getParameter("index")).getMax();
		
		useTDA = paramManager.getIntValue("tda") == 1 ? true : false;
		showFDR = paramManager.getIntValue("showFDR") == 1 ? true : false;
		showDecoy = paramManager.getIntValue("showDecoy") == 1 ? true : false;
		
		minPeptideLength = paramManager.getIntValue("minLength");
		maxPeptideLength = paramManager.getIntValue("maxLength");
		if(minPeptideLength > maxPeptideLength)
		{
			return "MinPepLength must not be larger than MaxPepLength";
		}
		
		minCharge = paramManager.getIntValue("minCharge");
		maxCharge = paramManager.getIntValue("maxCharge");
		if(minCharge > maxCharge)
		{
			return "MinCharge must not be larger than MaxCharge";
		}
		
		numThreads = paramManager.getIntValue("thread");
		replicateMergedResults = paramManager.getIntValue("replicate") == 1 ? true : false;
		doNotDseEdgeScore = paramManager.getIntValue("edgeScore") == 1 ? true : false;
		
		dbIndexDir = paramManager.getFile("dd");
		
		return null;
	}
	
}
