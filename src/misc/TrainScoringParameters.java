package misc;

import java.io.File;
import java.util.Calendar;

import msscorer.ScoringParameterGeneratorWithErrors;
import msutil.AminoAcidSet;

public class TrainScoringParameters {
	
	private static final String PARAM_DIR = System.getProperty("user.home")+"/Developments/MS_Java_Dev/src/resources/ionstat";
	private static final String BACKUP_DIR = System.getProperty("user.home")+"/Developments/ionstat";
	private static final String SPEC_DIR = System.getProperty("user.home")+"/Research/Data/IonStat/SpectraForTraining";
	
	public static void main(String argv[]) throws Exception
	{
		backup();
		createParamFiles();
	}
	
	public static void backup() throws Exception
	{
		File paramDir = new File(PARAM_DIR);
		boolean paramExists = false;
		for(File paramFile : paramDir.listFiles())
		{
			if(paramFile.getName().endsWith(".param"))
				paramExists = true;
		}
		if(!paramExists)
		{
			System.out.println("No param file to backup.");
			return;
		}
		Calendar calendar = Calendar.getInstance();
		String dateStr = calendar.get(Calendar.MONTH) + "_" + calendar.get(Calendar.DAY_OF_MONTH) + "_" + calendar.get(Calendar.YEAR);
		String backupDirName = "ParamBackup_" + dateStr;
		File backupDir = new File(BACKUP_DIR+"/"+backupDirName);
		if(backupDir.exists())
		{
			System.out.println("Backup directory already exists: " + backupDir.getPath());
			System.exit(-1);
		}
		backupDir.mkdir();
		System.out.println(backupDir.getPath()+ " is created.");
		
		boolean backupSuccess = true;
		for(File paramFile : paramDir.listFiles())
		{
			if(paramFile.getName().endsWith(".param"))
			{
				File newFile = new File(backupDir, paramFile.getName());
				boolean isBackupSuccessful = paramFile.renameTo(newFile);
				System.out.println("Moving " + paramFile.getPath() + " to " + newFile.getPath() + (isBackupSuccessful ? " succeeded." : " failed."));
				if(!isBackupSuccessful)
				{
					backupSuccess = false;
					break;
				}
			}
		}
		if(backupSuccess)
			System.out.println("Backup complete.");
		else
		{
			backupDir.delete();
			System.out.println(backupDir.getPath()+ " is deleted.");
			System.out.println("Backup failed.");
			System.exit(0);
		}
	}
	
	public static void createParamFiles() throws Exception
	{
		File specDir = new File(SPEC_DIR);
		if(!specDir.exists())
		{
			System.err.println("Training spectra directory doesn't exist:" + specDir.getPath());
			System.exit(-1);
		}
		
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		for(File specFile : specDir.listFiles())
		{
			String specFileName = specFile.getName();
			if(specFileName.endsWith(".mgf"))
			{
				String id = specFileName.substring(0, specFileName.lastIndexOf('.'));
				String paramFileName = id+".param";
				File outputFile = new File(PARAM_DIR, paramFileName);
				System.out.println("Generating " + outputFile.getPath());
				int numSpecsPerPeptide = 1;
				int errorScalingFactor = 10;
				if(id.contains("HCD"))
				{
					numSpecsPerPeptide = 3;
					errorScalingFactor = 100;
				}
				else if(id.contains("ETD"))
				{
					errorScalingFactor = 0;
				}
				ScoringParameterGeneratorWithErrors.generateParameters(specFile, numSpecsPerPeptide, errorScalingFactor, outputFile, aaSet, false, false);
			}
		}
		System.out.println("Success!");
	}
}
