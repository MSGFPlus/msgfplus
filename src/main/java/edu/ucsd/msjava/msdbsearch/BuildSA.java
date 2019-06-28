package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.ui.MSGFPlus;
import org.apache.commons.io.FilenameUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class BuildSA {

    /**
     * Constructor
     * @param argv
     */
    public static void main(String argv[]) {
        if (argv.length < 1)
            printUsageAndExit("");

        if (argv.length < 2 || argv.length % 2 != 0)
            printUsageAndExit("The number of parameters must be even. If a file path has a space, surround it with double quotes.");

        File dbPath = null;
        File outputDir = null;
        int mode = 2;
        String decoyProteinPrefix = MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX;

        for (int i = 0; i < argv.length; i += 2) {
            if (!argv[i].startsWith("-") || i + 1 >= argv.length)
                printUsageAndExit("Invalid parameters");
            else if (argv[i].equalsIgnoreCase("-d")) {
                dbPath = new File(argv[i + 1]);
                if (!dbPath.exists())
                    printUsageAndExit(argv[i + 1] + " doesn't exist.");
            } else if (argv[i].equalsIgnoreCase("-o")) {
                outputDir = new File(argv[i + 1]);
            } else if (argv[i].equalsIgnoreCase("-tda")) {
                if (argv[i + 1].equals("0"))
                    mode = 0;
                else if (argv[i + 1].equals("1"))
                    mode = 1;
                else if (argv[i + 1].equals("2"))
                    mode = 2;
                else
                    printUsageAndExit("Invalid parameter: -tda " + argv[i + 1]);
            } else if (argv[i].equalsIgnoreCase("-decoy")) {
                decoyProteinPrefix = argv[i + 1];
            }
        }
        if (dbPath == null)
            printUsageAndExit("Database must be specified!");

        buildSA(dbPath, outputDir, mode, decoyProteinPrefix);
    }

    /**
     * Show the syntax
     * @param message
     */
    public static void printUsageAndExit(String message) {
        System.out.println();
        if (!message.isEmpty()) {
            System.out.println("Error: " + message);
            System.out.println();
        }
        System.out.println("Usage: java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.msdbsearch.BuildSA");
        System.out.println("\t-d DatabaseFile (*.fasta or *.fa or *.faa; if a directory path, index all FASTA files)");
        System.out.println("\t[-tda 0/1/2] (0: Target database only, 1: Concatenated target-decoy database only, 2: Both (Default))");
        System.out.println("\t[-o OutputDir] (Directory to save index files; default is the same as the input file)");
        System.out.println("\t[-decoy DecoyPrefix] (Prefix for decoy protein names; default is " + MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX + ")");
        System.out.println();
        System.out.println("Documentation: https://github.com/MSGFPlus/msgfplus");

        System.exit(-1);
    }

    /**
     * Index a directory with several FASTA files, or the specified FASTA file
     * @param dbPath
     * @param outputDir
     * @param mode
     * @param decoyProteinPrefix
     */
    public static void buildSA(File dbPath, File outputDir, int mode, String decoyProteinPrefix) {
        if (dbPath.isDirectory()) {
            for (File f : dbPath.listFiles()) {
                if (isFastaFile(f.getName())) {
                    buildSAFiles(f, outputDir, mode, decoyProteinPrefix);
                }
            }
        } else {
            if (isFastaFile(dbPath.getName())) {
                buildSAFiles(dbPath, outputDir, mode, decoyProteinPrefix);
            }
        }
        System.out.println("Done");
    }

    /**
     * Index a protein database (FASTA file)
     * @param databaseFile       FASTA file path
     * @param outputDir          Output directory
     * @param mode               0: target only, 1: target-decoy only, 2: both
     * @param decoyProteinPrefix Decoy protein prefix
     */
    public static void buildSAFiles(File databaseFile, File outputDir, int mode, String decoyProteinPrefix) {
        if (outputDir == null) {
            outputDir = databaseFile.getAbsoluteFile().getParentFile();
        }

        if (!validateOutputDirectory(outputDir)) {
            System.exit(-1);
        }

        String dbFileName = databaseFile.getName();

        if (decoyProteinPrefix == null || decoyProteinPrefix.trim().isEmpty())
            decoyProteinPrefix = MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX;

        // Make sure that decoyProteinPrefix does not end in an underscore, since we add it below
        while (decoyProteinPrefix.endsWith("_")) {
            decoyProteinPrefix = decoyProteinPrefix.substring(0, decoyProteinPrefix.length() - 1);
        }

        if (decoyProteinPrefix.trim().isEmpty())
            decoyProteinPrefix = MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX;

        // decoy
        if (mode == 1 || mode == 2) {
            String concatDBFileName = dbFileName.substring(0, dbFileName.lastIndexOf('.')) + MSGFPlus.DECOY_DB_EXTENSION;
            File concatTargetDecoyDBFile = new File(Paths.get(outputDir.getPath(), concatDBFileName).toString());
            if (!concatTargetDecoyDBFile.exists()) {
                System.out.println("Creating " + concatDBFileName + ".");
                if (!ReverseDB.reverseDB(databaseFile.getPath(), concatTargetDecoyDBFile.getPath(), true, decoyProteinPrefix)) {
                    System.err.println("Cannot create decoy database file!");
                    System.out.println("Consider using -o to specify the output directory");
                    System.exit(-1);
                }
            }
            System.out.println("Building suffix array: " + concatTargetDecoyDBFile.getPath());
            CompactFastaSequence tdaSequence = new CompactFastaSequence(concatTargetDecoyDBFile.getPath());
            tdaSequence.setDecoyProteinPrefix(decoyProteinPrefix);

            float ratioUniqueProteins = tdaSequence.getRatioUniqueProteins();
            if (ratioUniqueProteins < 0.5f) {
                tdaSequence.printTooManyDuplicateSequencesMessage(concatTargetDecoyDBFile.getName(), "MS-GF+", ratioUniqueProteins);
                System.exit(-1);
            }

            float fractionDecoyProteins = tdaSequence.getFractionDecoyProteins();
            if (fractionDecoyProteins < 0.4f || fractionDecoyProteins > 0.6f) {
                System.err.println("Error while reading: " + databaseFile.getName() + " (fraction of decoy proteins: " + fractionDecoyProteins + ")");
                if (databaseFile.getName().toLowerCase().endsWith(".revCat.fasta".toLowerCase())) {
                    System.err.println("Delete " + databaseFile.getName() + " and run MS-GF+ (or BuildSA) again.");
                } else {
                    String baseName = FilenameUtils.removeExtension(databaseFile.getName());
                    System.err.println("Delete files starting with " + baseName +
                            " (but keep " + databaseFile.getName() + ") and run MS-GF+ (or BuildSA) again.");
                }
                System.err.println("Decoy protein names should start with " + tdaSequence.getDecoyProteinPrefix());
                System.exit(-1);
            }

            new CompactSuffixArray(tdaSequence);
        }

        if (mode == 0 || mode == 2) {
            File targetDBFile = new File(Paths.get(outputDir.getPath(), dbFileName).toString());
            if (!targetDBFile.exists()) {
                System.out.println("Creating " + targetDBFile.getName() + ".");
                if (!ReverseDB.copyDB(databaseFile.getPath(), targetDBFile.getPath())) {
                    System.err.println("Cannot create target database file!");
                    System.out.println("Consider using -o to specify the output directory");
                    System.exit(-1);
                }
            }
            System.out.println("Building suffix array: " + databaseFile.getPath());
            CompactFastaSequence sequence = new CompactFastaSequence(targetDBFile.getPath());
            sequence.setDecoyProteinPrefix(decoyProteinPrefix);

            new CompactSuffixArray(sequence);
        }

        System.out.println();
    }

    /**
     * Return True if the file path ends in .fasta, .fa, or .faa
     * @param filePath
     * @return
     */
    public static boolean isFastaFile(String filePath) {
        String fileNameLcase = filePath.toLowerCase();

        return fileNameLcase.endsWith(".fasta") ||
               fileNameLcase.endsWith(".fa") ||
               fileNameLcase.endsWith(".faa");
    }

    private static boolean validateOutputDirectory(File outputDir) {

        try {
            if (!outputDir.exists()) {
                // Attempt to create the output directory
                Boolean success = outputDir.mkdirs();
                if (!success) {
                    System.err.println("Error creating the output directory (access denied?): " + outputDir.getPath());
                    return false;
                }
            }
        }
        catch (Throwable ex) {
            System.err.println("Error validating / creating the output directory: " + outputDir.getPath());
            return false;
        }

        // Assure that we can create files in the output directory
        Path testFilePath = Paths.get(outputDir.getPath(), "WritePermTestFile.tmp");

        if (!Files.isWritable(testFilePath)) {

            Boolean accessDenied = true;

            try {
                // On Windows 10, Files.isWritable() returns false on a newly created directory where we _do_ have write permission
                // Try creating a test file

                File testFile = new File(testFilePath.toString());
                if (testFile.exists())
                    testFile.delete();

                BufferedWriter writer = new BufferedWriter(new FileWriter(testFile.getPath()));
                writer.write("test");
                writer.close();

                if (testFile.exists()) {
                    // Files.isWritable reports false, but we were able to create a test file
                    accessDenied = false;
                    testFile.delete();
                }

            } catch (Exception ex) {
                // Ignore exceptions here
            }

            if (accessDenied) {
                System.err.println("Write access denied to directory: " + outputDir.getPath());
                System.out.println("Consider using -o to specify the output directory");
                return false;
            }
        }

        return true;
    }

}
