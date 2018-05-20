package msgfplus;

import java.io.*;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import edu.ucsd.msjava.msdbsearch.CompactFastaSequence;
import edu.ucsd.msjava.msdbsearch.ReverseDB;
import edu.ucsd.msjava.ui.MSGFPlus;
import org.junit.Ignore;
import org.junit.Test;

import edu.ucsd.msjava.misc.ConvertToMgf;
import edu.ucsd.msjava.misc.VennDiagram;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScoredSpectrum;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.TSVParser;

public class TestMisc {

    @Test
    @Ignore
    public void testCleavageState() {

        Map<String, Integer> peptides = new HashMap<String, Integer>(){
            {
                // These test cases correspond to those in the UnitTests project of the
                // Peptide Hit Results Processor.  See:
                // https://github.com/PNNL-Comp-Mass-Spec/PHRP/blob/master/UnitTests/PeptideCleavageStateCalculatorTests.cs

                // Fully tryptic peptides
                put("K.ACDEFGR.S", 2); // Normal, fully tryptic peptide
                put("R.ACDEFGR.S", 2); // Normal, fully tryptic peptide
                put("-.ACDEFGR.S", 2); // Fully tryptic at the N-Terminus of the protein
                put("R.ACDEFGH.-", 1); // Fully tryptic at the C-Terminus of the protein; getNumCleavedTermini reports 1
                put("-.ACDEFG.-",  1); // Peptide spans the entire protein; getNumCleavedTermini reports 1

                // Partially tryptic peptides
                put("K.ACDEFGH.S", 1); // Normal, partially tryptic peptide
                put("L.ACDEFGR.S", 1); // Normal, partially tryptic peptide
                put("K.ACDEFGR.P", 2); // Would have been fully tryptic, but ends with R followed by P; getNumCleavedTermini reports 2
                put("K.PCDEFGR.S", 2); // Would have been fully tryptic, but starts with K followed by P; getNumCleavedTermini reports 2

                // Non-tryptic peptides
                put("L.ACDEFGH.S", 0); // Normal, non-tryptic peptide
                put("-.ACDEFGH.S", 1); // Normal, non-tryptic peptide that happens to be at the N-terminus; getNumCleavedTermini reports 1
                put("L.ACDEFGH.-", 0); // Normal, non-tryptic peptide that happens to be at the C-terminus
                put("L.ACDEFGR.P", 1); // Would have been partially tryptic, but ends with R followed by P; getNumCleavedTermini reports 1
                put("K.PCDEFGR.P", 2); // Would have been fully tryptic, but has a P after both the K and the R; getNumCleavedTermini reports 2
            }
        };

        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCysWithTerm();
        aaSet.registerEnzyme(Enzyme.TRYPSIN);

        Enzyme enzyme = Enzyme.getEnzymeByName("Tryp");

        for (Map.Entry<String, Integer> entry : peptides.entrySet()) {
            Integer computedTerminii = enzyme.getNumCleavedTermini(entry.getKey(), aaSet);

            Integer expectedterminii = entry.getValue();
            System.out.println("Peptide " + entry.getKey() + " has computedTerminii = " + computedTerminii + "; expected " + expectedterminii);
        }


    }

    @Test
    public void testMasses()
    {
        System.out.println(Composition.H - Composition.ChargeCarrierMass());
    }
    
    @Test
    public void testMisc()
    {
        String title = "Scan:25485 RT:62.983 PrecursorScan:25482 nMSN:19700 PrecursorMonoisoMZ:1134.2547 PEPMASS:Monoiso PrecursorMZ:1134.5891 PrecursorCharge:3 PrecursorScanFTMS:1 FTResolution:17500 IBP:3750861.85 ITot:126255817.49 max2med:29.41 InjTime:31.98 HCD=54.0063972473145eV IsolationMZ:1134.5900 PrecursorAb:8797869.00 MPY:1.00 ms1PrecursorTotAb:47857048671.30 ms1PrecursorInjTime:0.26 ms1PrecursorMZ:1134.5891 ms1PrecursorMzAvg:1134.8744 ms1PrecursorMzRMS:0.3998 ms1PrecursorIntens:8797869.00 ms1PrecursorRT:62.974 ms2IsolationWidth:2.50 ms1SelMZ:1134.2067-1135.5900 ms1SelAvgMZ:1134.7794 ms1SelRmsMZ:0.0636 PrecursorHasMax:1 ms1PrecursorAb:110982284.68 ms1PrecursorMax:111444926.50 numOCMF:270,270,0,7 PrecursorMaxMZ:1134.9262 PrecursorMaxAb:33904613.19 PrecursorMaxRT:62.988 PrecursorWayMMF:0.71 PrecursorMaxMMF:0.71 mzRmsMax:1.33 mzRmsMs2:1.32 maxDelMz:0.3310,5-0,99,99 ms2DelMz:0.3309,5-0,99,99 FilterMzPeakExists(25482):1 PCFD2,2;1,0,1134.2562,3,5,0.3331,0.0025,99,1,14,1.04,0.18,10,10,1.8,3.2,7.5,961;1,0,1134.5877,1,2,0.9994,0.0000,41,1,25,1.25,0.42,16,16,1.8,3.2,7.5,915 Precursor1HasMax:1 Precursor1MaxInjTime:0.26 Precursor1MaxTotAb:47857049600.00 Precursor1MaxAb:108699575.31 Precursor1MaxRT:62.974 Precursor1MaxWidth:0.982 Precursor1MaxWid50:0.174 Precursor1MaxRatio:1.0117 Precursor1MaxBkg:0.00 Precursor1AbuBkg:0.00 Precursor1MaxHW:0.50 Precursor1MaxSkew:-0.00 Precursor2HasMax:1 Precursor2MaxInjTime:0.26 Precursor2MaxTotAb:47857049600.00 Precursor2MaxAb:108699575.31 Precursor2MaxRT:62.974 Precursor2MaxWidth:0.960 Precursor2MaxWid50:0.174 Precursor2MaxRatio:1.0117 Precursor2MaxBkg:0.00 Precursor2AbuBkg:0.00 Precursor2MaxHW:0.50 Precursor2MaxSkew:-0.00 PrecursorMaxNoise:2.28 PrecursorRTStep:0.012 ConvVer:20120705a NumPeaks:472 Filter:FTMS + p NSI d Full ms2 1134.59@hcd28.00<mailto:1134.59@hcd28.00> [100.00-3495.00]";
        System.out.println(title.matches("^Scan:\\d+\\s.+"));
        System.out.println(title.matches("^Scan:\\d+\\sRT:\\d+\\.\\d+\\s.+"));
        System.out.println(title.matches("^Scan:\\d+\\sRT:\\d+\\.\\d+\\sPrecursorScan:\\d+\\??\\s.+"));
        String[] token = title.split("\\s+");
        int scanNum = Integer.parseInt(token[0].substring("Scan:".length()));
        System.out.println(scanNum);
    }
    
    @Test
    public void testVennDiagram()
    {
        File result1 = new File("/Users/kims336/Research/Data/Tao/Global/MSGFPlus_10ppm_TI1/CPTAC_OvC_JB5427_iTRAQ_01_9Apr12_Cougar_12-03-21_dta.tsv");
        File result2 = new File("/Users/kims336/Research/Data/Tao/Global/MSGFPlus_20ppm_TI2/CPTAC_OvC_JB5427_iTRAQ_01_9Apr12_Cougar_12-03-21_dta.tsv");
        
        try {
            VennDiagram.vennDiagram(result1, result2, 0.01f);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    @Test
    public void testMzMLParser() throws URISyntaxException, IOException {
        File mzMLFile = new File(TestMisc.class.getClassLoader().getResource("tiny.pwiz.mzML").toURI());
        File mgfFile = File.createTempFile("tiny.pwiz", "mgf");
        try {
            ConvertToMgf.convert(mzMLFile, mgfFile, false, null, null, -1, -1, -1, false);
        } catch (Exception e) {
            e.printStackTrace();
        }
        mgfFile.deleteOnExit();
    }
    
    @Test
    public void testTrypsinCredit()
    {
        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCysWithTerm();
        aaSet.registerEnzyme(Enzyme.TRYPSIN);
        System.out.println("PeptideCleavageCredit: " + aaSet.getPeptideCleavageCredit());
        System.out.println("PeptideCleavagePenalty: " + aaSet.getPeptideCleavagePenalty());
        System.out.println("NeighborCredit: " + aaSet.getNeighboringAACleavageCredit());
        System.out.println("NeighborPenalty: " + aaSet.getNeighboringAACleavagePenalty());
        
    }
    
    @Test
    @Ignore
    public void generateTRexPRMSpectrum()
    {
        File specFile = new File("D:\\Research\\Data\\TRex\\TRex_GLVGAPGLRGLPGK.mgf");
        SpectraAccessor accessor = new SpectraAccessor(specFile);
        Spectrum spec = accessor.getSpecItr().next();
        
        NewRankScorer scorer = NewScorerFactory.get(ActivationMethod.CID, InstrumentType.LOW_RESOLUTION_LTQ, Enzyme.TRYPSIN, Protocol.STANDARD);
        
        scorer.doNotUseError();
        NewScoredSpectrum<NominalMass> scoredSpec = scorer.getScoredSpectrum(spec);
        int maxNominalMass = NominalMass.toNominalMass(spec.getParentMass());
        
        // PRM spectrum
        System.out.println("BEGIN IONS");
        System.out.print("TITLE=PRM_SpecIndex="+spec.getSpecIndex());
        if(spec.getTitle() != null)
            System.out.println(" " + spec.getTitle());
        else
            System.out.println();
        if(spec.getAnnotation() != null)
            System.out.println("SEQ=" + spec.getAnnotationStr());
        System.out.println("PEPMASS=" + spec.getPrecursorPeak().getMz());
        System.out.println("SCANS=" + spec.getScanNum());
        System.out.println("CHARGE="+spec.getCharge()+"+");
        int peptideNominalMass = 1272;
        for(int m=1; m<maxNominalMass; m++)
        {
            NominalMass prm = new NominalMass(m);
            NominalMass srm = new NominalMass(peptideNominalMass-m);
            float prefixScore = scoredSpec.getNodeScore(prm, true);
            float suffixScore = scoredSpec.getNodeScore(srm, false);
            System.out.format("%d\t%d\n", m, Math.round(prefixScore+suffixScore));
            
        }
        System.out.println("END IONS");
    }        

    @Test
    @Ignore
    public void generateTRexPRMSpectra()
    {
        File outputFile = new File("D:\\Research\\Data\\TRex\\MaxCharge4\\TRex48216_Vectors.txt");
        PrintStream out = null;
        try {
            out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
        File idFile = new File("D:\\Research\\Data\\TRex\\MaxCharge4\\NoDecoy.tsv");
        HashMap<String, Integer> titleToNominalMass = new HashMap<String, Integer>();
        TSVParser parser = new TSVParser();
        parser.parse(idFile.getPath());
        ArrayList<String> titleList = parser.getList("Title");
        ArrayList<String> peptideList = parser.getList("Peptide");
        ArrayList<String> specEValueList = parser.getList("SpecEValue");
        for(int i=0; i<specEValueList.size(); i++)
        {
            double specEValue = Double.parseDouble(specEValueList.get(i));
            if(specEValue > 1E-10) continue;
            Peptide peptide = new Peptide(peptideList.get(i), aaSet);
            int nominalMass = peptide.getNominalMass();
            String title = titleList.get(i);
            titleToNominalMass.put(title, nominalMass);
        }
        
        NewRankScorer scorer = NewScorerFactory.get(ActivationMethod.CID, InstrumentType.LOW_RESOLUTION_LTQ, Enzyme.TRYPSIN, Protocol.STANDARD);
        scorer.doNotUseError();
        
        File specFile = new File("D:\\Research\\Data\\TRex\\TRex48216.mgf");
        SpectraAccessor accessor = new SpectraAccessor(specFile);
        Iterator<Spectrum> itr = accessor.getSpecItr();
        while(itr.hasNext())
        {
            Spectrum spec = accessor.getSpecItr().next();
            String title = spec.getTitle();
            int nominalMass;
            if(titleToNominalMass.containsKey(title)) nominalMass = titleToNominalMass.get(title);
            else nominalMass = NominalMass.toNominalMass(spec.getParentMass()) - 18;
            
            NewScoredSpectrum<NominalMass> scoredSpec = scorer.getScoredSpectrum(spec);
            
            // PRM spectrum
            //out.println("BEGIN IONS");
            out.println("SCAN="+spec.getScanNum());
//            if(spec.getTitle() != null)
//                out.println(" " + spec.getTitle());
//            else
//                out.println();
//            if(spec.getAnnotation() != null)
//                out.println("SEQ=" + spec.getAnnotationStr());
//            out.println("PEPMASS=" + spec.getPrecursorPeak().getMz());
            out.println("PEPTIDE_MASS=" + nominalMass);
//            out.println("SCANS=" + spec.getScanNum());
//            out.println("CHARGE="+spec.getCharge()+"+");
            
//            int peptideNominalMass = 1272;
            for(int m=1; m<nominalMass; m++)
            {
                NominalMass prm = new NominalMass(m);
                NominalMass srm = new NominalMass(nominalMass-m);
                float prefixScore = scoredSpec.getNodeScore(prm, true);
                float suffixScore = scoredSpec.getNodeScore(srm, false);
                out.println(m+"\t"+Math.round(prefixScore+suffixScore));
            }
            out.println(nominalMass+"\t"+0);
        }
        System.out.println("Done.");            
    }


    @Test
    public void testReverseDB() throws URISyntaxException, IOException {
        File dbFile = new File(TestSA.class.getClassLoader().getResource("ecoli.fasta").toURI());
        File dbDecoyFile = File.createTempFile("ecoli-reversed", ".fasta");
        ReverseDB.reverseDB(dbFile.getAbsolutePath(), dbDecoyFile.getAbsolutePath(), true, MSGFPlus.DECOY_PROTEIN_PREFIX);

        CompactFastaSequence tdaSequence = new CompactFastaSequence(dbDecoyFile.getPath());
        float ratioUniqueProteins = tdaSequence.getRatioUniqueProteins();
        if (ratioUniqueProteins < 0.5f) {
            tdaSequence.printTooManyDuplicateSequencesMessage(dbDecoyFile.getName(), "MS-GF+", ratioUniqueProteins);
        }

        dbDecoyFile.deleteOnExit();

    }
    
}
