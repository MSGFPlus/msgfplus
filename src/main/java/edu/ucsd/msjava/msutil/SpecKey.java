package edu.ucsd.msjava.msutil;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import edu.ucsd.msjava.parser.MzXMLSpectraIterator;

public class SpecKey extends Pair<Integer, Integer> {

    private ArrayList<Integer> specIndexList;
    private float precursorMz;

    public SpecKey(int specIndex, int charge) {
        super(specIndex, charge);
    }

    public void setPrecursorMz(float precursorMz) {
        this.precursorMz = precursorMz;
    }

    public int getSpecIndex() {
        return super.getFirst();
    }

    public int getCharge() {
        return super.getSecond();
    }

    public float getPrecursorMz() {
        return precursorMz;
    }

    public String getSpecKeyString() {
        return getSpecIndex() + ":" + getCharge();
    }

    public static SpecKey getSpecKey(String specKeyString) {
        String[] token = specKeyString.split(":");
        return new SpecKey(Integer.parseInt(token[0]), Integer.parseInt(token[1]));
    }

    public void addSpecIndex(int scanNum) {
        if (specIndexList == null) {
            specIndexList = new ArrayList<Integer>();
//			specIndexList.add(super.getFirst());
        }
        specIndexList.add(scanNum);
    }

    public ArrayList<Integer> getSpecIndexList() {
        return specIndexList;
    }

    public static ArrayList<SpecKey> getSpecKeyList(Iterator<Spectrum> itr, int startSpecIndex,
                                                    int endSpecIndex, int minCharge, int maxCharge, ActivationMethod activationMethod, int minNumPeaksPerSpectrum) {
        if (activationMethod == ActivationMethod.FUSION)
            return getFusedSpecKeyList(itr, startSpecIndex, endSpecIndex, minCharge, maxCharge);

        ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();

        int numProfileSpectra = 0;
        int numSpectraWithTooFewPeaks = 0;
        final int MAX_INFORMATIVE_MESSAGES = 10;
        int informativeMessageCount = 0;

        while (itr.hasNext()) {
            Spectrum spec = itr.next();
            int specIndex = spec.getSpecIndex();

            if (specIndex < startSpecIndex)
                continue;
            if (specIndex >= endSpecIndex)
                continue;

            spec.setChargeIfSinglyCharged();
            int charge = spec.getCharge();

            if (activationMethod != ActivationMethod.ASWRITTEN) {
                if (informativeMessageCount < MAX_INFORMATIVE_MESSAGES) {
                    System.out.println(
                            "Use spectrum " + spec.getID() +
                                    " since assumed activationMethod is " + activationMethod.toString());
                    informativeMessageCount++;
                } else {
                    if (informativeMessageCount == MAX_INFORMATIVE_MESSAGES) {
                        System.out.println(" ...");
                        informativeMessageCount++;
                    }
                }
            } else {

                ActivationMethod currentActivationMethod = spec.getActivationMethod();
                if (currentActivationMethod == null) {
                    if (informativeMessageCount < MAX_INFORMATIVE_MESSAGES) {
                        System.out.println("Skip spectrum " + spec.getID() + " since activationMethod is unknown");
                        informativeMessageCount++;
                    } else {
                        if (informativeMessageCount == MAX_INFORMATIVE_MESSAGES) {
                            System.out.println(" ...");
                            informativeMessageCount++;
                        }
                    }
                    continue;
                }

                if (currentActivationMethod != activationMethod) {
                    if (informativeMessageCount < MAX_INFORMATIVE_MESSAGES) {
                        System.out.println(
                                "Skip spectrum " + spec.getID() +
                                        " since activationMethod is " + currentActivationMethod.toString() +
                                        ", not " + activationMethod.toString());
                        informativeMessageCount++;
                    } else {
                        if (informativeMessageCount == MAX_INFORMATIVE_MESSAGES) {
                            System.out.println(" ...");
                            informativeMessageCount++;
                        }
                    }
                    continue;
                }
            }

            if (!spec.isCentroided()) {
                if (informativeMessageCount < MAX_INFORMATIVE_MESSAGES) {
                    System.out.println("Skip spectrum " + spec.getID() + " since it is not centroided");
                    informativeMessageCount++;
                } else {
                    if (informativeMessageCount == MAX_INFORMATIVE_MESSAGES) {
                        System.out.println(" ...");
                        informativeMessageCount++;
                    }
                }
                numProfileSpectra++;
                continue;
            }

            if (spec.size() < minNumPeaksPerSpectrum) {
//				System.out.println("Spectrum " + spec.getScanNum() + " has too few peaks (#Peaks: " + spec.size()+"): ignored.");
                numSpectraWithTooFewPeaks++;
                continue;
            }

            if (charge == 0) {
                for (int c = minCharge; c <= maxCharge; c++)
                    specKeyList.add(new SpecKey(specIndex, c));
            } else if (charge > 0) {
                specKeyList.add(new SpecKey(specIndex, charge));
            }
        }

        System.out.println("Ignoring " + numProfileSpectra + " profile spectra.");
        System.out.println("Ignoring " + numSpectraWithTooFewPeaks + " spectra having less than " + minNumPeaksPerSpectrum + " peaks.");
        return specKeyList;
    }

    public static ArrayList<SpecKey> getFusedSpecKeyList(Iterator<Spectrum> itr, int startSpecIndex, int endSpecIndex, int minCharge, int maxCharge) {
        HashMap<Peak, ArrayList<Integer>> precursorSpecIndexMap = new HashMap<Peak, ArrayList<Integer>>();

        while (itr.hasNext()) {
            Spectrum spec = itr.next();
            int specIndex = spec.getSpecIndex();
            if (specIndex < startSpecIndex || specIndex >= endSpecIndex)
                continue;
            Peak precursor = spec.getPrecursorPeak();
            if (spec.getActivationMethod() == null) {
                System.out.println("Error: activation method is not available: Scan=" + spec.getSpecIndex() + ", PrecursorMz=" + spec.getPrecursorPeak().getMz());
                System.exit(-1);
            }

            ArrayList<Integer> list = precursorSpecIndexMap.get(precursor);
            if (list == null) {
                list = new ArrayList<Integer>();
                precursorSpecIndexMap.put(precursor, list);
            }
            list.add(specIndex);
        }

        Iterator<Entry<Peak, ArrayList<Integer>>> mapItr = precursorSpecIndexMap.entrySet().iterator();
        ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();
        while (mapItr.hasNext()) {
            Entry<Peak, ArrayList<Integer>> entry = mapItr.next();
            Peak precursor = entry.getKey();
            ArrayList<Integer> list = entry.getValue();
            Collections.sort(list);

            int charge = precursor.getCharge();
            if (charge == 0) {
                for (int c = minCharge; c <= maxCharge; c++) {
                    SpecKey specKey = new SpecKey(list.get(0), c);
                    for (int specIndex : list)
                        specKey.addSpecIndex(specIndex);
                    specKeyList.add(specKey);
                }
            } else if (charge > 0) {
                SpecKey specKey = new SpecKey(list.get(0), charge);
                for (int specIndex : list)
                    specKey.addSpecIndex(specIndex);
                specKeyList.add(specKey);
            } else {
                System.out.println("Error: negative precursor charge: " + precursor);
                System.exit(-1);
            }
        }
        return specKeyList;
    }

    public static void main(String argv[]) throws Exception {
        test();
    }

    public static void test() throws Exception {
        String fileName = "/home/sangtaekim/Research/Data/HeckRevision/CIDETDPairs/mzXML/090121_NM_Trypsin_20.mzXML";
        int minCharge = 2, maxCharge = 3;
        MzXMLSpectraIterator itr = new MzXMLSpectraIterator(fileName);
        ArrayList<SpecKey> list = SpecKey.getFusedSpecKeyList(itr, 0, Integer.MAX_VALUE, minCharge, maxCharge);
        for (SpecKey specKey : list) {
            if (specKey.getSpecIndexList() == null || specKey.getSpecIndexList().size() != 2)
                System.out.println(specKey.getSpecKeyString() + "\t" + specKey.getSpecIndexList());
        }
        System.out.println("Size: " + list.size());
    }

//	public static ArrayList<SpecKey> getFusedSpecKeyListOld(Iterator<Spectrum> itr, int minCharge, int maxCharge)
//	{
//		ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();
//		
//		int prevSpecIndex = Integer.MIN_VALUE;
//		int prevCharge = Integer.MIN_VALUE;
//		float previousPrecursorMz = Float.MIN_VALUE;
//		
//		while(itr.hasNext())
//		{
//			Spectrum spec = itr.next();
//			int specIndex = spec.getSpecIndex();
//			int charge = spec.getCharge();
//			float precursorMz = spec.getPrecursorPeak().getMz();
//			if(spec.getActivationMethod() == null)
//			{
//				System.out.println("Error: activation method is not available: Scan=" + spec.getSpecIndex()+", PrecursorMz=" + precursorMz);
//				System.exit(-1);
//			}
//			
//			if(specIndex == prevSpecIndex+1 && charge == prevCharge && precursorMz == previousPrecursorMz)
//			{
//				if(charge == 0)
//				{
//					for(int i=0; i<=maxCharge-minCharge; i++)
//						specKeyList.get(specKeyList.size()-1-i).addSpecIndex(specIndex);
//				}
//				else if(charge > 0)
//				{
//					specKeyList.get(specKeyList.size()-1).addSpecIndex(specIndex);
//				}
//			}
//			else
//			{
//				if(charge == 0)
//				{
//					for(int c=minCharge; c<=maxCharge; c++)
//					{
//						SpecKey specKey = new SpecKey(specIndex, c);
//						specKey.addSpecIndex(specIndex);
//						specKeyList.add(specKey);
//					}
//				}
//				else if(charge > 0)
//				{
//					SpecKey specKey = new SpecKey(specIndex, charge);
//					specKey.addSpecIndex(specIndex);
//					specKeyList.add(specKey);
//				}
//			}
//			prevSpecIndex = specIndex;
//			prevCharge = charge;
//			previousPrecursorMz = precursorMz;
//		}
//		return specKeyList;
//	}

}
