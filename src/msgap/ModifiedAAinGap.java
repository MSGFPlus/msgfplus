package msgap;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Peptide;

public class ModifiedAAinGap {
	private class NDimBinaryMatrix{
		
		private BitSet matrix = null;
		private int totalDim, length, lengthPerDirection; // static?
		
		NDimBinaryMatrix(int totalDim, int maxLength, int maxLengthPerDirection){
			this.totalDim = totalDim;
			this.length = maxLength+1;
			this.lengthPerDirection = maxLengthPerDirection+1;
			matrix = new BitSet(pow(length, totalDim));
		}
		
		void setBit(int dim, int loc){ // dim < totalDim
			assert(dim < totalDim && loc < length);
			matrix.set(loc * pow(length, dim));
		}
		
		private void setBit(int[] coordinate, BitSet matrix){
			int bitSetIndex = 0;
			for(int dim = 0; dim < totalDim; dim++){
				int loc = coordinate[dim];
				//if(loc < length){
					bitSetIndex += loc * pow(length, dim);
				//}
			}
			matrix.set(bitSetIndex);
		}
		
		private int[] getCoordinate(int bitSetIndex){
			int[] cooradinate = new int[totalDim];
			
			for(int dim = totalDim-1; dim >=0; dim--){
				int divider = pow(length, dim);
				cooradinate[ dim ] = bitSetIndex / divider;
				bitSetIndex = bitSetIndex % divider; 
			 }
			
			return cooradinate;
		}
		
		void or(NDimBinaryMatrix other){
			matrix.or(other.matrix);
		}
		
		NDimBinaryMatrix getShiftedMatrix(int dim, int n){
			assert(dim < totalDim);
			//System.out.println(matrix);
			NDimBinaryMatrix ret = new NDimBinaryMatrix(totalDim, length-1, lengthPerDirection-1);
			
			for (int i = matrix.nextSetBit(0); i >= 0; i = matrix.nextSetBit(i+1)) {
				int[] cooradinate = getCoordinate(i);
				int newco = cooradinate[dim] + n;
				if(newco < lengthPerDirection){
					cooradinate[dim] += n;
					setBit(cooradinate, ret.matrix);
				}
			}
			return ret;
			//System.out.println(matrix);
		}
		
		
		ArrayList<int[]> getAllSetCoordinates(){
			ArrayList<int[]> ret = new ArrayList<int[]>();
			for (int i = matrix.nextSetBit(0); i >= 0; i = matrix.nextSetBit(i+1)) {
				int[] cor = getCoordinate(i);
				int sum=0;
				boolean write = true;
				for(int c:cor){
					if(c >= lengthPerDirection){
						write = false;
						break;
					}
					sum+=c;
				}
				if(write && sum < length) ret.add(cor);
			}
			
			return ret;
		}
		
		private int pow(int a, int b){
			int ret = 1;
			for(int i=0; i<b; i++)
				ret *= a;
			return ret;
		}
		
		public String toString(){
			String ret = "";
			for(int[] cor: getAllSetCoordinates()){
				ret+="(";
				for(int c : cor){
					ret +=" " + c + " ";
				}
				ret+= ") ";
			}
			return ret;
		}
	}
	
	/** TODO : maxModificationNumPerModification is not correctly working**/
	static private int maxModificationNum = 1; // inclusive
	static private int maxModificationNumPerModification = 1;
	static private AminoAcidSet aaSet = null;
	private static HashMap<AminoAcid, Integer> modifiedAAIndexMap = null;
	private static int[] massDiffs = null;
	private NDimBinaryMatrix modifiedAAMatrix = null;
	private static boolean isModifiedSearch = false;
	
//	private int[] modifiedAANum = null;
	
	public static AminoAcidSet getAASet() {return aaSet;}
	public static int getMaxModNum(){return maxModificationNum;}
	public static int getmaxModNumPerModification() {return maxModificationNumPerModification;}
	
	public ModifiedAAinGap(){
		//modifiedAANum = new int[modifiedAAIndexMap.size()];
		modifiedAAMatrix = new NDimBinaryMatrix(modifiedAAIndexMap.size(), maxModificationNum, maxModificationNumPerModification);
	}
	
	/*
	public HashMap<AminoAcid, ArrayList<Integer>> getModifiedAANumberMap(){
		HashMap<AminoAcid, ArrayList<Integer>> ret = new HashMap<AminoAcid, ArrayList<Integer>>();
		for(AminoAcid aa: modifiedAAIndexMap.keySet()){
			ArrayList<Integer> num = new ArrayList<Integer>();
			if(modifiedAAIndexMap.get(aa) != null){
				for(int i=0;i<maxModificationNum;i++){
					if(((modifiedAANum[modifiedAAIndexMap.get(aa)] >> i) & 0x01) == 0x01)
						num.add(i);
				}
			}
			ret.put(aa, num);
		}
		return ret;
	}
	*/
	
	/*
	// should return number of modifications
	private void subGetAllPossibleMassDifferences(int aaIndex, Integer mass, int numModification, HashMap<Integer, Integer> masses){
		if(mass != null && (!modifiedAAIndexMap.containsValue(aaIndex))){
			if(masses.containsKey(mass))
				numModification = numModification < masses.get(mass) ? numModification : masses.get(mass);
			masses.put(mass, numModification);
			return;
		}
			
		for(int i=0;i<=maxModificationNum;i++){
			if(((modifiedAANum[aaIndex] >> i) & 0x01) == 0x01){
			//	if(mass == null) mass = ;
			//	mass += i*massDiffs[aaIndex];
				if(numModification + i > maxModificationNum) continue;
				//subGetAllPossibleMassDifferences(aaIndex, i*massDiffs[aaIndex], numModification + i, masses);
				subGetAllPossibleMassDifferences(aaIndex+1, i*massDiffs[aaIndex], numModification + i, masses);
			}
				
		}
		
	}
	
	// return hashmap with massdif as key nummodification as value
	public  HashMap<Integer, Integer> getAllPossibleMassDifferences2(){
		 HashMap<Integer, Integer> masses = new HashMap<Integer, Integer>();
		subGetAllPossibleMassDifferences(0, null, 0, masses);
		
		return masses;
	}
	*/
	public  HashMap<Integer, Integer> getAllPossibleMassDifferences(){
		 HashMap<Integer, Integer> masses = new HashMap<Integer, Integer>();
		 
		
		 for(int[] coordinate : modifiedAAMatrix.getAllSetCoordinates()){

			 int massDiff = 0;
			 int numModification = 0;
			 for(int aaIndex = 0; aaIndex<coordinate.length; aaIndex++){
				 int aaNum = coordinate[aaIndex];
				 massDiff += aaNum*massDiffs[aaIndex];
				 numModification += aaNum;
			 }

			 if(masses.containsKey(massDiff))
				numModification = numModification < masses.get(massDiff) ? numModification : masses.get(massDiff);
			 
			 masses.put(massDiff, numModification);
				
		 }

		return masses;
	}
	
	static void initialize(AminoAcidSet as, int maxModNum, int maxModNumPerModification){
		int aaIndex = 0;
		modifiedAAIndexMap = new  HashMap<AminoAcid, Integer>();
		aaSet = as;
		for(AminoAcid aa : as){
			if(Character.isLowerCase(aa.getResidueStr())){ // change later
				modifiedAAIndexMap.put(aa, aaIndex++);
			}
		}

		massDiffs = new int[aaIndex];
		
		for (AminoAcid aa : modifiedAAIndexMap.keySet()){
			massDiffs[modifiedAAIndexMap.get(aa)] = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys().getAminoAcid(Character.toUpperCase(aa.getResidueStr())).getNominalMass() - aa.getNominalMass();
		}
		
		maxModificationNum = maxModNum;
		maxModificationNumPerModification = maxModNumPerModification;

		if(aaIndex > 0) isModifiedSearch =  true;
		else isModifiedSearch =  false;
		
	}
	
	static public boolean isAnyAAModified() {return isModifiedSearch;}
	
	public void setModifiedAANumber(AminoAcid aa){
		int aaIndex =  -1;
		
		if(modifiedAAIndexMap.get(aa) != null)
			aaIndex = modifiedAAIndexMap.get(aa);
		
		/*for(int i=0; i<modifiedAANum.length; i++){
			if(i == aaIndex) modifiedAANum[i] |= 0x01<<1;
			else modifiedAANum[i] |= 0x01;
		}
		*/

		if(aaIndex >= 0)
			modifiedAAMatrix.setBit(aaIndex, 1);
		else modifiedAAMatrix.setBit(0, 0); // all 0
		
	}
	
	public void addModifiedAADist(ModifiedAAinGap other, AminoAcid aa){
		int aaIndex =  -1;
		
		if(modifiedAAIndexMap.get(aa) != null)
			aaIndex = modifiedAAIndexMap.get(aa);
		
		/*
		for(int i=0; i<modifiedAANum.length; i++){
			if(i == aaIndex) modifiedAANum[i] |= other.modifiedAANum[i]<<1;
			else modifiedAANum[i] |= other.modifiedAANum[i];
		}
		
		*/
		NDimBinaryMatrix toOR = null;
		
		if(aaIndex >= 0){
			 toOR = other.modifiedAAMatrix.getShiftedMatrix(aaIndex, 1);
		}else
			 toOR  = other.modifiedAAMatrix;
		
		modifiedAAMatrix.or(toOR);
		
	}
	
	static private void subGetModifiedPeptides(ArrayList<AminoAcid> modifiedPeptide, int modNum, ArrayList<ArrayList<AminoAcid>> peptides){
		int pepLength = modifiedPeptide.size();
		
		if(pepLength >= peptides.get(0).size()){
			if(modNum>0) 
				peptides.add(modifiedPeptide);
			return;
		}
			
		int[] aaIndices = new int[modifiedAAIndexMap.size()];
		
		for(AminoAcid aa : modifiedPeptide){
			if(modifiedAAIndexMap.get(aa) == null) continue;
			int n = ++aaIndices[modifiedAAIndexMap.get(aa)];
			if(n > maxModificationNumPerModification)
				return;
		}
		
		AminoAcid nextAA = peptides.get(0).get(pepLength);
		
		if(modNum >= maxModificationNum){	
			ArrayList<AminoAcid> nextModifiedPeptide = new ArrayList<AminoAcid>(modifiedPeptide);
			nextModifiedPeptide.add(nextAA);
			subGetModifiedPeptides(nextModifiedPeptide, modNum, peptides);
		}else{
			for(AminoAcid aa : aaSet){
				if(Character.toUpperCase(aa.getResidueStr()) == nextAA.getResidueStr()){
					ArrayList<AminoAcid> nextModifiedPeptide = new ArrayList<AminoAcid>(modifiedPeptide);
					nextModifiedPeptide.add(aa);
					if(!aa.equals(nextAA))
						subGetModifiedPeptides(nextModifiedPeptide, modNum+1, peptides);
					else subGetModifiedPeptides(nextModifiedPeptide, modNum, peptides);
				}
			}
		}
	}

	// modified + unmodified
  static public ArrayList<Peptide> getModifiedPeptides(Peptide originalPeptide){
		ArrayList<Peptide> ret = new ArrayList<Peptide>();
		ArrayList<ArrayList<AminoAcid>> peptides = new ArrayList<ArrayList<AminoAcid>>();
		peptides.add(originalPeptide);
		
		subGetModifiedPeptides(new ArrayList<AminoAcid>(), 0, peptides);
		
		for(ArrayList<AminoAcid> pep : peptides){
			ret.add(new Peptide(pep));
		}

		return ret;
	}
  
	/*public String toString(){
		String out = new String();
		for(AminoAcid aa : modifiedAAIndexMap.keySet()){
			out+=aa+"\t"+ String.format("%x", modifiedAANum[modifiedAAIndexMap.get(aa)])+"\t";
		}
		return out;
	}*/
}
