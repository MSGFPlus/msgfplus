package edu.ucsd.msjava.msutil;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;

public class Composition extends Matter
{
	public static final double C = 12.0f;
	public static final double C13 = 13.00335483;
	public static final double C14 = 14.003241;
	public static final double H = 1.007825035;
	public static final double DEUTERIUM = 2.014101779;	
	public static final double N = 14.003074;
	public static final double N15 = 15.000108898;
	public static final double O = 15.99491463;
	public static final double S = 31.9720707;
	public static final double P = 30.973762;
	public static final double Br = 78.9183361;
	public static final double Cl = 34.96885272;
	public static final double Fe = 55.9349393;
	public static final double Se = 79.9165196;
	
	public static final double H2 = H*2;
	public static final double NH = N+H;
	public static final double NH2 = N+2*H;
	public static final double H2O = H*2+O;
	public static final double NH3 = N+H*3;
	public static final double CO = C+O;
	public static final double ISOTOPE = C13 - C;
	public static final double ISOTOPE2 = C14 - C;
//	public static final double PROTON = 1.0072697;
	public static final double PROTON = 1.007825035;
	public static final double NEUTRON = 1.0086650;
	public static final double OFFSET_Y = H*3+O;
	public static final double OFFSET_B = PROTON;

	public static final Composition NIL = new Composition(0,0,0,0,0);

	int number; // MSB 8 (C) 8 (H) 6 (N) 6 (O) 4 (S)
	static final double[] monoMass = new double[] { C, H, N, O, S};
	static final float[] avgMass = new float[] { 12.011f, 1.00794f, 14.00674f, 15.9994f, 32.066f };


	//  private Composition() {}

	public Composition(int C, int H, int N, int O, int S)
	{
		number = C*0x01000000 + H*0x00010000 + N*0x00000400 + O*0x00000010 + S; 
	}

	public Composition(int number)
	{
		this.number = number;
	}

	public Composition(Composition c)
	{
		this.number = c.number;
	}
	
	public Composition(String compositionStr)
	{
		HashMap<Character, Integer> compTable = new HashMap<Character, Integer>();
		compTable.put('C', 0);
		compTable.put('H', 0);
		compTable.put('N', 0);
		compTable.put('O', 0);
		compTable.put('S', 0);

		int number = 0;
		char element = '*';
		int i=0;
		while(i<compositionStr.length())
		{
			char c = compositionStr.charAt(i);
			if(Character.isLetter(c))
			{
				if(number > 0)
					compTable.put(element, number);
				element = c;
				number = 0;
			}
			else if(Character.isDigit(c))
			{
				number = 10*number + Integer.parseInt(String.valueOf(c));
			}
			i++;
		}
		if(number > 0)
			compTable.put(element, number);
		this.number = new Composition(compTable.get('C'),
				compTable.get('H'),compTable.get('N'),compTable.get('O'),compTable.get('S')).number;
		
	}
	//  public static final Composition getInstance(int C, int H, int N, int O, int S)
	// {
	//	  int number = C*0x01000000 + H*0x00010000 + N*0x00000400 + O*0x00000010 + S; 
	// }

	// public static final Composition getInstance(int number)
	// {
	//  }

	public int getC() { return (number & 0xFF000000) >>> 24; }
	public int getH() { return (number & 0x00FF0000) >> 16; }
	public int getN() { return (number & 0x0000FC00) >> 10; } 
	public int getO() { return (number & 0x000003F0) >> 4; }
	public int getS() { return (number & 0x0000000F); }
	public int getNumber()  { return number; }
	//  public int getIndex()	{ return number; }
	
	@Override
	public int hashCode()
	{
		return number;
	}

	public static float getMonoMass(int number)
	{
		return (float)(
				((number & 0xFF000000) >>> 24) * Composition.C +
				((number & 0x00FF0000) >> 16) * Composition.H +
				((number & 0x0000FC00) >> 10) * Composition.N +   
				((number & 0x000003F0) >> 4) * Composition.O + 
				(number & 0x0000000F) * Composition.S);
	}

	public static float getAvgMass(int number)
	{
		return 
		((number & 0xFF000000) >>> 24)*avgMass[0] + 
		((number & 0x00FF0000) >> 16)*avgMass[1] + 
		((number & 0x0000FC00) >> 10)*avgMass[2] + 
		((number & 0x000003F0) >> 4)*avgMass[3] + 
		(number & 0x0000000F)*avgMass[4];
	}

	@Override
	public float getMass()
	{
		return (float)(getC()*Composition.C+getH()*Composition.H+getN()*Composition.N+getO()*Composition.O+getS()*Composition.S);
	}

	@Override
	public double getAccurateMass()
	{
		return (getC()*Composition.C+getH()*Composition.H+getN()*Composition.N+getO()*Composition.O+getS()*Composition.S);
	}

	public int getNominalMass()
	{
		return getC()*12+getH()*1+getN()*14+getO()*16+getS()*32;
	}

	public float getAvgMass()
	{
		return getC()*avgMass[0]+getH()*avgMass[1]+getN()*avgMass[2]+getO()*avgMass[3]+getS()*avgMass[4];
	}

	public String toString()
	{
		return new String(getC()+" "+getH()+" "+getN()+" "+getO()+" "+getS());
	}

	public void add(Composition c)
	{
		number += c.number;
	}
	
	public Composition getAddition(Composition c)
	{
		return new Composition(number+c.number);
	}

	public Composition getSubtraction(Composition c)
	{
		int newC = getC()-c.getC();
		int newH = getH()-c.getH();
		int newN = getN()-c.getN();
		int newO = getO()-c.getO();
		int newS = getS()-c.getS();
		
		if(newC < 0 || newH < 0 || newN < 0 || newO < 0 || newS < 0)
			return null;
		return new Composition(newC, newH, newN, newO, newS);
	}
	
	public boolean equals(Object o)
	{
		if(o instanceof Composition)
		{
			Composition c = (Composition)o;
			if(number == c.number)
				return true;
		}
		return false;
	}

	public static Double getMass(String compositionStr)
	{
		if(!compositionStr.matches("(([A-Z][a-z]?([+-]\\d+|\\d*)))+"))
			return null;
		
		HashMap<String, Integer> compTable = new HashMap<String, Integer>();
		compTable.put("C", 0);
		compTable.put("H", 0);
		compTable.put("N", 0);
		compTable.put("O", 0);
		compTable.put("S", 0);
		compTable.put("P", 0);
		compTable.put("Br", 0);
		compTable.put("Cl", 0);
		compTable.put("Fe", 0);
		compTable.put("Se", 0);

		int i=0;
		while(i<compositionStr.length())
		{
			int j = i;
			String atom;
			if(i+1<compositionStr.length() && Character.isLowerCase(compositionStr.charAt(i+1)))
				j += 2;
			else
				j += 1;
			
			atom = compositionStr.substring(i, j);
			
			i = j;
			
			Integer number = compTable.get(atom);
			if(number == null || !number.equals(0))
				return null;
			
			while(j < compositionStr.length())
			{
				char c = compositionStr.charAt(j);
				if(c != '+' && c != '-' && !Character.isDigit(c))
					break;
				else
					j++;
			}
			
			int n;
			if(j == i)
				n = 1;
			else
				n = Integer.parseInt(compositionStr.substring(i, j));
					
			compTable.put(atom, n);
			i = j;
		}
		
		double modMass = 
			compTable.get("C")*Composition.C
			+ compTable.get("H")*Composition.H
			+ compTable.get("N")*Composition.N
			+ compTable.get("O")*Composition.O
			+ compTable.get("S")*Composition.S
			+ compTable.get("P")*Composition.P
			+ compTable.get("Br")*Composition.Br
			+ compTable.get("Cl")*Composition.Cl
			+ compTable.get("Fe")*Composition.Fe
			+ compTable.get("Se")*Composition.Se;
		return modMass;
	}
	
	/*
  public int compareTo(Composition c) 
  {
    float diff = getMass() - c.getMass();
    if(diff == 0)
      return this.number - c.number;
    else if(diff > 0)
      return 1;
    else
      return -1;
  }
	 */

	public static class CompositionComparator implements Comparator<Integer> {
		public int compare(Integer c1, Integer c2)
		{
			double mass1 = Composition.getMonoMass(c1);
			double mass2 = Composition.getMonoMass(c2);
			if(mass1 > mass2)
				return 1;
			else if(mass1 < mass2)
				return -1;
			else
			{
				return c1-c2;
			}
		}
		public boolean equals(Integer c1, Integer c2)
		{
			return (c1 == c2);
		}
	}


	/**
	 * Comparator method for 2 edges in composition representation. The order is
	 * defined by the mass of the edges and then by the composition itself.
	 * @param comp1 the composition of the first edge.
	 * @param comp2 the composition of the second edge.
	 * @return positive if the second edge is greater than the first one, negative
	 *                  if the reverse is true and 0 if they are equal.
	 */
	public static int compareCompositions(int comp1, int comp2) {
		double mass1 = Composition.getMonoMass(comp1), mass2 = Composition.getMonoMass(comp2);    
		if (mass1 < mass2)       return -1;
		if (mass2 < mass1)       return 1;
		if (comp1 < comp2)       return -1;
		if (comp2 < comp1)       return 1;
		return 0;
	}


	/**
	 * Check compositions sums that are equal to another standard composition.
	 */
	private static void checkEquality() {
		// check which compositions are the sum of any two standard composition
		AminoAcid[] stdAa = AminoAcid.getStandardAminoAcids(); 
		Composition[] stdComp = new Composition[stdAa.length];
		for (int i=0; i<stdAa.length; i++) {
			stdComp[i] = stdAa[i].getComposition();
		}

		System.out.println("Composition equalities: ");
		for (int i=0; i<stdAa.length; i++) {
			for (int j=i; j<stdAa.length; j++) {
				Composition sum = stdComp[i].getAddition(stdComp[j]);
				for (int k=0; k<stdAa.length; k++) {
					if (sum.equals(stdComp[k])) {
						System.out.println(stdAa[i].toString() + " plus " + stdAa[j].toString() + " equals " + stdAa[k].toString());
					}
				}
			}
		}

		int[] singleMasses = new int[stdAa.length];
		for (int i=0; i<stdAa.length; i++) {
			singleMasses[i] = stdAa[i].getNominalMass();
		}
		System.out.println("Integer equalities: ");
		for (int i=0; i<stdAa.length; i++) {
			for (int j=i; j<stdAa.length; j++) {
				int sum = stdComp[i].getNominalMass() + stdComp[j].getNominalMass();
				for (int k=0; k<stdAa.length; k++) {
					if (sum==stdComp[k].getNominalMass()) {
						System.out.println(stdAa[i].toString() + " plus " + stdAa[j].toString() + " equals " + stdAa[k].toString());
					}
				}
			}
		}

	}

	public static void main(String argv[])
	{
		/*
    Composition[] aa = {
      new Composition(2,3,1,1,0),
      new Composition(3,5,1,1,0),
      new Composition(3,5,1,2,0),
      new Composition(5,7,1,1,0),
      new Composition(5,9,1,1,0),
      new Composition(4,7,1,2,0),
      new Composition(3,5,1,1,1),
      new Composition(6,11,1,1,0),
      new Composition(4,6,2,2,0),
      new Composition(4,5,1,3,0),
      new Composition(5,8,2,2,0),
      new Composition(6,12,2,1,0),
      new Composition(5,7,1,3,0),
      new Composition(5,9,1,1,1),Serializable, 
      new Composition(6,7,3,1,0),
      new Composition(9,9,1,1,0),
      new Composition(6,12,4,1,0),
      new Composition(9,9,1,2,0),
      new Composition(11,10,2,1,0),
    };
		 */
		checkEquality();
	}
}