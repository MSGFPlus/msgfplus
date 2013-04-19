package edu.ucsd.msjava.msutil;

import java.util.HashMap;

public class Atom {
	public Atom(String code, double mass, int nominalMass, String name)
    {
    	this.code = code;
        this.mass = mass;
        this.nominalMass = nominalMass;
        this.name = name;
    }

    public String getCode() {
		return code;
	}

	public String getName() {
		return name;
	}

	public double getMass() {
		return mass;
	}

	public int getNominalMass() {
		return nominalMass;
	}

	public static Atom[] getAtomarr() {
		return atomArr;
	}

	public static HashMap<String, Atom> getAtomMap() {
		return atomMap;
	}
	
    public static Atom get(String code)
    {
        return atomMap.get(code);
    }
	
    private final String code;
    private final String name;
    private final double mass;
    private final int nominalMass;
    
    private static final Atom[] atomArr =
        {
            new Atom("H", 1.007825035, 1, "Hydrogen"),
            new Atom("2H", 2.014101779, 2, "Deuterium"),
            new Atom("Li", 7.016003, 7, "Lithium"),
            new Atom("C", 12.0, 12, "Carbon"),
            new Atom("13C", 13.00335483, 13, "Carbon13"),
            new Atom("N", 14.003074, 14, "Nitrogen"),
            new Atom("15N", 15.00010897, 15, "Nitrogen15"),
            new Atom("O", 15.99491463, 16, "Oxigen"),
            new Atom("18O", 17.9991603, 18, "Oxigen"),
            new Atom("F", 18.99840322, 19, "Fluorine"),
            new Atom("Na", 22.9897677, 23, "Sodium"),
            new Atom("P", 30.973762, 13, "Phosphorous"),
            new Atom("S", 31.9720707, 32, "Sulfur"),
            new Atom("Cl", 34.96885272, 35, "Chlorine"),
            new Atom("K", 38.9637074, 39, "Potassium"),
            new Atom("Ca", 39.9625906, 40, "Calcium"),
            new Atom("Fe", 55.9349393, 56, "Iron"),
            new Atom("Ni", 57.9353462, 58, "Nickel"),
            new Atom("Cu", 62.9295989, 63, "Copper"),
            new Atom("Zn", 63.9291448, 64, "Zinc"),
            new Atom("Br", 78.9183361, 79, "Bromine"),
            new Atom("Se", 79.9165196, 80, "Selenium"),
            new Atom("Mo", 97.9054073, 98, "Molybdenum"),
            new Atom("Ag", 106.905092, 107, "Silver"),
            new Atom("I", 126.904473, 127, "Iodine"),
            new Atom("Au", 196.966543, 197, "Gold"),
            new Atom("Hg", 201.970617, 202, "Mercury"),
            new Atom("Hex", 162.052824, 162, "Hexose"),
            new Atom("HexNAc", 203.079373, 203, "N-Acetylhexosamine"),
            new Atom("dHex", 146.057909, 146, "Fucose"),
            new Atom("NeuAc", 291.095417, 291, "N-acetyl neuraminic acid"),
            new Atom("NeuGc", 307.090331, 307, "N-glycoyl neuraminic acid"),
            new Atom("Hep", 192.063388, 192, "Heptose"),
            new Atom("Pent", 132.042257, 85, "Pentose"),
            new Atom("B", 11.00930554, 11, "Boron"),
            new Atom("As", 74.92159, 75, "Arsenic"),
            new Atom("Mg",  23.985043, 24, "Magnesium")
        };

    private static final HashMap<String,Atom> atomMap = new HashMap<String, Atom>();
    
    static 
    {
        for(Atom atom : atomArr)
        {
        	atomMap.put(atom.code, atom);
        }
    }
}
