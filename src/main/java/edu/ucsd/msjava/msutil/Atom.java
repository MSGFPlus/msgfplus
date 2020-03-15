package edu.ucsd.msjava.msutil;

import java.util.HashMap;

public class Atom {
    public Atom(String code, double mass, int nominalMass, String name) {
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

    public static Atom get(String code) {
        return atomMap.get(code);
    }

    private final String code;
    private final String name;
    private final double mass;
    private final int nominalMass;

    private static final Atom[] atomArr =
            {
            /*
            Most of the following data can be automatically parsed out of the
            unimod.xml file from http://www.unimod.org/xml/unimod.xml by using the
            following regular expression and backreference replacement. It will
            not output correct nominal masses, those will need to be corrected by hand.
            
            (copy the entire contents of <umod:mod_bricks> to a separate text file; also, need to make each element use only one line)
            regex search: ^<umod:brick title=("[a-zA-Z0-9\-_]+") full_name=("[a-zA-Z0-9\-_ ]*") mono_mass="(\d+\.?\d*)" avge_mass="(\d+\.?\d*)"/?>$
            regex replace: new Atom(\1, \3, \3, \2),
            */
                    new Atom("-", 0, 0, ""), // Empty, should not be encountered, but we also don't want to error on it.
                    new Atom("H", 1.007825035, 1, "Hydrogen"),
                    new Atom("2H", 2.014101779, 2, "Deuterium"),
                    new Atom("Li", 7.016003, 7, "Lithium"),
                    new Atom("B", 11.0093055, 11, "Boron"),
                    new Atom("C", 12.0, 12, "Carbon"),
                    new Atom("13C", 13.00335483, 13, "Carbon 13"),
                    new Atom("N", 14.003074, 14, "Nitrogen"),
                    new Atom("15N", 15.00010897, 15, "Nitrogen 15"),
                    new Atom("O", 15.99491463, 16, "Oxygen"),
                    new Atom("18O", 17.9991603, 18, "Oxygen 18"),
                    new Atom("F", 18.99840322, 19, "Fluorine"),
                    new Atom("Na", 22.9897677, 23, "Sodium"),
                    new Atom("Mg", 23.9850423, 24, "Magnesium"),
                    new Atom("Al", 26.9815386, 27, "Aluminium"),
                    new Atom("P", 30.973762, 31, "Phosphorus"),
                    new Atom("S", 31.9720707, 32, "Sulfur"),
                    new Atom("Cl", 34.96885272, 35, "Chlorine"),
                    new Atom("K", 38.9637074, 39, "Potassium"),
                    new Atom("Ca", 39.9625906, 40, "Calcium"),
                    new Atom("Cr", 51.9405098, 52, "Chromium"),
                    new Atom("Mn", 54.9380471, 55, "Manganese"),
                    new Atom("Fe", 55.9349393, 56, "Iron"),
                    new Atom("Ni", 57.9353462, 58, "Nickel"),
                    new Atom("Co", 58.9331976, 59, "Cobalt"),
                    new Atom("Cu", 62.9295989, 63, "Copper"),
                    new Atom("Zn", 63.9291448, 64, "Zinc"),
                    new Atom("As", 74.9215942, 75, "Arsenic"),
                    new Atom("Br", 78.9183361, 79, "Bromine"),
                    new Atom("Se", 79.9165196, 80, "Selenium"),
                    new Atom("Mo", 97.9054073, 98, "Molybdenum"),
                    new Atom("Ru", 101.9043485, 102, "Ruthenium"),
                    new Atom("Pd", 105.903478, 106, "Palladium"),
                    new Atom("Ag", 106.905092, 107, "Silver"),
                    new Atom("Cd", 113.903357, 114, "Cadmium"),
                    new Atom("I", 126.904473, 127, "Iodine"),
                    new Atom("Pt", 194.964766, 195, "Platinum"),
                    new Atom("Au", 196.966543, 197, "Gold"),
                    new Atom("Hg", 201.970617, 202, "Mercury"),
                    // Unimod mod bricks, definitions from http://www.unimod.org/xml/unimod.xml
                    new Atom("Hex", 162.0528235, 162, "Hexose"),
                    new Atom("HexNAc", 203.079372605, 203, "N-Acetyl Hexosamine"),
                    new Atom("Ac", 42.0105647, 42, "Acetate"), // WARNING: SAME SYMBOL AS ACTINIUM!!!!
                    new Atom("dHex", 146.05790887, 146, "Deoxy-hexose"),
                    new Atom("HexA", 176.03208806, 176, "Hexuronic acid"),
                    new Atom("Kdn", 250.06886753, 250, "3-deoxy-d-glycero-D-galacto-nonulosonic acid"),
                    new Atom("Kdo", 220.05830283, 220, "2-keto-3-deoxyoctulosonic acid"),
                    new Atom("Me", 14.01565007, 14, "Methyl"),
                    new Atom("NeuAc", 291.095416635, 291, "N-acetyl neuraminic acid"),
                    new Atom("NeuGc", 307.09033126500003, 307, "N-glycoyl neuraminic acid"),
                    new Atom("Water", 18.0105647, 18, "Water"),
                    new Atom("Phos", 79.96633092500001, 80, "Phosphate"),
                    new Atom("Sulf", 79.95681459000001, 80, "Sulfate"),
                    new Atom("Pent", 132.0422588, 132, "Pentose"),
                    new Atom("Hep", 192.06338820000002, 192, "Heptose"),
                    new Atom("HexN", 161.068807905, 161, "Hexosamine"),
            };

    private static final HashMap<String, Atom> atomMap = new HashMap<String, Atom>();

    static {
        for (Atom atom : atomArr) {
            atomMap.put(atom.code, atom);
        }
    }
}
