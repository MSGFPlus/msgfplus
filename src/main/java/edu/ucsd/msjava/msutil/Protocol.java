package edu.ucsd.msjava.msutil;

import edu.ucsd.msjava.params.ParamObject;
import edu.ucsd.msjava.params.UserParam;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;


public class Protocol implements ParamObject {
    private String name;
    private String description;

    private Protocol(String name, String description) {
        this.name = name;
        this.description = description;
    }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public String getParamDescription() {
        return name;
    }

    // static members
    public static Protocol get(String name) {
        return table.get(name);
    }

    public static final Protocol AUTOMATIC;
    public static final Protocol PHOSPHORYLATION;
    public static final Protocol ITRAQ;
    public static final Protocol ITRAQPHOSPHO;
    public static final Protocol TMT;
    public static final Protocol STANDARD;

    public static Protocol[] getAllRegisteredProtocols() {
        return protocolList.toArray(new Protocol[0]);
    }

    private static HashMap<String, Protocol> table;
    private static ArrayList<Protocol> protocolList;

    private static void add(Protocol prot) {
        if (table.put(prot.name, prot) == null)
            protocolList.add(prot);
    }

    static {
        AUTOMATIC = new Protocol("Automatic", "Automatic");
        PHOSPHORYLATION = new Protocol("Phosphorylation", "Phospho-enriched");
        ITRAQ = new Protocol("iTRAQ", "iTRAQ");
        ITRAQPHOSPHO = new Protocol("iTRAQPhospho", "iTRAQPhospho");
        TMT = new Protocol("TMT", "TMT");
        STANDARD = new Protocol("Standard", "Standard");

        table = new HashMap<String, Protocol>();
        protocolList = new ArrayList<Protocol>();

        protocolList.add(AUTOMATIC);
        add(PHOSPHORYLATION);
        add(ITRAQ);
        add(ITRAQPHOSPHO);
        add(TMT);
        add(STANDARD);

        // Parse activation methods defined by a user
        File protocolFile = Paths.get("params", "protocols.txt").toFile();
        if (protocolFile.exists()) {
            ArrayList<String> paramLines = UserParam.parseFromFile(protocolFile.getPath(), 2);
            for (String paramLine : paramLines) {
                String[] token = paramLine.split(",");
                String shortName = token[0];
                String description = token[1];
                Protocol newProt = new Protocol(shortName, description);
                add(newProt);
            }
        }
    }

}
