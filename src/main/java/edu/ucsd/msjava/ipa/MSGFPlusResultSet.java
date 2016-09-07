package edu.ucsd.msjava.ipa;

import edu.ucsd.msjava.parser.BufferedLineReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MSGFPlusResultSet {
    private String header;
    private List<PSM> psmList;

    public MSGFPlusResultSet(File msgfPlusResultFile) {
        parse(msgfPlusResultFile);
    }

    public String getHeader() {
        return header;
    }

    public List<PSM> getPSMList() {
        return psmList;
    }

    private void parse(File msgfPlusResultFile) {
        psmList = new ArrayList<PSM>();
        BufferedLineReader in = null;
        try {
            in = new BufferedLineReader(msgfPlusResultFile.getPath());
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        String s;

        header = in.readLine();    // header
        while ((s = in.readLine()) != null) {
            PSM psm = new PSM(s);
            psmList.add(psm);
        }

        try {
            in.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
