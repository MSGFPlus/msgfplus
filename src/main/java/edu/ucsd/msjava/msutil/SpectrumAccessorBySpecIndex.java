package edu.ucsd.msjava.msutil;

import java.util.ArrayList;

public interface SpectrumAccessorBySpecIndex {
    Spectrum getSpectrumBySpecIndex(int specIndex);

    Spectrum getSpectrumById(String specId);

    String getID(int specIndex);

    Float getPrecursorMz(int specIndex);

    String getTitle(int specIndex);

    ArrayList<Integer> getSpecIndexList();
}
