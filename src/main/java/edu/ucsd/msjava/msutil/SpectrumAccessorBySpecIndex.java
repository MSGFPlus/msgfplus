package edu.ucsd.msjava.msutil;

import java.util.ArrayList;

public interface SpectrumAccessorBySpecIndex {
	public	Spectrum			getSpectrumBySpecIndex(int specIndex);
	public	String				getID(int specIndex);
	public	Float				getPrecursorMz(int specIndex); 
	public 	String				getTitle(int specIndex);
	public	ArrayList<Integer>	getSpecIndexList();
}
