package msutil;

import java.util.ArrayList;

public interface SpectrumAccessorBySpecIndex {
	public Spectrum getSpectrumBySpecIndex(int specIndex);
	public ArrayList<Integer> getSpecIndexList();
}
