package msutil;

import java.util.ArrayList;

public interface SpectrumAccessorByScanNum {
	public Spectrum getSpectrumByScanNum(int scanNum);
	public ArrayList<Integer> getScanNumList();
}
