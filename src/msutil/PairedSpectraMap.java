package msutil;

import java.util.ArrayList;

public class PairedSpectraMap implements SpectrumAccessorByScanNum {
	private SpectrumAccessorByScanNum specMap;
	private ArrayList<Integer> scanNumList;
	
	public PairedSpectraMap(SpectraIterator itr)
	{
		this.specMap = specMap;
	}
	
	private void makeScanNumList()
	{
		int prevScanNum = Integer.MAX_VALUE;
		for(int scanNum : specMap.getScanNumList())
		{
			
			Spectrum spec = specMap.getSpectrumByScanNum(scanNum);
		}
	}

	@Override
	public ArrayList<Integer> getScanNumList() {
		return scanNumList;
	}

	@Override
	public Spectrum getSpectrumByScanNum(int scanNum) {
		return specMap.getSpectrumByScanNum(scanNum);
	}
}
