package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.misc.CountPSMs;

public class TestMSUtils {
	@Test
	public void countPSMs()
	{
		File tsvFile = new File("/Users/sangtaekim/Research/Data/Viktor/103111-Yeast-2hr-01.ANNOTATED.tsv");
		try {
			CountPSMs.countID(tsvFile.getPath(), 0.02f);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
