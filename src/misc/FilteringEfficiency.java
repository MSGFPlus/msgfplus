package misc;

import msutil.AminoAcid;
import msutil.AminoAcidSet;

public class FilteringEfficiency {
	public static void main(String atgv[])
	{
		float[] gf = new float[1000];
		
		int[] aaMass = new int[20];
		int aaIndex = -1;
		for(AminoAcid aa : AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys())
		{
			aaIndex++;
			aaMass[aaIndex] = aa.getNominalMass();
		}
		
		gf[0] = 1;
		for(int i=1; i<gf.length; i++)
		{
			for(int j=0; j<aaMass.length; j++)
			{
				int prevMass = i-aaMass[j];
				if(prevMass < 0)
					continue;
				gf[i] += 0.05*gf[prevMass];
			}
		}
		for(int i=0; i<gf.length; i++)
			System.out.println(i+"\t"+gf[i]);
		
		int sum = 0;
		for(int m : aaMass)
			sum += m;
		System.out.println("AverageAAMass: " + sum/(float)20);
	}
}
