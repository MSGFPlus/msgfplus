package msgf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Map.Entry;

import msutil.AminoAcid;
import msutil.Composition;
import msutil.Matter;
import msutil.Peptide;
import msutil.Sequence;

public class Profile<T extends Matter> extends ArrayList<ProfilePeak<T>>{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	// return profile peaks whose probability is equal or larger than threshold
	public Sequence<T> getNodesWithProbEqualOrHigherThan(float threshold)
	{
		Sequence<T> seq = new Sequence<T>();
		for(ProfilePeak<T> p : this)
		{
			if(p.getProbability() >= threshold)
				seq.add(p.getNode());
		}
		return seq;
	}
	
	public static Profile<Composition> getCompositionProfile(ArrayList<Peptide> dictionary, boolean prefix)
	{
		Hashtable<Composition, Integer> hist = new Hashtable<Composition, Integer>();

		for(Peptide peptide : dictionary)
		{ 
			Composition composition = new Composition(0);
			for(int i=0; i<peptide.size(); i++)
			{
				AminoAcid aa;
				if(prefix)
				  aa = peptide.get(i);
				else
				  aa = peptide.get(peptide.size()-1-i);
			composition = composition.getAddition(aa.getComposition());
			Integer occ = hist.get(composition);
			if(occ == null)
			hist.put(composition, 1);
			else
			hist.put(composition, occ+1);
			}
		}

		Profile<Composition> profile = new Profile<Composition>();
		for(Composition c : hist.keySet())
			profile.add(new ProfilePeak<Composition>(c, hist.get(c)/(float)dictionary.size()));

		Collections.sort(profile);
		
		return profile;
	}
	
	public Profile<NominalMass> toNominalMasses()
	{
		Profile<NominalMass> nominalMassProfile = new Profile<NominalMass>();
		Hashtable<Integer, Float> summedProfile = new Hashtable<Integer, Float>();
		for(ProfilePeak<T> p : this)
		{
			int mass = p.getNode().getNominalMass();
			float prob = p.getProbability();
			Float prevProb = summedProfile.get(mass);
			if(prevProb == null)
				summedProfile.put(mass, prob);
			else
				summedProfile.put(mass, prevProb+prob);
		}
		
		for(Integer mass : summedProfile.keySet())
		{
			float prob = summedProfile.get(mass);
			nominalMassProfile.add(new ProfilePeak<NominalMass>(NominalMassFactory.getInstanceFor(mass), prob));
		}
		
		Collections.sort(nominalMassProfile);
		return nominalMassProfile;
	}
	
	public String toString()
	{
		StringBuffer buf = new StringBuffer();
		for(ProfilePeak<T> p : this)
			buf.append(p.getNode().getMass() + "\t" + p.getProbability() + "\n");
		return buf.toString();
	}
	
	public Hashtable<T, Float> getHashtable()
	{
		Hashtable<T, Float> hashtable = new Hashtable<T, Float>();
		for(ProfilePeak<T> peak : this)
			hashtable.put(peak.getNode(), peak.getProbability());
		return hashtable;
	}
	
	public float getSumProbabilities()
	{
		float sumProb = 0;
		for(ProfilePeak<T> peak : this)
			sumProb += peak.getProbability();
		return sumProb;
	}
	
	public float getEuclideanDistance()
	{
		float dist = 0;
		for(ProfilePeak<T> peak : this)
			dist += peak.getProbability()*peak.getProbability();
		return (float)Math.sqrt(dist);
	}
	
	public Profile<T> getSubtraction(Profile<T> prof)
	{
		Profile<T> subtraction = new Profile<T>();
		Hashtable<T, Float> table = prof.getHashtable();
		for(ProfilePeak<T> peak : prof)
		{
			Float prob = table.get(peak.getNode()); 
			if(prob == null)	// only in prof
				table.put(peak.getNode(), peak.getProbability());
			else
				table.put(peak.getNode(), prob-peak.getProbability());
		}
		for(Entry<T, Float> entry : table.entrySet())
			subtraction.add(new ProfilePeak<T>(entry.getKey(), entry.getValue()));
		Collections.sort(subtraction);
		return subtraction;
	}
	
	public static<T extends Matter> float getDotProduct(Profile<T> prof1, Profile<T> prof2)
	{
		float dotProduct = 0; 
		Hashtable<T, Float> table1 = prof1.getHashtable();
		for(ProfilePeak<T> peak : prof2)
		{
			Float prob = table1.get(peak.getNode()); 
			if(prob != null)
				dotProduct += prob*peak.getProbability();
		}
		return dotProduct;
	}

	
	public static<T extends Matter> float getCosine(Profile<T> prof1, Profile<T> prof2)
	{
		return getDotProduct(prof1, prof2)/(prof1.getEuclideanDistance()*prof2.getEuclideanDistance());
	}
}
