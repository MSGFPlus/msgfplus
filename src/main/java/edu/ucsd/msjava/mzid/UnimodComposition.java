package edu.ucsd.msjava.mzid;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;

import edu.ucsd.msjava.msutil.Atom;
import edu.ucsd.msjava.msutil.Composition;

public class UnimodComposition {
	
	
	public UnimodComposition()
	{
		compMap = new LinkedHashMap<String,Integer>();
		compMap.put("H", 0);
		compMap.put("C", 0);
		compMap.put("N", 0);
		compMap.put("O", 0);
		compMap.put("P", 0);
		compMap.put("S", 0);
	}
	
	public void add(Composition comp)
	{
		add("C", comp.getC());
		add("H", comp.getH());
		add("N", comp.getN());
		add("O", comp.getO());
		add("S", comp.getS());
	}
	
	public void add(String deltaComposition)
	{
		String[] token = deltaComposition.split("\\s+");
		for(String e : token)
		{
			if(e.matches("\\d*?[a-zA-Z]+(\\(-?\\d+\\))?"))
			{
				String element;
				int num;
				if(e.matches("\\d*?[a-zA-Z]+"))
				{
					element = e;
					num = 1;
				}
				else
				{
					element = e.substring(0, e.indexOf('('));
					num = Integer.parseInt(e.substring(e.indexOf('(')+1, e.lastIndexOf(')')));
				}
				add(element, num);
			}
			else if(e.matches("\\d+\\.?\\d*"))
			{
				double mass = Double.parseDouble(e);
				add(mass);
			}
			else
			{
				System.err.println("Wrong Unimod delta_composition: " + deltaComposition);
				System.exit(-1);
			}
		}
	}
	
	public void add(String element, int number)
	{
		Integer num = compMap.get(element);
		if(num == null)
			compMap.put(element, number);
		else
			compMap.put(element, num+number);
	}
	
	public void add(double deltaMass)
	{
		if(this.deltaMass == null)
			this.deltaMass = deltaMass;
		else
			this.deltaMass += deltaMass;
	}
	
	public Double getMass()
	{
		double mass = 0;
		Iterator<Entry<String, Integer>> itr = compMap.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Integer> entry = itr.next();
			String element = entry.getKey();
			int num = entry.getValue();
			if(num == 0)
				continue;
			Atom atom = Atom.get(element);
			if(atom == null)
				return null;
			mass += atom.getMass()*num;
		}
		
		if(deltaMass != null)
			mass += deltaMass;
		return mass;
	}
	
	public static Double getMass(String unimodCompositionStr)
	{
		UnimodComposition comp = new UnimodComposition();
		comp.add(unimodCompositionStr);
		return comp.getMass();
	}
	
	@Override
	public String toString()
	{
		StringBuffer buf = new StringBuffer();
		Iterator<Entry<String, Integer>> itr = compMap.entrySet().iterator();
		boolean first = true;
		while(itr.hasNext())
		{
			Entry<String, Integer> entry = itr.next();
			String element = entry.getKey();
			int num = entry.getValue();
			if(num == 0)
				continue;
			else if(num == 1)
			{
				if(!first)
					buf.append(" ");
				else
					first = false;
				buf.append(element);
			}
			else
			{
				if(!first)
					buf.append(" ");
				else
					first = false;
				buf.append(element+"("+num+")");
			}
		}
		
		if(deltaMass != null)
			buf.append(" " + deltaMass);
		return buf.toString();
	}
	
	private Map<String,Integer> compMap;
	private Double deltaMass = null;
	
}
