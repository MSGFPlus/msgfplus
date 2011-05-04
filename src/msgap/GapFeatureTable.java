package msgap;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import msgf.ScoreDist;
import msutil.Matter;

public class GapFeatureTable<T  extends Matter> {
	private HashMap<T, ScoreDist[]> scoreDistMap = null;
	private HashMap<T, ModifiedAAinGap[]> modifiedAAMap = null;
	
	private int numHub;
	
	public GapFeatureTable(int numHub){
		scoreDistMap = new HashMap<T, ScoreDist[]>(1000);
		modifiedAAMap = new HashMap<T, ModifiedAAinGap[]>();
		this.numHub = numHub;
	}
	
	public boolean containsNode(T node){
		return scoreDistMap.containsKey(node);
	}
		
	public ScoreDist getDistBetween(int hubIndex, T node){
		if(scoreDistMap.get(node) == null) return null;
		return scoreDistMap.get(node)[hubIndex];
	}
	
	public void putDist(int hubIndex, T node, ScoreDist dist){
		ScoreDist[] d = scoreDistMap.get(node);
		boolean isnew = false;
		if(d == null){
			d = new ScoreDist[numHub];
			isnew = true;
		}
		
		d[hubIndex] = dist;
		if(isnew) scoreDistMap.put(node, d);
	}
	
	public void putModifiedAADist(int hubIndex, T node,  ModifiedAAinGap aaDist){ 
		ModifiedAAinGap[] d = modifiedAAMap.get(node);
		boolean isnew = false;
		if(d == null){
			d = new ModifiedAAinGap[numHub];
			isnew = true;
		}
		
		d[hubIndex] = aaDist;
		if(isnew) modifiedAAMap.put(node, d);
	}
	
	public ModifiedAAinGap getModifiedAADistBetween(int hubIndex, T node){
		if(modifiedAAMap.get(node) == null) return null;
		return modifiedAAMap.get(node)[hubIndex];
	}
	
	
	public void remove(T node){
		scoreDistMap.remove(node);
		modifiedAAMap.remove(node);
	}
	
	
	public ArrayList<Integer> getConnectedHubIndicesFrom(T node){
		ScoreDist[] d = scoreDistMap.get(node);
		if(d == null) return null;
		
		ArrayList<Integer>  indices = new ArrayList<Integer>();
		
		for(int i=0;i<d.length;i++)
			if(d[i] != null) indices.add(i);
				
		return indices;
	}
	
}