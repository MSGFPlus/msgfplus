package msutil;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

import msgf.DeNovoGraph;
import msgf.IntMassFactory;
import msgf.MassFactory;
import msgf.Tolerance;

/**
 * A factory class instantiate compositions.
 * @author sangtaekim
 *
 */
public class CompositionFactory extends MassFactory<Composition>{
	
	private static final int arraySize = 1 << 27;
	private static final int indexMask = 0xFFFFFFE0;
	private static final int offsetMask = 0x0000001F;
	
	private int[] map;
	private ArrayList<Composition> tempData;	// temporary
	private int[] data;
	
	public CompositionFactory(AminoAcidSet aaSet, Enzyme enzyme, int maxLength)
	{
		super(aaSet, enzyme, maxLength);
		this.map = new int[arraySize];
		tempData = new ArrayList<Composition>();
		makeAllPossibleMasses();
	}
	
	// private class for getIntermediateCompositions, don't generate all possible nodes
	private CompositionFactory(AminoAcidSet aaSet, int maxLength)
	{
		super(aaSet, null, maxLength);
		this.map = new int[arraySize];
		tempData = new ArrayList<Composition>();
	}
	
	@Override
	public Composition getZero() {
		return Composition.NIL;
	}

	public Composition getNextNode(Composition curNode, AminoAcid aa) {
		int num = curNode.number + aa.getComposition().number;
		return new Composition(num);
	}

	public Composition getComplementNode(Composition srm, Composition pmNode) {
		return pmNode.getSubtraction(srm);
	}
	
	public ArrayList<DeNovoGraph.Edge<Composition>> getEdges(Composition curNode) {
		// prevNode, score, prob, index
		int curNum = curNode.number;
		ArrayList<DeNovoGraph.Edge<Composition>> edges = new ArrayList<DeNovoGraph.Edge<Composition>>();
		for(AminoAcid aa : aaSet)
		{
			int prevNum = curNum - aa.getComposition().number;
			DeNovoGraph.Edge<Composition> edge = new DeNovoGraph.Edge<Composition>(new Composition(prevNum), aa.getProbability(), aaSet.getIndex(aa), aa.getMass()); 
			if(prevNum == 0 && enzyme != null)
			{
				if(enzyme.isCleavable(aa))
					edge.setCleavageScore(aaSet.getPeptideCleavageCredit());
				else
					edge.setCleavageScore(aaSet.getPeptideCleavagePenalty());
			}
			edges.add(edge);
		}
		return edges;
	}
	
	@Override
	public int size()	
	{ 
		if(data == null)	// not finalized yet
		{
			if(tempData == null)
				return -1;
			else
				return tempData.size();
		}
		else 
			return data.length; 
	}
	
	public int[] getData()	{ return data; }
	
	public ArrayList<Composition> getNodes(float mass, Tolerance tolerance)
	{
		ArrayList<Composition> compositions = new ArrayList<Composition>();
	
		float toleranceDa = tolerance.getToleranceAsDa(mass);
		float minMass = mass-toleranceDa;
		float maxMass = mass+toleranceDa;
		// binary search
		int minIndex=0, maxIndex=data.length, i=-1;
		while(true)
		{
			i = (minIndex+maxIndex)/2; 
			double m = Composition.getMonoMass(data[i]);
			if(m < minMass)
				minIndex = i;
			else if(m > maxMass)
				maxIndex = i;
			else
				break;
			if(maxIndex - minIndex <= 1)
				break;
		}
		for(int cur=i; cur>=0; cur--)
		{
			double m = Composition.getMonoMass(data[cur]);
			if(m >= minMass && m <= maxMass)
				compositions.add(new Composition(data[cur]));
			else if(m < minMass)
				break;
		}
		for(int cur=i+1; cur<data.length; cur++)
		{
			double m = Composition.getMonoMass(data[cur]);
			if(m >= minMass && m <= maxMass)
				compositions.add(new Composition(data[cur]));
			else if(m > maxMass)
				break;
		}
		Collections.sort(compositions);
		return compositions;
	}

	public Composition getNode(float mass) {
		// binary search
		int minIndex=0, maxIndex=data.length, i=-1;
		while(true)
		{
			i = (minIndex+maxIndex)/2; 
			double m = Composition.getMonoMass(data[i]);
			if(m < mass)
				minIndex = i;
			else if(m > mass)
				maxIndex = i;
			else
				break;
			if(maxIndex - minIndex <= 1)
				break;
		}

		if(minIndex == maxIndex)
			return new Composition(data[minIndex]);
		else
		{
			Composition compMin = new Composition(data[minIndex]);
			Composition compMax = new Composition(data[maxIndex]);
			float min = compMin.getMass();
			float max = compMax.getMass();
			if(Math.abs(mass-min) < Math.abs(mass-max))
				return compMin;
			else
				return compMax;
		}
	}
	
	@Override
	public ArrayList<Composition> getLinkedNodeList(Collection<Composition> destCompositionList)
	{
		return getIntermediateCompositions(new Composition(0), destCompositionList);
	}
	
	// return set of compositions contained in paths from (0,0,0,0,0) to despCompositions
	public ArrayList<Composition> getIntermediateCompositions(Composition source, Collection<Composition> destCompositionList)
	{
		CompositionFactory intermediateCompositions = new CompositionFactory(this.aaSet, maxLength);
		
		for(Composition c : destCompositionList)
		{
			intermediateCompositions.setAndAddIfNotExist(c.number);
		}
		
		int start = 0;
		while(true)
		{
			int end = intermediateCompositions.size();
			for(int i=start; i<end; i++)
			{
				int number = intermediateCompositions.tempData.get(i).getNumber();
				for(AminoAcid aa : aaSet)
				{
					Composition aaComp = aa.getComposition(); 
					int prevNumber = number-aaComp.getNumber();
					if(this.isSet(prevNumber) && !intermediateCompositions.isSet(prevNumber))
					{
						intermediateCompositions.setAndAddIfNotExist(prevNumber);
					}
				}
			}
			if(end == intermediateCompositions.size())
				break;
			start = end;
		}
		
		Collections.sort(intermediateCompositions.tempData);
		return intermediateCompositions.tempData;
	}
	
	public boolean contains(Composition node) {
		return isSet(node.number);
	}
	
	private boolean isSet(int number)
	{
		int index = (number & indexMask) >>> 5;
		int offset = number & offsetMask;
		if((map[index] & (1 << offset)) == 0)
			return false;
		else
			return true;
	}
	
	protected void set(int number)
	{
		int index = (number & indexMask) >>> 5;
		int offset = number & offsetMask;
		map[index] |= (1 << offset);
	}
	
	protected void clear(int number)
	{
		int index = (number & indexMask) >>> 5;
		int offset = number & offsetMask;
		map[index] &= ~(1 << offset);
	}
	
	protected void add(int number)
	{
		tempData.add(new Composition(number)); 
	}
	
	private void setAndAddIfNotExist(int number)
	{
		int index = (number & indexMask) >>> 5;
		int offset = number & offsetMask;
		if((map[index] & (1 << offset)) == 0)	// nonexistant 
		{
			map[index] |= (1 << offset);	// set
			tempData.add(new Composition(number));	// add
		}
	}
	
	private CompositionFactory finalizeCompositionSet()	
	{ 
		if(tempData != null) Collections.sort(tempData);
		data = new int[tempData.size()];
		for(int i=0; i<tempData.size(); i++)
			data[i] = tempData.get(i).getNumber();
		tempData = null;
		return this;
	}
	

	protected void makeAllPossibleMasses()
	{
		setAndAddIfNotExist(0);

		Composition[] aaComposition = new Composition[aaSet.size()];
		int index=0;
		for(AminoAcid aa : aaSet)
			aaComposition[index++] = aa.getComposition();

		int start = 0;
		for(int l=0; l<maxLength; l++)
		{
			int end = tempData.size();
			for(int i=start; i<end; i++)
			{
				for(int j=0; j<aaComposition.length; j++)
					setAndAddIfNotExist(tempData.get(i).getNumber()+aaComposition[j].getNumber());
			}
			start = end;
		}
		finalizeCompositionSet();
	}	

	public static void main(String[] argv)
	{
	}
}