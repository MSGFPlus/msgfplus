package msgf;

import java.util.ArrayList;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.Matter;

public class IntMassFactory extends MassFactory<IntMassFactory.IntMass> {
	private float rescalingConstant;
	private IntMass[] factory;
	private IntMass zero;
	private int[] aaMassIndex;
	
	public IntMassFactory(AminoAcidSet aaSet, Enzyme enzyme, int maxLength, float rescalingConstant, boolean preComputeEdges)
	{
		super(aaSet, enzyme, maxLength);
		this.rescalingConstant = rescalingConstant;
		int heaviestAAIndex = this.getMassIndex(aaSet.getHeaviestAA().getMass());
		int maxIndex = heaviestAAIndex*maxLength;
		factory = new IntMass[maxIndex+2];
		zero = factory[0] = new IntMass(0);
		aaMassIndex = new int[128];
		for(AminoAcid aa : aaSet)
			aaMassIndex[aa.getResidue()] = getMassIndex(aa.getMass());
		makeAllPossibleMasses(preComputeEdges);
	}

	public IntMassFactory(AminoAcidSet aaSet, Enzyme enzyme, int maxLength, float rescalingConstant)
	{
		this(aaSet, enzyme, maxLength, rescalingConstant, true);
	}
	
	public IntMass getInstance(float mass)
	{
		int massIndex = getMassIndex(mass);
		return getInstanceOfIndex(massIndex);
	}

	public float getRescalingConstant()	{ return rescalingConstant; }

	// returns instance exists in the factory
	public IntMass getInstanceOfIndex(int index)
	{
		if(index < factory.length)
			return factory[index];
		else
			return null;
	}
	
	public int getMassIndex(float mass)
	{
		return Math.round(mass*rescalingConstant);
	}
	
	public float getMassFromIndex(int massIndex)
	{
		return massIndex/rescalingConstant;
	}
	
	public ArrayList<DeNovoGraph.Edge<IntMass>> getEdges(IntMass curNode)
	{
		if(edgeMap != null)
			return edgeMap.get(curNode);
		int curIndex = curNode.massIndex;
		ArrayList<DeNovoGraph.Edge<IntMass>> edges = new ArrayList<DeNovoGraph.Edge<IntMass>>();
		for(AminoAcid aa : aaSet)
		{
			int prevIndex = curIndex - aaMassIndex[aa.getResidue()];
			IntMass prevNode = new IntMass(prevIndex);
			DeNovoGraph.Edge<IntMass> edge = new DeNovoGraph.Edge<IntMass>(prevNode, aa.getProbability(), aaSet.getIndex(aa), aa.getMass());
			int cleavageScore = 0;
			if(prevIndex == 0 && enzyme != null)
			{
				if(enzyme.isCleavable(aa))
					cleavageScore += aaSet.getPeptideCleavageCredit();
				else
					cleavageScore += aaSet.getPeptideCleavagePenalty();
			}
			edge.setCleavageScore(cleavageScore);
			edges.add(edge);
		}
		return edges;
	}
	
	@Override
	public DeNovoGraph.Edge<IntMass> getEdge(IntMass curNode, IntMass prevNode)
	{
		return null;
	}
	
	@Override
	public IntMass getPreviousNode(IntMass curNode, AminoAcid aa) {
		int index = curNode.getMassIndex() - getMassIndex(aa.getMass());
		if(index < 0)
			return null;
		else
			return factory[index];
	}
	
	public IntMass getNextNode(IntMass curNode, AminoAcid aa) {
		int index = curNode.getMassIndex() + getMassIndex(aa.getMass());
		if(factory[index] == null)
			factory[index] = new IntMass(index);
		return factory[index];
	}

	public IntMass getComplementNode(IntMass srm, IntMass pmNode) {
		int index = pmNode.massIndex - srm.massIndex;
		if(factory[index] != null)
			return factory[index];
		else
			return new IntMass(index);
	}

	public ArrayList<IntMass> getNodes(float peptideMass, Tolerance tolerance) {
		ArrayList<IntMass> nodes = new ArrayList<IntMass>();
		float tolDa = tolerance.getToleranceAsDa(peptideMass);
		int minIndex = getMassIndex(peptideMass-tolDa);
		int maxIndex = getMassIndex(peptideMass+tolDa);
		for(int index = minIndex; index<=maxIndex; index++)
		{
			if(factory[index] != null)
				nodes.add(factory[index]);
			else
				nodes.add(new IntMass(index));
		}
		return nodes;
	}

	public IntMass getNode(float peptideMass) {
		int index = getMassIndex(peptideMass);
		if(factory[index] != null)
			return factory[index];
		else
			return new IntMass(index);
	}
	
	public class IntMass extends Matter {
		private int massIndex;
		protected IntMass(int massIndex)
		{
			this.massIndex = massIndex;
		}
		@Override
		public float getMass() {
			return massIndex/rescalingConstant;
		}

		@Override
		public int getNominalMass() {
			return massIndex;
		}

		public int getMassIndex() {
			return massIndex;
		}
		
		@Override
		public int hashCode()	{ return massIndex; }
		
		@Override
		public boolean equals(Object obj)	
		{
			if(!(obj instanceof IntMass))
				return false;
			return (massIndex == ((IntMass)obj).massIndex);
		}
		
		@Override
		public String toString()
		{
			return String.valueOf(massIndex);
		}
	}

	@Override
	public IntMass getZero() {
		return zero;
	}

	public boolean contains(IntMass node) {
		int index = node.massIndex;
		if(index < 0 || index >= factory.length)
			return false;
		return factory[node.massIndex] != null;
	}
}

