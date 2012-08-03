package edu.ucsd.msjava.msgf;

import java.util.ArrayList;

import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Constants;
import edu.ucsd.msjava.msutil.Enzyme;

public class NominalMassFactory extends MassFactory<NominalMass> {
	private float rescalingConstant = Constants.INTEGER_MASS_SCALER;
	private NominalMass[] factory;
	private NominalMass zero;
	
	public NominalMassFactory(AminoAcidSet aaSet, Enzyme enzyme, int maxLength)
	{
		super(aaSet, enzyme, maxLength);
		int heaviestNominalMass = aaSet.getHeaviestAA().getNominalMass();
		int maxIndex = heaviestNominalMass*maxLength;
		factory = new NominalMass[maxIndex+2];
		zero = factory[0] = new NominalMass(0);
		makeAllPossibleMasses(true);
	}

	private NominalMassFactory(int maxLength)
	{
		super(null, null, maxLength);
	}
	
	public NominalMass getInstance(float mass)
	{
		int massIndex = getMassIndex(mass);
		return getInstanceOfIndex(massIndex);
	}

	public float getRescalingConstant()	{ return rescalingConstant; }

	// returns instance exists in the factory
	public NominalMass getInstanceOfIndex(int index)
	{
		if(index < factory.length)
		{
			return factory[index];
		}
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
	
	public ArrayList<DeNovoGraph.Edge<NominalMass>> getEdges(NominalMass curNode)
	{
		return edgeMap.get(curNode);
	}
	
	@Override
	public NominalMass getPreviousNode(NominalMass curNode, AminoAcid aa) 
	{
		int index = curNode.getNominalMass() - aa.getNominalMass();
		if(index < 0)
			return null;
		return factory[index];
	}
	
	public NominalMass getNextNode(NominalMass curNode, AminoAcid aa) {
		int index = curNode.getNominalMass() + aa.getNominalMass();
		if(factory[index] == null)
			factory[index] = new NominalMass(index);
		return factory[index];
	}

	public NominalMass getComplementNode(NominalMass srm, NominalMass pmNode) {
		int index = pmNode.getNominalMass() - srm.getNominalMass();
		if(factory[index] != null)
			return factory[index];
		else
			return new NominalMass(index);
	}

	public ArrayList<NominalMass> getNodes(float peptideMass, Tolerance tolerance) {
		ArrayList<NominalMass> nodes = new ArrayList<NominalMass>();
		float tolDa = tolerance.getToleranceAsDa(peptideMass);
		int minIndex = getMassIndex(peptideMass-tolDa);
		int maxIndex = getMassIndex(peptideMass+tolDa);
		for(int index = minIndex; index<=maxIndex; index++)
		{
			if(factory[index] != null)
				nodes.add(factory[index]);
			else
				nodes.add(new NominalMass(index));
		}
		return nodes;
	}

	public NominalMass getNode(float peptideMass) {
		int index = getMassIndex(peptideMass);
		if(factory[index] != null)
			return factory[index];
		else
			return new NominalMass(index);
	}
	
	@Override
	public NominalMass getZero() {
		return zero;
	}
	
	public boolean contains(NominalMass node) {
		int index = node.getNominalMass();
		if(index < 0 || index >= factory.length)
			return false;
		return factory[index] != null;
	}
	
//	// static methods
	private static NominalMassFactory defaultNominalMassFactory = new NominalMassFactory(50);
	public static NominalMass getInstanceFor(float mass)
	{
		return defaultNominalMassFactory.getInstance(mass);
	}
}

