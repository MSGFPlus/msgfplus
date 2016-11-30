package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.msutil.Modification.Location;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class FlexAminoAcidGraph extends DeNovoGraph<NominalMass> {
    public static final int MODIFIED_EDGE_PENALTY = 0;
    private ScoredSpectrum<NominalMass> scoredSpec;
    private Enzyme enzyme;
    private boolean direction;    // true: forward (e.g. Lys-C), false: reverse (e.g. Trypsin)
    private AminoAcidSet aaSet;
    private boolean useProtNTerm;
    private boolean useProtCTerm;

    private HashMap<NominalMass, ArrayList<DeNovoGraph.Edge<NominalMass>>> edgeMap;
    private HashMap<NominalMass, Integer> nodeScore;

    private static AtomicInteger negativeCompNodeMassWarnCount;
    private static AtomicInteger negativeNodeMassWarnCount;

    private static AtomicInteger getNodeScoreNullNodeCount;
    private static AtomicInteger getNodeScoreExceptionCount;

    public FlexAminoAcidGraph(
            AminoAcidSet aaSet,
            int peptideMass,
            Enzyme enzyme,
            ScoredSpectrum<NominalMass> scoredSpec
    ) {
        this(aaSet, peptideMass, enzyme, scoredSpec, false, false);
    }

    public FlexAminoAcidGraph(
            AminoAcidSet aaSet,
            int peptideMass,
            Enzyme enzyme,
            ScoredSpectrum<NominalMass> scoredSpec,
            boolean useProteinNTerm,
            boolean useProteinCTerm
    ) {
        this.enzyme = enzyme;
        this.direction = scoredSpec.getMainIonDirection();
        this.scoredSpec = scoredSpec;
        this.aaSet = aaSet;
        this.useProtNTerm = useProteinNTerm;
        this.useProtCTerm = useProteinCTerm;

        super.source = new NominalMass(0);

        super.pmNode = new NominalMass(peptideMass);

        if (negativeNodeMassWarnCount == null) {
            negativeNodeMassWarnCount = new AtomicInteger();
        }

        if (negativeCompNodeMassWarnCount == null) {
            negativeCompNodeMassWarnCount = new AtomicInteger();
        }

        if (getNodeScoreNullNodeCount == null) {
            getNodeScoreNullNodeCount = new AtomicInteger();
        }

        if (getNodeScoreExceptionCount == null) {
            getNodeScoreExceptionCount = new AtomicInteger();
        }

        edgeMap = new HashMap<NominalMass, ArrayList<DeNovoGraph.Edge<NominalMass>>>();
        edgeMap.put(source, new ArrayList<DeNovoGraph.Edge<NominalMass>>());
        setForwardEdgesFromSource();
        setForwardEdgesFromIntermediateNodes();
        super.intermediateNodes = new ArrayList<NominalMass>(edgeMap.keySet());
        Collections.sort(super.intermediateNodes);
        this.setBackwardEdgesFromSink();

        super.sinkNodes = new ArrayList<NominalMass>();
        sinkNodes.add(pmNode);

        computeNodeScores();
    }

    @Override
    public NominalMass getComplementNode(NominalMass node) {
        return new NominalMass(pmNode.getNominalMass() - node.getNominalMass());
    }

    @Override
    public ArrayList<DeNovoGraph.Edge<NominalMass>> getEdges(NominalMass curNode) {

        return edgeMap.get(curNode);
    }

    @Override
    public int getNodeScore(NominalMass node) {

        if (node == null) {
            int errorCount = getNodeScoreNullNodeCount.addAndGet(1);
            if (notifyError(errorCount)) {
                System.out.println("Note: null node encountered in getNodeScore");
            }
            return 0;
        }

        try {
            return nodeScore.get(node);
        } catch (Exception ex) {
            int errorCount = getNodeScoreExceptionCount.addAndGet(1);
            if (notifyError(errorCount)) {
                System.out.println("Note: Exception in getNodeScore retrieving node at nominal mass " +
                        node.getNominalMass() + ": " + ex.getMessage());
            }
            return 0;
        }

    }

    @Override
    public int getScore(Peptide pep) {
        int score = 0;

        NominalMass prevNode = source;
        int nominalMass = 0;
        for (int i = 0; i < pep.size() - 1; i++) {
            AminoAcid aa;
            if (direction == true)
                aa = pep.get(i);
            else
                aa = pep.get(pep.size() - 1 - i);

            nominalMass += aa.getNominalMass();
            NominalMass curNode = new NominalMass(nominalMass);
            int nodeScore = getNodeScore(curNode);
            int edgeScore = scoredSpec.getEdgeScore(curNode, prevNode, aa.getMass());
            if (prevNode == source && direction == false && enzyme != null) {
                if (enzyme.isCleavable(aa))
                    edgeScore += aaSet.getPeptideCleavageCredit();
                else
                    edgeScore += aaSet.getPeptideCleavagePenalty();
            }
            prevNode = curNode;
            score += nodeScore + edgeScore;
        }
        if (direction == true && enzyme != null) {
            if (enzyme.isCleavable(pep.get(pep.size() - 1)))
                score += aaSet.getPeptideCleavageCredit();
            else
                score += aaSet.getPeptideCleavagePenalty();
        }
        if (direction == true)
            nominalMass += pep.get(pep.size() - 1).getNominalMass();
        else
            nominalMass += pep.get(0).getNominalMass();

        if (nominalMass != pmNode.getNominalMass())
            return Integer.MIN_VALUE;
        else
            return score;
    }

    @Override
    public int getScore(Annotation annotation) {
        int score = getScore(annotation.getPeptide());
        if (enzyme != null) {
            AminoAcid neighboringAA;
            if (enzyme.isCTerm())
                neighboringAA = annotation.getPrevAA();
            else
                neighboringAA = annotation.getNextAA();
            if (neighboringAA == null || enzyme.isCleavable(neighboringAA))
                score += aaSet.getNeighboringAACleavageCredit();
            else
                score += aaSet.getNeighboringAACleavagePenalty();
        }
        return score;
    }

    @Override
    public boolean isReverse() {
        return !direction;
    }

    @Override
    public AminoAcidSet getAASet() {
        return aaSet;
    }

    private void computeNodeScores() {
        nodeScore = new HashMap<NominalMass, Integer>();
        nodeScore.put(source, 0);

        boolean warnNegativeNodeMass = false;
        boolean warnNegativeCompNodeMass = false;

        for (int i = 1; i < intermediateNodes.size(); i++) {
            NominalMass node = intermediateNodes.get(i);
            NominalMass compNode = this.getComplementNode(node);
            if (node.getNominalMass() < 0 && !warnNegativeNodeMass) {
                warnNegativeNodeMass = true;
                // Mass of the node is negative
                // This can happen if we have a negative dynamic mod at the C-terminus, for example Lys-Loss
                int warnCount = negativeNodeMassWarnCount.addAndGet(1);
                if (notifyError(warnCount)) {
                    System.out.println("Note: negative node mass in computeNodeScores " +
                            "(count = " + Integer.toString(warnCount) + ")");
                }
            }
            if (compNode.getNominalMass() < 0 && !warnNegativeCompNodeMass) {
                warnNegativeCompNodeMass = true;
                int warnCount = negativeCompNodeMassWarnCount.addAndGet(1);
                if (notifyError(warnCount)) {
                    System.out.println("Note: negative compnode mass in computeNodeScores " +
                            "(count = " + Integer.toString(warnCount) + ")");
                }
            }
            int score;
            if (isReverse())
                score = scoredSpec.getNodeScore(compNode, node);
            else
                score = scoredSpec.getNodeScore(node, compNode);
            nodeScore.put(node, score);
        }
        for (NominalMass node : this.sinkNodes)
            nodeScore.put(node, 0);
    }

    private boolean notifyError(int errorCount) {
        if (errorCount < 5 || errorCount == 100 || errorCount == 1000 || errorCount % 10000 == 0) {
            return true;
        } else {
            return false;
        }
    }

    private void setForwardEdgesFromSource() {
        Location location;
        if (direction) {
            if (!this.useProtNTerm)
                location = Location.N_Term;
            else
                location = Location.Protein_N_Term;
        } else {
            if (!this.useProtCTerm)
                location = Location.C_Term;
            else
                location = Location.Protein_C_Term;
        }

        ArrayList<AminoAcid> aaList = aaSet.getAAList(location);
        makeForwardEdges(source, aaList, enzyme != null && direction == enzyme.isNTerm());
    }

    private void setForwardEdgesFromIntermediateNodes() {
        ArrayList<AminoAcid> aaList = aaSet.getAAList(Location.Anywhere);
        for (int i = 1; i < pmNode.getNominalMass(); i++)
            makeForwardEdges(new NominalMass(i), aaList, false);
    }

    private void setBackwardEdgesFromSink() {
        Location location;
        if (direction) {
            if (!this.useProtCTerm)
                location = Location.C_Term;
            else
                location = Location.Protein_C_Term;
        } else {
            if (!this.useProtNTerm)
                location = Location.N_Term;
            else
                location = Location.Protein_N_Term;
        }

        ArrayList<AminoAcid> aaList = aaSet.getAAList(location);

//		if(enzymaticCleavageOnly && direction != enzyme.isNTerm())
//			aaList = aaSet.getEnzymeAAList();
//		else
//			aaList = aaSet.getAAList(location);

        int peptideNominalMass = pmNode.getNominalMass();
        ArrayList<DeNovoGraph.Edge<NominalMass>> edges = new ArrayList<DeNovoGraph.Edge<NominalMass>>();
        for (AminoAcid aa : aaList) {
            NominalMass prevNode = new NominalMass(peptideNominalMass - aa.getNominalMass());
            if (edgeMap.containsKey(prevNode)) {
                DeNovoGraph.Edge<NominalMass> edge = new DeNovoGraph.Edge<NominalMass>(prevNode, aa.getProbability(), aaSet.getIndex(aa), aa.getMass());
                edges.add(edge);
                if (enzyme != null && direction != enzyme.isNTerm()) {
                    if (enzyme.isCleavable(aa))
                        edge.setCleavageScore(aaSet.getPeptideCleavageCredit());
                    else
                        edge.setCleavageScore(aaSet.getPeptideCleavagePenalty());
                }
                if (aa.isModified())
                    edge.setErrorScore(MODIFIED_EDGE_PENALTY);
            }
        }
        edgeMap.put(pmNode, edges);
    }

    private void makeForwardEdges(NominalMass curNode, ArrayList<AminoAcid> aaList, boolean addCleavageScore) {
        if (edgeMap.get(curNode) == null)
            return;
        int curNominalMass = curNode.getNominalMass();
        for (AminoAcid aa : aaList) {
            int nextNodeNominalMass = curNominalMass + aa.getNominalMass();
            if (nextNodeNominalMass >= pmNode.getNominalMass())
                continue;
            NominalMass nextNode = new NominalMass(nextNodeNominalMass);
            ArrayList<DeNovoGraph.Edge<NominalMass>> edges = edgeMap.get(nextNode);
            if (edges == null) {
                edges = new ArrayList<DeNovoGraph.Edge<NominalMass>>();
                edgeMap.put(nextNode, edges);
            }

            DeNovoGraph.Edge<NominalMass> edge = new DeNovoGraph.Edge<NominalMass>(
                    curNode,
                    aa.getProbability(),
                    aaSet.getIndex(aa),
                    aa.getMass());
//			if(curNode.getNominalMass() == 57 && nextNode.getNominalMass() == 114)
//				System.out.println("Debug");
            int errorScore = scoredSpec.getEdgeScore(nextNode, curNode, aa.getMass());
//			if(aa.isModified())
//				errorScore += MODIFIED_EDGE_PENALTY;
            if (errorScore < -100 || errorScore > 100) {
                System.err.println("Invalid ErrorScore: " + errorScore);
                System.exit(-1);
            }
            edge.setErrorScore(errorScore);
            if (addCleavageScore) {
                if (enzyme.isCleavable(aa))
                    edge.setCleavageScore(aaSet.getPeptideCleavageCredit());
                else
                    edge.setCleavageScore(aaSet.getPeptideCleavagePenalty());
            }

            edges.add(edge);
        }
    }

//	private void computeEdgeScores()
//	{
//		for(NominalMass curNode : intermediateNodes)
//		{
//			ArrayList<DeNovoGraph.Edge<NominalMass>> edges = edgeMap.get(curNode);
//			for(DeNovoGraph.Edge<NominalMass> edge : edges)
//			{
//				NominalMass prevNode = edge.getPrevNode();
//				int errorScore = scoredSpec.getEdgeScore(curNode, prevNode, edge.getEdgeMass());
//				assert(errorScore == edge.getErrorScore());
//				edge.setErrorScore(errorScore);
//			}
//			edgeMap.put(curNode, edges);
//		}
//	}	
}
