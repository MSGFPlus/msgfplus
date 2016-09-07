package edu.ucsd.msjava.msgf;

import java.util.ArrayList;
import java.util.Collection;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Sequence;

public interface DeNovoNodeFactory<T extends Matter> {
    AminoAcidSet getAASet();

    T getZero();

    ArrayList<T> getNodes(float mass, Tolerance tolerance);

    T getNode(float mass);    // get the closest node from the mass

    T getComplementNode(T srm, T pmNode);

    ArrayList<T> getLinkedNodeList(Collection<T> destNodes);

    ArrayList<DeNovoGraph.Edge<T>> getEdges(T curNode);

    DeNovoGraph.Edge<T> getEdge(T curNode, T prevNode);

    Sequence<T> toCumulativeSequence(boolean isPrefix, Peptide pep);

    T getPreviousNode(T curNode, AminoAcid aa);

    T getNextNode(T curNode, AminoAcid aa);

    int size();

    boolean contains(T node);

    boolean isReverse();

    Enzyme getEnzyme();
}
