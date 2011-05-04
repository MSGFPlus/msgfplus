package msgap;

import java.util.ArrayList;
import java.util.HashMap;

import msgf.DeNovoGraph;
import msgf.Profile;
import msgf.ProfilePeak;
import msutil.Matter;

public class ProfileNormalizerTest <T extends Matter>{

	private HashMap<T, Float> normalizerTable = null;
	
	private void updateNormalizerTable(DeNovoGraph<T> graph){
		
		normalizerTable = new HashMap<T, Float>();
		HashMap<T, Float> forwardTable = new HashMap<T, Float>();
		HashMap<T, Float> backwardTable = new HashMap<T, Float>();
		
		forwardTable.put(graph.getSource(), 1f);
		
		ArrayList<T> allNodes = new ArrayList<T>(graph.getIntermediateNodeList());
		allNodes.addAll(graph.getSinkList());
		
		for(int i=1; i<allNodes.size(); i++) // forward
		{
			T curNode = allNodes.get(i);
			float value = 0;
			
			// modified by Sangtae
			for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode)){
				T prevNode = edge.getPrevNode();
				if(forwardTable.containsKey(prevNode)){
					value += forwardTable.get(prevNode);
				}
			}
			
			if(value > 0) forwardTable.put(curNode, value);
		}
		
		for(T sink : graph.getSinkList()){
			backwardTable.put(sink, 1f);
		}
		
		for(int i=allNodes.size()-1; i>0; i--){ // backward
			T curNode = allNodes.get(i);
			if(!forwardTable.containsKey(curNode)) continue;
			
			// modified by Sangtae
			for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode)){
				T prevNode = edge.getPrevNode();
				if(backwardTable.containsKey(curNode)){
					float value = backwardTable.get(curNode);
					if(backwardTable.containsKey(prevNode)){
						value += backwardTable.get(prevNode);
					}
					backwardTable.put(prevNode, value);
						
				}
			}
		}
		
		float sourceValue = backwardTable.get(graph.getSource());
		
		for(int i=1; i<allNodes.size(); i++){
			T curNode = allNodes.get(i);
			if(backwardTable.containsKey(curNode) && forwardTable.containsKey(curNode)){
				float value = backwardTable.get(curNode) * forwardTable.get(curNode);
				normalizerTable.put(curNode, sourceValue / value);
			}
		}
	}
	
	public void normalize(Profile<T> profile, DeNovoGraph<T> graph){
		if(normalizerTable == null) updateNormalizerTable(graph);
	//	System.out.println(normalizerTable);
		for(ProfilePeak<T> peak : profile){
			float probability = peak.getProbability()*normalizerTable.get(peak.getNode());
			peak.setProbability(probability);
		}
	}
}
