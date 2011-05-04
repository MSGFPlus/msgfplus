package suffixgraph.misc;

import java.util.HashMap;
import java.util.HashSet;

import suffixgraph.graphs.AbstractSuffixGraph;
import suffixgraph.nodes.AbstractNode;

public class GraphvizGraph {
  
  final static int NORMAL_NODE = compressColor(0, 500, 0);
  
  HashMap<Integer, Integer> nodeColors;
  HashMap<Long, Integer> normalEdgeMultiplicity;
  HashMap<Long, HashSet<String>> edgeLabels;
  
  public GraphvizGraph() {
    nodeColors = new HashMap<Integer, Integer>();
    normalEdgeMultiplicity = new HashMap<Long, Integer>();
    edgeLabels = new HashMap<Long, HashSet<String>>();
  }
  
  // 10 bit colors (1024) steps
  static private int compressColor(int r, int g, int b) {
    r = Math.min(r, 1023); b = Math.min(b, 1023); g = Math.min(g, 1023);
    return (r<<20) | (g<<10) | b;
  }
 
  /*
  static private String rgbColorString(int color) {
    return String.format("#%02x%02x%02x", (color>>>20)/4, ((color>>>10)&0x3FF)/4, (color&0x3FF)/4);  
  }
  
  static private int normalizeCount(int count) {
    //return 1024-count/200;
    return 500;
  }
  */
  
  public void setNormalNode(int index) {
    nodeColors.put(index, NORMAL_NODE);
  }
  
  public void setNormalEdge(int source, int sink) {
    long edge = (((long)(source))<<32)|sink;
    if (!normalEdgeMultiplicity.containsKey(edge)) normalEdgeMultiplicity.put(edge, 1);
    else normalEdgeMultiplicity.put(edge, normalEdgeMultiplicity.get(edge)+1);
  }
  
  public void setEdgeLabel(int source, int sink, String label) {
    long edge = (((long)(source))<<32)|sink;
    if (!edgeLabels.containsKey(edge)) edgeLabels.put(edge, new HashSet<String>());
    edgeLabels.get(edge).add(label);
  }
  
  public String toString() {
    StringBuffer sb = new StringBuffer();
    sb.append("digraph G {\nsize=\"40,30\"; ratio = fill; edge [fontsize=100]\n");
    // print all the nodes
    for (int nodeId : nodeColors.keySet()) {
      sb.append("node" + nodeId + " [width=1,height=1,style=filled,color=\"forestgreen\",label=\"\"];\n");
    }

    // print the edges
    for (long edge : normalEdgeMultiplicity.keySet()) {
      int nodeId1 = (int)(edge>>>32);
      int nodeId2 = (int)edge;
      StringBuffer label = new StringBuffer();
      for (String item : edgeLabels.get(edge)) {
        label.append(item + "|");
      }
      String finalLabel = label.substring(0, label.length()-1).toString();
      sb.append("node" + nodeId1 + " -> node" + nodeId2 + " [style=\"setlinewidth(10)\",color=\"forestgreen\",label=\"" + finalLabel + "\"];\n");
    }
    sb.append("}");
    return sb.toString();
  }
  
  
  public String toString(AbstractSuffixGraph<? extends AbstractNode> csg) {
    StringBuffer sb = new StringBuffer();
    sb.append("digraph G {\nsize=\"40,30\"; ratio = fill; edge [fontsize=100]\n");
    // print all the nodes
    for (int nodeId : nodeColors.keySet()) {
      if (csg.getNodeAt(nodeId).getPositions().length==0)
        sb.append("node" + nodeId + " [width=1,height=1,style=filled,color=\"forestgreen\",label=\"\"];\n");
      else
        sb.append("node" + nodeId + " [shape=box,width=1,height=1,style=filled,color=\"forestgreen\",label=\"\"];\n");
    }

    // print the edges
    for (long edge : normalEdgeMultiplicity.keySet()) {
      int nodeId1 = (int)(edge>>>32);
      int nodeId2 = (int)edge;
      StringBuffer label = new StringBuffer();
      for (String item : edgeLabels.get(edge)) {
        label.append(item + "|");
      }
      String finalLabel = label.substring(0, label.length()-1).toString();
      sb.append("node" + nodeId1 + " -> node" + nodeId2 + " [style=\"setlinewidth(10)\",color=\"forestgreen\",label=\"" + finalLabel + "\"];\n");
    }
    sb.append("}");
    return sb.toString();
  }
  
  public String toString(AbstractSuffixGraph<? extends AbstractNode> cgsg, AbstractSuffixGraph<? extends AbstractNode> csg) {
    
    
    StringBuffer sb = new StringBuffer();
    sb.append("digraph G {\nsize=\"40,30\"; ratio = fill;\n");
    // print all the normal nodes
    for (int nodeId : nodeColors.keySet()) {
      if (cgsg.getNodeAt(nodeId).getPositions().length==0)
        sb.append("node" + nodeId + " [width=1,height=1,style=filled,color=\"forestgreen\",label=\"\"];\n");
      else
        sb.append("node" + nodeId + " [shape=box,width=1,height=1,style=filled,color=\"forestgreen\",label=\"\"];\n");
      for (int j=0; j<cgsg.getNodeAt(nodeId).getEdgeCount(); j++) {
        int sinkNode = cgsg.getNodeAt(nodeId).getNodeIdAt(j);
        long edge = ((long)nodeId)<<32 | sinkNode;
        if (!normalEdgeMultiplicity.containsKey(edge)) {
          if (nodeColors.containsKey(sinkNode)) {
            sb.append("node" + nodeId + " -> node" + sinkNode + " [style=\"setlinewidth(10)\",color=\"forestgreen\"];\n");
          }
          else {
            sb.append("node" + nodeId + " -> node" + sinkNode + " [style=\"setlinewidth(10)\",color=\"black\"];\n");
          }
        }
      }
    }
    
    for (int i=csg.size(); i<cgsg.size(); i++) {
      if (cgsg.getNodeAt(i).getPositions().length==0)
        sb.append("node" + i + " [width=1,height=1,style=filled,color=\"red\",label=\"\"];\n");
      else 
        sb.append("node" + i + " [width=1,height=1,shape=box,style=filled,color=\"red\",label=\"\"];\n");
      for (int j=0; j<cgsg.getNodeAt(i).getEdgeCount(); j++) {
        int sinkNode = cgsg.getNodeAt(i).getNodeIdAt(j);
        long edge = ((long)i)<<32 | sinkNode;
        if (!normalEdgeMultiplicity.containsKey(edge)) {
          if (!nodeColors.containsKey(sinkNode)) {
            sb.append("node" + i + " -> node" + sinkNode + " [style=\"setlinewidth(10)\",color=\"red\"];\n");
          }
          else {
            sb.append("node" + i + " -> node" + sinkNode + " [style=\"setlinewidth(10)\",color=\"blue\"];\n");
          }
        }
      }
    }

    // print the edges
    for (long edge : normalEdgeMultiplicity.keySet()) {
      int nodeId1 = (int)(edge>>>32);
      int nodeId2 = (int)edge;
      StringBuffer label = new StringBuffer();
      for (String item : edgeLabels.get(edge)) {
        label.append(item + "|");
      }
      String finalLabel = label.substring(0, label.length()-1).toString();
      sb.append("node" + nodeId1 + " -> node" + nodeId2 + " [style=\"setlinewidth(10)\",color=\"forestgreen\",label=\"" + finalLabel + "\"];\n");
    }
    sb.append("}");
    return sb.toString();
  }
  
  
  
}