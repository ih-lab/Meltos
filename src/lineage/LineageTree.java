package lineage;


import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class LineageTree implements Comparable<LineageTree> {
	protected HashMap<Integer, TreeNode> nodes;
	protected HashMap<Integer, ArrayList<Integer>> edges;
	
	//protected HashMap<Integer, >
	
	protected int numSamples;
	protected final static double VAF_ERROR_MARGIN = 0.1;
	
	public LineageTree(int numInputSamples) {
		numSamples = numInputSamples;
		nodes = new HashMap<Integer, TreeNode>();
		edges = new HashMap<Integer, ArrayList<Integer>>();
	}
	
	public void addNode(TreeNode n) {
		if(!nodes.containsKey(n.getNodeId())) {
			nodes.put(n.getNodeId(), n);
		}
	}
	
	public void addEdge(Integer from, Integer to) {
		if(!edges.containsKey(from)) {
			edges.put(from, new ArrayList<Integer>());
		}
		if(!edges.get(from).contains(to)) {
			edges.get(from).add(to);
		}
	}
	
	public void removeEdge(Integer from, Integer to) {
		ArrayList<Integer> nbrs = edges.get(from);
		if(nbrs != null) {
			for(Integer n : nbrs) {
				if(n.equals(to)) {
					nbrs.remove(n);
					break;
				}
			}
		}	
	}
	
	public boolean containsEdge(Integer from, Integer to) {
		if(edges.get(from) == null) return false;
		for(Integer n : edges.get(from)) {
			if(n.equals(to)) {
				return true;
			}
		}
		return false;
	}
	
	public TreeNode getRoot() {
		return nodes.get(0);
	}
	
	public boolean checkConstraint(TreeNode n) {
		ArrayList<Integer> nbrs = edges.get(n.getNodeId());			
		if(nbrs == null) return true;
				
		for(int i = 0; i < numSamples; i++) {
			double affSum = 0;
			double errMargin = 0.0;
			for(Integer n2 : nbrs) {
				affSum += nodes.get(n2).getVAF(i);
			}
			errMargin = VAF_ERROR_MARGIN;
			if(affSum > n.getVAF(i) + errMargin) {
				return false;
			}
		}
		return true;
	}
	
	public String toString() {
		String graph = "";
		for(Integer n1 : edges.keySet()) {
			ArrayList<Integer> nbrs = edges.get(n1);
			for(Integer n2 : nbrs) {
				graph += n1 + " -> " + n2 + "\n";
			}
		}
		return graph;
	}
	
	/** 
	 * Returns the error score associated with the tree, 
	 * which is the sqrt of the sum of the children AAF sum deviation from the parent AAF
	 */
	public double getErrorScore() {
		return computeErrorScore();
	}
	
	public double computeErrorScore() {
		ArrayList<Integer> nodeIds = new ArrayList<Integer>(edges.keySet());
		Collections.sort(nodeIds);
		
		double err = 0;
		for(Integer n : nodeIds) {
			ArrayList<Integer> nbrs = edges.get(n);			
			for(int i = 0; i < numSamples; i++) {
				double affSum = 0;
				for(Integer n2 : nbrs) {
					affSum += nodes.get(n2).getVAF(i);
				}
				if(affSum > nodes.get(n).getVAF(i)) {
					err += Math.pow(affSum - nodes.get(n).getVAF(i), 2);
				}
			}
		}
		return Math.sqrt(err);
	}
	
	public int compareTo(LineageTree t) {
		return new Double(this.getErrorScore()).compareTo(t.getErrorScore());
	}
}	
