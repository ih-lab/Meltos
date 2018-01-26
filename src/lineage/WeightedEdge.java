package lineage;



import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.function.ToLongFunction;

//Represents a potential directed edge to be added to the tree.
public class WeightedEdge implements Comparable<WeightedEdge>{
	private double weight;
	private TreeNode parent;
	private TreeNode child;
	
	public WeightedEdge(TreeNode p, TreeNode c, double w){
		weight = w;
		parent = p;
		child = c;
	}
	public double getWeight(){
		return weight;
	}
	public TreeNode getParent(){
		return parent;
	}
	public TreeNode getChild(){
		return child;
	}
	public int getParentId(){
		return parent.getNodeId();
	}
	public int getChildId(){
		return child.getNodeId();
	}
	

	@Override
	public int compareTo(WeightedEdge arg0) {
		if(getWeight()<arg0.getWeight())
			return 1;
		else
			return -1;
	}
	
	@Override
	public boolean equals(Object e){
		if(e instanceof WeightedEdge){
		    WeightedEdge edge = (WeightedEdge) e;
			return (getParentId()==edge.getParentId() && getChildId()==edge.getChildId());
		
		} else
			return false;
	}
}
	
	
