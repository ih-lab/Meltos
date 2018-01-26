package lineage;


import java.util.ArrayList;
import java.util.Comparator;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.function.ToLongFunction;


public class WeightedEdgeSystem implements Comparable<WeightedEdgeSystem>{

	private double fullWeight;
	//private int type;//(1 is old node->new node, 2 is old node->new node->old node, 3 is a more complex
	ArrayList<WeightedEdge> toBreak;
	ArrayList<WeightedEdge> toMake;
	private TreeNode newNode;
	private TreeNode parentNode;
	public WeightedEdgeSystem(TreeNode kernel, TreeNode top, ArrayList<WeightedEdge> broken, ArrayList<WeightedEdge> made){
		toBreak = broken;
		toMake = made;
		newNode = kernel;
		parentNode = top;
		fullWeight = 0;
		for(WeightedEdge w: toMake)
			fullWeight+= w.getWeight();
		
		for(WeightedEdge w: toBreak)
			fullWeight-=w.getWeight();
		//if(made.size()==0)
		//	System.out.println("WHY??!");
		//type = kind;
	}

	
	public double getWeight(){
		return fullWeight;
	}
	
	public TreeNode getNewNode(){
		return newNode;
	}
	
	public TreeNode getParentNode(){
		return parentNode;
	}

	@Override
	public int compareTo(WeightedEdgeSystem arg0) {
		if(getWeight()<arg0.getWeight())
			return 1;
		else
			return -1;
	}
	
	/**
	
	boolean isType(String kind){
		return kind.equals(type);
	}
	**/
	
	
}
