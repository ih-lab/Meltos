package lineage;


public class ParametersVolatile {
	
	// Input type
	protected enum Format { VCF, SNV, SNV_WITH_PROFILE}
	protected Format INPUT_FORMAT = Format.SNV;
	protected boolean CP = false;
	
	// SNV group partitioning
	protected double MAX_ALLOWED_VAF;
	protected double MIN_VAF_PRESENT; 
	protected double MAX_VAF_ABSENT;
	protected int MIN_SNVS_PER_GROUP;  
	protected int MIN_ROBUST_SNVS_PER_GROUP;  
	protected double MIN_GROUP_PROFILE_SUPPORT;  
	protected double MIN_VAF_TARGET_RATIO_PER_SAMPLE;
	
	// Clusters
	/** Minimum size a sub-population cluster must have to be a considered a node in the network */
	protected int MIN_CLUSTER_SIZE;
	protected int MIN_PRIVATE_CLUSTER_SIZE; 
	protected double MIN_ROBUST_CLUSTER_SUPPORT;
	
	/** Maximum centroid difference up to which two clusters can be collapsed */
	protected double MAX_COLLAPSE_CLUSTER_DIFF; 
	
	// Constraint graph and spanning tree generation
	/** Maximum VAF (used for the root node) */
	protected double VAF_MAX;
	/** Error margin used for comparing VAF centroid values when adding edges in the network */
	protected double VAF_ERROR_MARGIN;	
	protected boolean ALL_EDGES;
	
	/** Stop tree search once this many valid trees have been found */
	protected int MAX_NUM_TREES;
	protected int MAX_NUM_GROW_CALLS;
	protected int NUM_TREES_FOR_CONSISTENCY_CHECK;
	//public ParametersVolatile{
		
	//}
}