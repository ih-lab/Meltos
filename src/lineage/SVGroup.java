package lineage;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import lineage.SVReClusterer.Cluster;
import lineage.SVReClusterer.DistanceMetric;


//NOTE: Many of the comments here are relics of when it was LICHeE code at the moment. I'll fix it later. They should still be accurate if you replace SSNV with SV.


/**
 * An SNV group is a set of SNVs occurring in a given subset of samples.
 * All SNVs are partitioned into SNV groups based on their occurrence across samples.
 * Given S samples, there can be at most 2^S different groups.
 * An SNV group is uniquely identified by an S-bit binary tag (each bit corresponding to
 * a given sample), where a bit is set if that sample contains the SNVs in this group.
 *
 * @autor viq
 */
public class SVGroup implements Serializable {

	private static final long serialVersionUID = 1L;

	/** Binary tag identifying the group 
	 * (the length of the tag is equal to the number of input samples) */
	private String tag;
	
	/** Number of samples represented by this group */
	private int numSamples;
	
	/** Indices of the samples represented in this group (from 0 to |tag|-1 MSF order) */
	private int[] sampleIndex;
	
	/** Alternative allele frequency data matrix (numSNVs x numSamples) */
	protected transient double[][] alleleFreqBySample;
	
	/** SubPopulation clusters */
	Cluster[] subPopulations;
	
	/** SNVs assigned to this group */
	private transient ArrayList<SVEntry> svs;
	
	/** Flag indicating whether this group is robust */
	private transient boolean isRobust;
	
	public transient boolean isOriginal;
	
	//private static Logger logger = LineageEngine.logger;
	
	public SVGroup(String groupTag, ArrayList<SVEntry> groupSNVs, boolean isGroupRobust) {
		tag = groupTag;
		isRobust = isGroupRobust;
		numSamples = 0;		
		sampleIndex = new int[tag.length()];
		for(int i = 0; i < tag.length(); i++) {
			if(tag.charAt(i) == '1') {
				sampleIndex[numSamples] = i;
				numSamples++;
			}
		}
		svs = groupSNVs;
		alleleFreqBySample = new double[svs.size()][numSamples];
		for(int i = 0; i < svs.size(); i++) {
			SVEntry snv = svs.get(i);
			for(int j = 0; j < numSamples; j++) {
				alleleFreqBySample[i][j] = snv.getVAF(sampleIndex[j]);
			}
		}
		//if(groupTag.equals("00010100"))
			//System.out.println("test");
		
		isOriginal=false;
	}
	
	public void setOriginal(){
		isOriginal=true;
	}

	public SVGroup(String groupTag, double[] centroid, int size) {
		tag = groupTag;
		isRobust = true;
		numSamples = 0;		
		sampleIndex = new int[tag.length()];
		
		for(int i = 0; i < tag.length(); i++) {
			if(tag.charAt(i) == '1') {
				sampleIndex[numSamples] = i;
				numSamples++;
			}
		}
		double[] c = new double[numSamples];
		int idx = 0;
		for(int i = 0; i < tag.length(); i++) {
			if(tag.charAt(i) == '1') {
				c[idx] = centroid[i];
				idx++;
			}
		}
		
		svs = new ArrayList<SVEntry>();
		alleleFreqBySample = new double[svs.size()][numSamples];
		subPopulations = new Cluster[1];
		SVReClusterer aafc = new SVReClusterer();
		subPopulations[0] = aafc.new Cluster(c, 0);
		
		
	}
	
	// Getters/Setters
	
	public ArrayList<SVEntry> getSVs() {
		return svs;
	}
	
	public double[][] getAlleleFreqBySample() {
		return alleleFreqBySample;
	}
	
	public int getNumSamples() {
		return numSamples;
	}
	
	public int getNumSamplesTotal() {
		return tag.length();
	}
	
	public int[] getSampleIds() {
		return sampleIndex;
	}
	
	public int getNumSVs() {
		return svs.size();
	}
	
	public Cluster[] getSubPopulations() {
		return subPopulations;
	}
	
	public String getTag() {
		return tag;
	}
	
	public boolean isRobust() {
		return isRobust;
	}
	
	/**
	 * Returns the index of this sample in the centroid/AAF data of the group
	 * @return -1 if this sample is not represented in the group
	 */
	public int getSampleIndex(int sampleId) {
		for(int i = 0; i < numSamples; i++) {
			if(sampleIndex[i] == sampleId) {
				return i;
			}
		}
		return -1;
	}
	
	/**
	 * Returns true if the given sample contains the mutations of this group
	 */
	public boolean containsSample(int sampleId) {
		return (getSampleIndex(sampleId) != -1);
	}
	
	public boolean equals(Object o) {
		if(!(o instanceof SVGroup)) {
			return false;
		}
		SVGroup g = (SVGroup) o;
		if(this.tag == g.tag) {
			return true;
		} else {
			return false;
		}
	}
	
	// --- Sub-population Cluster Filtering / Collapse ---
	
	/** Entry in the cluster centroid distance minimum priority queue */
	protected class ClusterPairDistance {
		/** Cluster pair */
		protected int clusterId1;
		protected int clusterId2;
		/** Distance between cluster centroids */
		protected double distance;
		
		public ClusterPairDistance(int cluster1, int cluster2, double clusterDistance) {
			clusterId1 = cluster1;
			clusterId2 = cluster2;
			distance = clusterDistance;
		}
	}
	
	
	/**
	 * Set the sub-populations of this group based on clustering results
	 * Performs filtering based on cluster size, as well as collapses clusters
	 * with similar centroids
	 * @param clusters - clustering algorithm results
	 */
	public void setSubPopulations(Cluster[] clusters) {
		
		// 1. filter out clusters that are too small 
		ArrayList<Cluster> filteredClusters = new ArrayList<Cluster>();
		for(Cluster c : clusters) {
			//if((c.getMembership().size() >= Parameters.MIN_CLUSTER_SIZE) || ((numSamples == 1) && c.getMembership().size() >= Parameters.MIN_PRIVATE_CLUSTER_SIZE)) { // don't filter out private mutations
			
			filteredClusters.add(c);
			//} else {
				//logger.log(Level.INFO, "**Filtered due to cluster size constraint (" + tag + " size " + c.getMembership().size() + "):");
			//	for(Integer snv : c.getMembership()) {
			//		SVEntry entry = svs.get(snv);
					//logger.log(Level.INFO, entry.toString());
			//	}
			//}
		}
		
		if(filteredClusters.size() < 1) {
			//logger.log(Level.WARNING, "All clusters in group " + tag + " have been filtered out. The group is now empty.");
			subPopulations = new Cluster[filteredClusters.size()];
			subPopulations = filteredClusters.toArray(subPopulations);
			return;
		}
		
		// 2. collapse clusters that have similar centroids
		
		// compute the distance matrix between clusters
		// as long as there are clusters to collapse (i.e. cluster centroid distance
		// is less than MAX_COLLAPSE_CLUSTER_DIFF), collapse clusters with smallest distance first
		ArrayList<ClusterPairDistance> minDistQueue = new ArrayList<ClusterPairDistance>();
		int numClusters = filteredClusters.size();
		for(int i = 0; i < numClusters; i++) {
			for(int j = i+1; j < numClusters; j++) {
				Cluster c1 = filteredClusters.get(i);
				Cluster c2 = filteredClusters.get(j);
				//double dist = c1.getDistanceToCluster(c2.getCentroid(), DistanceMetric.EUCLIDEAN);
				double dist = c1.getDistanceToCluster(c2.getCentroid(), DistanceMetric.AVG_PER_SAMPLE);
				ClusterPairDistance pd = new ClusterPairDistance(c1.getId(), c2.getId(), dist);
				
				int k = 0;
				for(k = 0; k < minDistQueue.size(); k++) {
					if(minDistQueue.get(k).distance > dist) {
						break;
					}
				}
				minDistQueue.add(k, pd);
			}
		}
		
		while((minDistQueue.size() > 0) && 
				((minDistQueue.get(0).distance < Parameters.MAX_COLLAPSE_CLUSTER_DIFF))) {
			ClusterPairDistance pd = minDistQueue.remove(0);
			Cluster c1 = clusters[pd.clusterId1];
			Cluster c2 = clusters[pd.clusterId2];
			
			// collapse into c1
			if(pd.clusterId1<pd.clusterId2)
				for(Integer obs : c2.getMembership()) {
					c1.addMember(obs);
					svs.get(obs).setCluster(pd.clusterId1);
				}
			else{
				for(Integer obs : c1.getMembership()) {
					c2.addMember(obs);
					svs.get(obs).setCluster(pd.clusterId2);
				}
			}
			numClusters--;
			
			c1.recomputeCentroidAndStdDev(alleleFreqBySample, svs.size(), numSamples);
			filteredClusters.remove(c2);
			//logger.log(Level.FINE, "Collapse clusters: group = " + tag + " cluster " + pd.clusterId1 + " and " + pd.clusterId2 + 
			//		" distance = " + pd.distance + ". New cluster: " + c1);
			
			// remove distances from c1 and c2 from the queue
			ArrayList<ClusterPairDistance> toRemove = new ArrayList<ClusterPairDistance>();
			for(ClusterPairDistance cpd : minDistQueue) {
				if(cpd.clusterId1 == c1.getId() || cpd.clusterId2 == c1.getId() 
				   || cpd.clusterId1 == c2.getId() || cpd.clusterId2 == c2.getId()) {
					toRemove.add(cpd);
				}
			}
			
			//Trying to fix indexing problems later on
			for(SVEntry svE: svs){
				if(svE.clusterID>pd.clusterId2)
					svE.setCluster(svE.clusterID-1);
			}
			
			minDistQueue.removeAll(toRemove);
			
			// compute the distance from c1 to all the other clusters
			for(Cluster c : filteredClusters) {
				if(c.getId() == c1.getId() || c.getId() == c2.getId()) {
					continue;
				}
				double dist = c1.getDistanceToCluster(c.getCentroid(), DistanceMetric.AVG_PER_SAMPLE);
				ClusterPairDistance cpd = new ClusterPairDistance(c1.getId(), c.getId(), dist);
				
				int k = 0;
				for(k = 0; k < minDistQueue.size(); k++) {
					if(minDistQueue.get(k).distance > dist) {
						break;
					}
				}
				minDistQueue.add(k, cpd);
			}
		}
		
		// compute the robustness of each cluster
		for(Cluster c : filteredClusters) {
			int numRobust = 0;
			for(Integer snv : c.getMembership()) {
				SVEntry entry = svs.get(snv);
				if(entry.isRobust()) {
					numRobust++;
					if(numRobust >= Parameters.MIN_ROBUST_CLUSTER_SUPPORT) {
						c.setRobust();
						break;
					}
				}
			}
		}
		
		subPopulations = new Cluster[filteredClusters.size()];
		subPopulations = filteredClusters.toArray(subPopulations);
	}
	
	public String toString() {
		String group = "";
		group += "(tag = " + this.tag + ", ";
		group += "numSamples = " + this.numSamples + ", ";
		group += "numSVs = " + this.svs.size() + ") ";
		if(this.subPopulations != null) group += "numSubPopulations = " + this.subPopulations.length;
		return group;
	}
	
	public void removeCluster(Cluster c) {
		if(subPopulations == null || subPopulations.length == 0) {
			return;
		}
		Cluster[] clusters = new Cluster[subPopulations.length - 1];
		int j = 0;
		for(int i = 0; i < subPopulations.length; i++) {
			if(!subPopulations[i].equals(c)) {
				clusters[j] = subPopulations[i];
				j++;
			} 
		}
		subPopulations = clusters;
	}
	
	public void addCluster(Cluster c) {
		if(subPopulations == null) {
			return;
		}
		Cluster[] clusters = new Cluster[subPopulations.length + 1];
		int j = 0;
		for(int i = 0; i < subPopulations.length; i++) {
			clusters[j] = subPopulations[i];
			j++; 
		}
		clusters[subPopulations.length] = c;
		subPopulations = clusters;
	}
	
	
	//Method for use in retrieving the other SVs assigned to the same cluster as a given node from the LICHeE tree input.
	public ArrayList<SVEntry> retrieveClusterEntries(SVEntry oldEntry){
		Cluster relevantCluster = subPopulations[oldEntry.getCluster()];
		ArrayList<SVEntry> otherEntries = new ArrayList<SVEntry>();
		for(SVEntry sv: svs)
			if(sv!=oldEntry)
				otherEntries.add(sv);
		return otherEntries;
		
		
	}
}
