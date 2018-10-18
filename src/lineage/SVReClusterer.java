package lineage;


import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import weka.clusterers.ClusterEvaluation;
import weka.clusterers.EM;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

/**
 * This class is a combination of a large number of different classes for accomplishing tasks from LICHeE's clustering algorithm, modified for 
 * the purposes of using clustering to assign SVs to SSNV nodes.
 */
public class SVReClusterer {

	/**
	 * Expectation Maximization clustering
	 * @param data - matrix of observations (numObs x numFeatures)
	 * @param k - number of clusters
	 */
	
	
	//Basic information required for the clustering and grouping processes
	private ArrayList<SVEntry> toHandle;
	private ArrayList<SVEntry> oldSSNVs;
	private ArrayList<SVEntry> combinedEntries;
	
	//Data from SNVDataStore:
	private int numSamples;
	/** Input sample names */
	private ArrayList<String> sampleNames;
	/** List of all input somatic SNVs */
	protected ArrayList<SVEntry> somaticSVs;
	/** List of SNVs with ambiguous VAFs in at least one sample */
	protected ArrayList<SVEntry> ambiguousSVs;
	
	HashMap<String, ArrayList<SVEntry>> tag2SVs;
	private HashMap<String, ArrayList<Cluster>> tag2Clusters;
	/** Index of the normal sample in the input sample list*/
	private int normalSample;
	
	private int filterCount0s;
	
	private double[] averageVar;
	//Part of a bizarre workaround for syntax problems in using the cluster class
	public SVReClusterer(){	}
	
	public SVReClusterer(ArrayList<SVEntry>  svs, ArrayList<SVEntry> ssnvNodes, ArrayList<SVEntry> ssnvs, int numS){
		toHandle = svs;
		combinedEntries = new ArrayList<SVEntry>();
		combinedEntries.addAll(svs);
		//combinedEntries.addAll(ssnvs);
		//combinedEntries.addAll(ssnvNodes);
		oldSSNVs = ssnvNodes;
		numSamples = numS;
		filterCount0s = 0;
		averageVar = new double[numS];
		
	}
	
	//Creates both the groups and the clusters, then returns the groups. I'll change the name later for readability purposes.
	public HashMap<String, SVGroup> createClusters(){
		// 1.Creates the groups based on the binary profiles of SVs and old SSNV nodes.
		createPreliminarySVGroups();
		
		// 2. get the SVs partitioned by group tag and create the appropriate SV group objects
		HashMap<String, ArrayList<SVEntry>> snvsByTag = tag2SVs;
		HashMap<String, SVGroup>  groups = new HashMap<String, SVGroup> ();
		for(String groupTag : snvsByTag.keySet()) {
			groups.put(groupTag, new SVGroup(groupTag, snvsByTag.get(groupTag), isRobustGroup(groupTag)));
		}
		if(groups.size() == 0) {
			//logger.warning("All SNV groups have been filtered out.");
			return groups;
		}
		//checkMissing(40,groups);
			
		
		for(int x = 0; x<numSamples; x++)
			averageVar[x]=0;
		
		// 3. cluster SNVs in each group
		ArrayList<String> toRemove = new ArrayList<String>();
		for(String groupId : groups.keySet()) {
			SVGroup group = groups.get(groupId);
			if(getHammingWeight(groupId)>0){//getHammingWeight(groupId)>0){
				//if(groupId.equals("00010100"))
					//System.out.println("test");
				Cluster[] clusters = em(group.getAlleleFreqBySample(), group.getNumSVs(), group.getNumSamples(), group.getSVs(), group.getSampleIds());
				//logger.fine("Clustering results for group: " + group.getTag());
				//for(Cluster c : clusters) {
					//System.out.println(c.getMembership().size());
				//	logger.fine(c.toString());
				//}
				group.setSubPopulations(clusters);
				//System.out.println(group.toString());
			}
			else{
				toRemove.add(groupId);
				for(SVEntry sv0: group.getSVs()){
					//System.out.println(sv0.getId()+"Hamming Weight too low");
					sv0.assignedNodeId=-2;
				}
				
				filterCount0s++;
			}
			
		}
		
		for(String s: toRemove)
			groups.remove(s);
		//checkMissing(40,groups);
		makeAverageVars(groups);
		
		assignAmbiguousSVs2ndHeuristic(groups);
		combinedEntries.addAll(oldSSNVs);
		createSVGroups();
		//checkMissing(40,groups);
		snvsByTag = tag2SVs;
		groups = new HashMap<String, SVGroup> ();
		
		for(String groupTag : snvsByTag.keySet()) {
			groups.put(groupTag, new SVGroup(groupTag, snvsByTag.get(groupTag), isRobustGroup(groupTag)));
		}
		//checkMissing(40,groups);
		if(groups.size() == 0) {
			//logger.warning("All SNV groups have been filtered out.");
			return groups;
		}
		
		toRemove = new ArrayList<String>();
		for(String groupId : groups.keySet()) {
			SVGroup group = groups.get(groupId);
			if(getHammingWeight(groupId)>0){
				
				Cluster[] clusters = em(group.getAlleleFreqBySample(), group.getNumSVs(), group.getNumSamples(), group.getSVs(), group.getSampleIds());
				//logger.fine("Clustering results for group: " + group.getTag());
				//for(Cluster c : clusters) {
				//	logger.fine(c.toString());
				//}
				
				group.setSubPopulations(clusters);
				//checkMissing(40,groups);
				//System.out.println(group.toString());
			}
			else{
				toRemove.add(groupId);
				for(SVEntry sv0: group.getSVs()){
					sv0.assignedNodeId=-2;
					//System.out.println(sv0.getId()+"Hamming Weight too low");
				}
				filterCount0s++;
			}
			
		}
		//checkMissing(40,groups);
		for(int x = 0; x<numSamples; x++)
			averageVar[x]=0;
		makeAverageVars(groups);
		return groups;
		
	}
	
	//Unused. Meant to make an arraylist of all clusters throughout all groups.
	public ArrayList<Cluster> extractClusters(){
		HashMap<String, SVGroup>  groups = createClusters();
		ArrayList<Cluster> clusters = new ArrayList<Cluster>();
		for(String groupId : groups.keySet()) {
			SVGroup svg = groups.get(groupId);
			for(Cluster c: svg.getSubPopulations()){
					clusters.add(c);
			}
		}
		
		return clusters;
	}
	
	//Generates clusters for a particular group based on EM algorithm on VAFs.
	public Cluster[] em(double[][] data, int numObs, int numFeatures, ArrayList<SVEntry> svEntries, int[] sampleIds) {
		Instances ds = convertMatrixToWeka(data, numObs, numFeatures);
		EM clusterer = new EM();
		try {
			//clusterer.setPreserveInstancesOrder(true);
			clusterer.buildClusterer(ds);
			ClusterEvaluation eval = new ClusterEvaluation();                                        
			eval.setClusterer(clusterer);                                  
			eval.evaluateClusterer(new Instances(ds));                               
			int numClusters = eval.getNumClusters();
			//System.out.println(eval.clusterResultsToString());//TEST PURPOSES
			
			Cluster[] clusters = new Cluster[numClusters];
			double[][] clusterCentroids = new double[numClusters][numFeatures];
			int[] clusterCount = new int[numClusters];
			
			double[] assignments = eval.getClusterAssignments();
			for(int i = 0; i < ds.numInstances(); i++) {
				Instance inst = ds.instance(i);
				int clusterId = (int) assignments[i];
				svEntries.get(i).setCluster((int)assignments[i]);///My attempt at identifying correct clustering. May not work as is. Unable to figure out if ordering is arbitrary after EM.
				for(int j = 0; j < numFeatures; j++) {
					clusterCentroids[clusterId][j] += inst.value(j);
				}
				clusterCount[clusterId]++;
			}
			
			for(int i = 0; i < numClusters; i++) {
				double[] mean = new double[numFeatures];
				for(int j = 0; j < numFeatures; j++) {
					mean[j] = clusterCentroids[i][j]/clusterCount[i];
				}
				clusters[i] = new Cluster(mean, i);
			}
			
			// cluster members & std dev
			double[][] clusterStdDev = new double[numClusters][numFeatures];
			for(int i = 0; i < ds.numInstances(); i++) {
				int clusterId = (int)assignments[i];
				clusters[clusterId].addMember(i);
				if(svEntries.get(i).description.equals("snv"))
					clusters[clusterId].setOriginal();
					
				for(int j = 0; j < numFeatures; j++) {
					clusterStdDev[clusterId][j] += Math.pow(ds.instance(i).value(j) - clusters[clusterId].getCentroid()[j], 2);
				}
			}
			
			for(int i = 0; i < numClusters; i++) {
				double[] dev = new double[numFeatures];
				for(int j = 0; j < numFeatures; j++) {
					dev[j] = Math.sqrt(clusterStdDev[i][j]/clusterCount[i]);
				}
				clusters[i].setStdDev(dev);
				clusters[i].makeStats(numSamples, sampleIds);
			}
			
			return clusters;
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
			return null;
		}
	}
	
	
	// ---- WEKA-related Utilities ----
		public Instances convertMatrixToWeka(double[][] data, int numObs, int numFeatures) {
			// convert the data to WEKA format
			FastVector atts = new FastVector();
			for(int i = 0; i < numFeatures; i++) {
				atts.addElement(new Attribute("Feature" + i, i));
			}
			
			
			
			Instances ds = new Instances("AAF Data", atts, numObs);
			for(int i = 0; i < numObs; i++) {
				ds.add(new Instance(numFeatures));
			}
			for(int i = 0; i < numFeatures; i++) {
				for(int j = 0; j < numObs; j++) {
					ds.instance(j).setValue(i, data[j][i]);
				}
			}
			return ds;
		}
		
		
		/**
		 * Cluster of observation points
		 * Each cluster has an associated centroid point and a list of members
		 */
		public class Cluster implements Serializable {
			
			private static final long serialVersionUID = 1L;

			/** Cluster id, unique per group */
			private int id;
			
			/** Cluster centroid */
			private double[] centroid;
			
			/** Cluster standard deviation */
			private double[] stdDev = null;
			
			/** List of observations assigned to this cluster */
			private ArrayList<Integer> members;
			
			private boolean robust = false;
			
			private NodeStats nodestat;
			
			boolean isOriginal;
			
			public Cluster(double[] clusterCentroid, int clusterId) {
				centroid = clusterCentroid;
				members = new ArrayList<Integer>();
				id = clusterId;
				isOriginal=false;
				
				//System.out.println(Arrays.toString(centroid));
				//System.out.println(this.toString());
			}
			
			public Cluster(double[] clusterCentroid, ArrayList<Integer> assignments, int clusterId) {
				centroid = clusterCentroid;
				members = assignments;
				id = clusterId;
				
			}
			
			public void setOriginal(){
				isOriginal=true;
			}
			
			public void makeStats(int numSamples, int[] indices){
				nodestat = new NodeStats(centroid, stdDev, numSamples, indices);
			}
			
			/**
			 * Compute the distance of a given observation to this cluster
			 * Currently the method computes the distance to the centroid
			 * @param x - observation point
			 * @param d - distance function to be used (e.g. Euclidean)
			 * @return distance of observation to cluster
			 */
			public double getDistanceToCluster(double[] x, DistanceMetric d) {
				return getDistance(x, centroid, d);
			}
			
			/**
			 * Add a new observation to the cluster
			 * @param obsId - Id of the observation (index in the data matrix)
			 */
			public void addMember(int obsId) {
				members.add(new Integer(obsId));
			}
			
			/** 
			 * Computes the mean of all the members of the cluster
			 */
			public void recomputeCentroidAndStdDev(double[][] data, int numObs, int numFeatures) {
				double[] newCentroid = new double[numFeatures];
				for(int i = 0; i < members.size(); i++) {
					for(int j = 0; j < numFeatures; j++) {
						newCentroid[j] += data[members.get(i)][j];
					}
				}
				for(int j = 0; j < numFeatures; j++) {
					newCentroid[j] = newCentroid[j]/members.size();
				}
				centroid = newCentroid;
				
				double[] clusterStdDev = new double[numFeatures];
				for(int i = 0; i < members.size(); i++) {
					for(int j = 0; j < numFeatures; j++) {
						clusterStdDev[j] += Math.pow(data[members.get(i)][j] - centroid[j], 2);
					}
				}
				for(int j = 0; j < numFeatures; j++) {
					clusterStdDev[j] =  Math.sqrt(clusterStdDev[j]/members.size());
				}
				stdDev = clusterStdDev;
			}
			
			/**
			 * Returns the cluster centroid (mean) per sample
			 */
			public double[] getCentroid() {
				return centroid;
			}
			
			/**
			 * Returns the standard deviation per sample
			 * @requires setStdDev() method to have been called (currently implemented for EM only),
			 * will return null otherwise
			 */
			public double[] getStdDev() {
				return stdDev;
			}
			
			public void setStdDev(double[] dev) {
				stdDev = dev;
			}
			
			public ArrayList<Integer> getMembership() {
				return members;
			}
			
			public int getId() {
				return id;
			}
			
			public boolean isRobust() {
				return robust;
			}
			
			public void setRobust() {
				robust = true;
			}
			
			public String toString() {
				String c = "";
				c += "Size: " + members.size() + "\n";
				DecimalFormat df = new DecimalFormat("#.##");
				c += "VAF Mean: [";
				for(int i = 0; i < centroid.length; i++) {
					c += " " + df.format(centroid[i]) + " ";
				}
				c += "] \n";
				c += "       Stdev:";
				if(stdDev != null) {
					c += " [";
					for(int i = 0; i < stdDev.length; i++) {
						c += " " + df.format(stdDev[i]) + " ";
					}
					c += "]";
				}
				return c;
			}
		}
		
		private double getDistance(double[] x, double[] y, DistanceMetric d) {
			switch(d) {
			case EUCLIDEAN:
				return getEuclideanDistance(x, y);
			case AVG_PER_SAMPLE:
				return getAvgSampleDistance(x, y);
			default:
				return 0;	
			}
		}
		
		/**
		 * Computes the Euclidean distance (2-norm distance) between 2 vectors
		 * @requires arrays x and y to have the same length
		 * @return Euclidean distance between x and y
		 */
		private double getEuclideanDistance(double[] x, double[] y) {
			double diffSum = 0;
			for(int i = 0; i < x.length; i++) {
				diffSum += Math.pow(Math.abs(x[i] - y[i]), 2);
			}
			return Math.sqrt(diffSum);
		}
		
		private double getAvgSampleDistance(double[] x, double[] y) {
			double diffSum = 0;
			for(int i = 0; i < x.length; i++) {
				diffSum += Math.abs(x[i] - y[i]);
			}
			return diffSum/x.length;
		}
		
		/** Distance measures */
		public enum DistanceMetric {
			EUCLIDEAN,
			AVG_PER_SAMPLE
		}
		
		public boolean isRobustGroup(String groupTag) {
			int numRobust = 0;
			for(SVEntry entry : tag2SVs.get(groupTag)) {
				if(entry.isRobust()) {
					numRobust++;
					if(numRobust >= Parameters.MIN_GROUP_PROFILE_SUPPORT) {
						return true;
					}
				}
			}
			return false;
		}
		
		/****** IO ******/
		
		//Creates groups out of SVs based on binary profiles. 
		public void createSVGroups(){
			//normalSample = normalSampleId;
			somaticSVs = new ArrayList<SVEntry>();
			tag2SVs = new HashMap<String, ArrayList<SVEntry>>();
			tag2Clusters = new HashMap<String, ArrayList<Cluster>>();
			ambiguousSVs = new ArrayList<SVEntry>();	
			loadSVs2();		
			
			


			//reportSVGroups();
			//if(Parameters.INPUT_FORMAT == Parameters.Format.SNV_WITH_PROFILE || clusterInputFile != null) return;
			
			// handle mutations from small groups as ambiguous
			ArrayList<String> smallGroups = new ArrayList<String>();
			for(String tag : tag2SVs.keySet()) {
				if(tag2SVs.get(tag).size() < Parameters.MIN_GROUP_PROFILE_SUPPORT) {
					for(SVEntry entry : tag2SVs.get(tag)) {
						ambiguousSVs.add(entry);
					}
					smallGroups.add(tag);
					
				}
			}
			for(String group : smallGroups) {
				//for(SVEntry sv: tag2SVs.get(group))
					//System.out.println(sv.getId()+" Group Too Small");
				tag2SVs.remove(group);
			}
			// assign ambiguous SNVs to existing groups or create new groups
			//assignAmbiguousSVs();
			// process resulting groups
			//filterGroups();
			//reportSVGroups();
		}
		
		
		public void createPreliminarySVGroups(){
			//normalSample = normalSampleId;
			somaticSVs = new ArrayList<SVEntry>();
			tag2SVs = new HashMap<String, ArrayList<SVEntry>>();
			tag2Clusters = new HashMap<String, ArrayList<Cluster>>();
			ambiguousSVs = new ArrayList<SVEntry>();	
			loadSVs();		
			
			


			//reportSVGroups();
			//if(Parameters.INPUT_FORMAT == Parameters.Format.SNV_WITH_PROFILE || clusterInputFile != null) return;
			
			// handle mutations from small groups as ambiguous
			ArrayList<String> smallGroups = new ArrayList<String>();
			for(String tag : tag2SVs.keySet()) {
				if(tag2SVs.get(tag).size() < Parameters.MIN_GROUP_PROFILE_SUPPORT) {
					for(SVEntry entry : tag2SVs.get(tag)) {
						ambiguousSVs.add(entry);
					}
					smallGroups.add(tag);
				}
			}
			for(String group : smallGroups) {
				//for(SVEntry sv: tag2SVs.get(group))
					//System.out.println(sv.getId()+" Group Too Small");
				tag2SVs.remove(group);
				
			}
			// assign ambiguous SNVs to existing groups or create new groups
			//assignAmbiguousSVs();
			// process resulting groups
			//filterGroups();
			//reportSVGroups();
		}
		
		//puts SVs into the data structures to be used by the remaining functions. Bridge between LICHeE code and lapsi code.
		public void loadSVs(){
			for(SVEntry entry: combinedEntries){
				somaticSVs.add(entry);	
				
				if(entry.isRobust()) {
					String tag = entry.getProfile();
					if (!tag2SVs.containsKey(tag)){
						tag2SVs.put(tag, new ArrayList<SVEntry>());	
					}
					tag2SVs.get(tag).add(entry); 
					
				} else {
					ambiguousSVs.add(entry);
				}
			}
		}
		
		//puts SVs into the data structures to be used by the remaining functions. Bridge between LICHeE code and lapsi code.
		public void loadSVs2(){
			for(SVEntry entry: combinedEntries){
				somaticSVs.add(entry);	
				
				//if(entry.isRobust()) {
				String tag = entry.getProfile();
				if (!tag2SVs.containsKey(tag)){
					tag2SVs.put(tag, new ArrayList<SVEntry>());	
				}
				tag2SVs.get(tag).add(entry); 
					
				//} else {
				//	ambiguousSVs.add(entry);
			}
		}
		
	
		/**
		 * Assign ambiguous SNVs to existing groups; 
		 * create new groups for SNVs with no suitable matches
		 */
		private void assignAmbiguousSVs() {
			boolean holdZero = false;//This may screw things up.
			if (ambiguousSVs.size() == 0) return;
			
			ArrayList<String> targetTags = new ArrayList<String>(tag2SVs.keySet());
			String all1s = "", all0s = "";
			for (int i = 0; i < numSamples; i++) { 
				all1s += "1";
				all0s += "0";
			}
			if(!tag2SVs.containsKey(all0s)) {
				targetTags.add(all0s);
			}
			else
				holdZero=true;//This may screw things up.
			targetTags.add(all1s);

			ArrayList<SVEntry> toRemove = new ArrayList<SVEntry>();
			for(SVEntry sv : ambiguousSVs) {	
				double bestDistToTarget = 0;
				String bestTarget = "";
					
				for(String target : targetTags) {				
					if(!canConvert(sv, target)) continue;
					// if the target is germline, move to germline regardless of distance
					if(target.equals(all1s)) {
						toRemove.add(sv);
						break;
					}				
					if(target.equals(all0s)) {
						double dist = distToZero(sv);
						if(dist > bestDistToTarget) {
							bestDistToTarget = dist;
							bestTarget = target;
						}
						continue;
					}
					// compare to each mutation in the target group
					for(SVEntry targetSNV : tag2SVs.get(target)) {
						// compute the distance to the target sv
						double dist = distToTargetVAF(sv, targetSNV);
						if(dist > bestDistToTarget) {
							bestDistToTarget = dist;
							bestTarget = target;
						}
					}
				}
				if(bestTarget.equals(all0s)) {
					sv.updateGroup(bestTarget);
					sv.assignedNodeId=-2;
					filterCount0s++;
					
					continue;
				}
				
				if(bestDistToTarget != 0 && bestDistToTarget >= getHammingWeight(bestTarget)*Parameters.MIN_VAF_TARGET_RATIO_PER_SAMPLE) {
					// found a valid match
					sv.updateGroup(bestTarget);
					toRemove.add(sv);
					tag2SVs.get(bestTarget).add(sv);
					
					continue;
				}
				
			
			}
			filterCount0s+=toRemove.size();
			ambiguousSVs.removeAll(toRemove);
			
			// remaining snvs had no suitable matches and potentially represent true branches
			// we minimize the number of additional nodes by merging the groups
			HashMap<String, ArrayList<SVEntry>> ambiguousGroups = mergeAmbiguousSVs(ambiguousSVs);
			
			for(String tag : ambiguousGroups.keySet()) {
				if(tag.equals(all0s)) continue;
				if(!tag2SVs.containsKey(tag)) {
					tag2SVs.put(tag, ambiguousGroups.get(tag));
				}
			}
			//if(!holdZero)//This may screw things up.
			tag2SVs.remove(all0s);
		}
		
		
		/**
		 * Assign ambiguous SNVs to existing groups; 
		 * create new groups for SNVs with no suitable matches
		 */
		private void assignAmbiguousSVs2ndHeuristic(HashMap<String, SVGroup> groups) {
			boolean holdZero = false;//This may screw things up.
			if (ambiguousSVs.size() == 0) return;
			
			ArrayList<String> targetTags = new ArrayList<String>(tag2SVs.keySet());
			String all1s = "", all0s = "";
			for (int i = 0; i < numSamples; i++) { 
				all1s += "1";
				all0s += "0";
			}
			if(!tag2SVs.containsKey(all0s)) {
				targetTags.add(all0s);
			}
			else
				holdZero=true;//This may screw things up.
			targetTags.add(all1s);

			ArrayList<SVEntry> toRemove = new ArrayList<SVEntry>();
			for(SVEntry sv : ambiguousSVs) {	
				double bestDistToTarget = 0;
				String bestTarget = "";
					
				for(String target : targetTags) {				
					if(!canConvert(sv, target)) continue;
					// if the target is germline, move to germline regardless of distance
					if(target.equals(all1s)) {
						System.out.println(sv.getId()+" Was germline");
						toRemove.add(sv);
						break;
					}				
					if(target.equals(all0s)) {
						double dist = probGenerate0(sv);
						
						
						if(dist > bestDistToTarget) {
							bestDistToTarget = dist;
							bestTarget = target;
						}
						continue;
					}
					// compare to each mutation in the target group
					for(Cluster cluster : groups.get(target).subPopulations) {
						// compute the distance to the target sv
						double dist = cluster.nodestat.probGenerateSV(sv);
						if(dist > bestDistToTarget) {
							bestDistToTarget = dist;
							bestTarget = target;
						}
					}
				}
				if(bestTarget.equals(all0s)) {
					sv.updateGroup(bestTarget);
					sv.assignedNodeId=-2;
					filterCount0s++;
					//System.out.println(sv.getId()+"Too Close to 0");
					continue;
				}
				
				if(bestDistToTarget != 0 && bestDistToTarget >= getHammingWeight(bestTarget)*Parameters.MIN_VAF_TARGET_RATIO_PER_SAMPLE) {
					// found a valid match
					sv.updateGroup(bestTarget);
					toRemove.add(sv);
					tag2SVs.get(bestTarget).add(sv);
					
					continue;
				}
				
				//System.out.println(sv.getId()+"Ambiguous + No Valid Match");
			}
			filterCount0s+=toRemove.size();
			ambiguousSVs.removeAll(toRemove);
			
			// remaining snvs had no suitable matches and potentially represent true branches
			// we minimize the number of additional nodes by merging the groups
			HashMap<String, ArrayList<SVEntry>> ambiguousGroups = mergeAmbiguousSVs(ambiguousSVs);
			
			for(String tag : ambiguousGroups.keySet()) {
				if(tag.equals(all0s)) continue;
				if(!tag2SVs.containsKey(tag)) {
					tag2SVs.put(tag, ambiguousGroups.get(tag));
				}
			}
			//if(!holdZero)//This may screw things up.
			tag2SVs.remove(all0s);
		}
		
		
		/**
		 * Group filtering based on:
		 * - group size
		 * - absence in all samples
		 * - minimum number of robust SNVs
		 */
		private void filterGroups() {
			// apply minimum size and robust size constraint
			ArrayList<String> filteredOut = new ArrayList<String>();
			for(String tag : tag2SVs.keySet()) {
				if(isAbsent(tag) || (tag2SVs.get(tag).size() < Parameters.MIN_SNVS_PER_GROUP)) {
					//System.out.println(isAbsent(tag));
					filteredOut.add(tag);
					continue;
				}
		
				if(Parameters.MIN_ROBUST_SNVS_PER_GROUP > 0) {
					int numRobust = 0;
					for(SVEntry entry : tag2SVs.get(tag)) {
						if(entry.isRobust()) {
							numRobust++;
							if(numRobust >= Parameters.MIN_ROBUST_SNVS_PER_GROUP) {
								break;
							}
						}
					}
					if(numRobust < Parameters.MIN_ROBUST_SNVS_PER_GROUP) {
						filteredOut.add(tag);
					}
				}
			}
			for(String group : filteredOut) {
				tag2SVs.remove(group);
			}
		}
		
		

		//Part of LICHeE's logging system. Not currently functional, but may be added back in later.
		public void reportSVGroups() {
			ArrayList<String> tags = new ArrayList<String>(tag2SVs.keySet());
			Collections.sort(tags);
			Collections.reverse(tags);

			
			for (String tag : tags){
				if (tag2SVs.get(tag).size() == 0) {
					continue;
				}
				
				double[] mean = new double[numSamples];	
				double[] robustMean = new double[numSamples];	
				int robustCount = 0;
				for (SVEntry entry : tag2SVs.get(tag)) {
					for (int i = 0; i < numSamples; i++){
						mean[i] += entry.getVAF(i);
						if(entry.isRobust()) {
							robustMean[i] += entry.getVAF(i);
						}
					}
					if(entry.isRobust()) {
						robustCount++;
					}
				}
				for (int i = 0; i < numSamples; i++) {
					mean[i] = mean[i]/(tag2SVs.get(tag).size());
					if(robustCount > 0) {
						robustMean[i] = robustMean[i]/robustCount;
					}
				}
				
				String m = tag + "\t" + tag2SVs.get(tag).size()+ "\t" + robustCount + "\t";
				DecimalFormat df = new DecimalFormat("#.##");
				for (int i = 0; i < numSamples; i++) {
					if(robustCount != 0) {
						m += df.format(mean[i]) + "(" + df.format(robustMean[i]) + ") \t";
					} else {
						m += df.format(mean[i]) + "(0) \t";
					}
				}

			}	
		}
		
		
		
		private double distToZero(SVEntry sv) {
			String tag = sv.getProfile();
			double dist = 0;
			for(int i = 0; i < numSamples; i++) {
				if(tag.charAt(i) == '1' || sv.evidenceOfPresence(i)) {
					dist += 0.0001/sv.getVAF(i);
				}
			}
			return dist;
		}
		
		private double distToTargetVAF(SVEntry snv, SVEntry targetSNV) {
			String tag = snv.getProfile();
			String targetTag = targetSNV.getProfile();
			double dist = 0;
			for(int i = 0; i < numSamples; i++) {
				if(targetTag.charAt(i) == '0' && tag.charAt(i) == '0' && !snv.evidenceOfPresence(i)) continue;
				double min = targetSNV.getVAF(i) < snv.getVAF(i) ? targetSNV.getVAF(i) : snv.getVAF(i);
				double max = targetSNV.getVAF(i) > snv.getVAF(i) ? targetSNV.getVAF(i) : snv.getVAF(i);
				if(min == 0) {
					min = 0.0001;
				}
				dist += min/max;
			}
			return dist;
		}
		
		private boolean canConvert(SVEntry sv, String target) {
			String tag = sv.getProfile();		
			for (int i = 0; i < tag.length(); i++) {
				if (tag.charAt(i) != target.charAt(i)) {
					if(tag.charAt(i) == '1') {
						return false;
					}
					if(!sv.evidenceOfPresence(i)) {
						return false;
					}
				}
			}
			return true;
		}
		
		private int getHammingWeight(String tag) {
			int w = 0;
			for(int i = 0; i < tag.length(); i++) {
				if(tag.charAt(i) == '1') {
					w++;
				}
			}
			return w;
		}
		
		
		private HashMap<String, ArrayList<SVEntry>> mergeAmbiguousSVs(ArrayList<SVEntry> svs) {
			HashMap<String, ArrayList<SVEntry>> groups = new HashMap<String, ArrayList<SVEntry>>();
			HashMap<String, ArrayList<SVEntry>> adj = new HashMap<String, ArrayList<SVEntry>>();
			
			// generate all possible target profiles 
			ArrayList<String> targets = generateAllPossibleTargets(svs);
			for(String t : targets) {
				//System.out.println("Target " + t);
				if(!adj.keySet().contains(t)) {
					adj.put(t, new ArrayList<SVEntry>());
				}
			}
			
			// min vertex cover
			for(SVEntry entry : svs) {
				for(String target : adj.keySet()) {
					if(canConvert(entry, target)) {
						adj.get(target).add(entry);
					}
				}
			}
			while(svs.size() > 0) {
				// find the largest set
				int maxSize = 0;
				String maxSet = "";
				for(String target : adj.keySet()) {
					if(adj.get(target).size() > maxSize) {
						maxSize = adj.get(target).size();
						maxSet = target;
					}
				}
				if(maxSize == 1) break;
				
				groups.put(maxSet, adj.get(maxSet));
				for(SVEntry entry : adj.get(maxSet)) {
					entry.updateGroup(maxSet);
					for(String target : adj.keySet()) {
						if(target.equals(maxSet)) continue;
						ArrayList<SVEntry> l = adj.get(target);
						if(l.contains(entry)) {
							l.remove(entry);
						}
					}
					svs.remove(entry);
				}
				adj.remove(maxSet);
			}
			
			// the remaining targets are supported by only 1 SNV
			// decide the groups of the remaining SNVs based on their proximity to the thresholds
			for(SVEntry snv : svs) {
				String tag = "";
				for(int i = 0; i < numSamples; i++) {
					if(snv.getProfile().charAt(i) == '0' && snv.evidenceOfPresence(i)) {
						double delta0 = snv.getVAF(i) - Parameters.MAX_VAF_ABSENT;
						double delta1 = Parameters.MIN_VAF_PRESENT - snv.getVAF(i);
						tag += (delta1 < delta0) ? '1' : snv.getProfile().charAt(i);
					} else {
						tag += snv.getProfile().charAt(i);
					}
				}
				
				if(!groups.containsKey(tag)) {
					groups.put(tag, new ArrayList<SVEntry>());
				}
				groups.get(tag).add(snv);
			}
			
			return groups;
		}
		
		
		private ArrayList<String> generateAllPossibleTargets(ArrayList<SVEntry> svs) {
			ArrayList<String> targets = new ArrayList<String>();
			for(SVEntry sv: svs) {
				extendTarget("", sv, targets);
			}
			return targets;
		}
		
		private void extendTarget(String partialProfile, SVEntry sv, ArrayList<String> targets) {
			int sample = partialProfile.length();
			if(sample == numSamples) {
				targets.add(partialProfile);
				return;
			}
			extendTarget(partialProfile + ""+ sv.getProfile().charAt(sample), sv, targets);
			if(sv.getProfile().charAt(sample) == '0' && sv.evidenceOfPresence(sample)) {
				extendTarget(partialProfile + '1', sv, targets);
			}
		}
		
		private boolean isAbsent(String tag) {
			for(int i = 0; i < numSamples; i++) {
				if(tag.charAt(i) == '1') {
					return false;
				}
			}
			return true;
		}
		
		public void makeAverageVars(HashMap<String, SVGroup> groups){
			int[] averageVarCounts = new int[numSamples];
			for(String s: groups.keySet()){
				SVGroup g = groups.get(s);
				
				for(int x = 0; x<numSamples; x++){
					if(g.containsSample(x))
						for(Cluster c: g.subPopulations){
							averageVar[x]+=c.stdDev[g.getSampleIndex(x)];
							averageVarCounts[x]++;
						}
				}
			}
			for(int x = 0; x<numSamples; x++)
				averageVar[x]=averageVar[x]/averageVarCounts[x];
		}
		
		public double probGenerate0(SVEntry sv){
			double prob = 1.0;
			for(int x = 0; x<numSamples; x++)
				prob = prob*pdf(sv.approxVAF[x], 0, Math.sqrt(averageVar[x]));
			return prob;
		}
		
		public static double pdf(double x) {
	        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
	    }

	    // return pdf(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
	    public static double pdf(double x, double mu, double sigma) {
	        return pdf((x - mu) / sigma) / sigma;
	    }
	    
	    public void checkMissing(int count,	HashMap<String, SVGroup> groups){
			boolean[] notMissing = new boolean[count];
			for(String s: groups.keySet()){
				if(groups.get(s).subPopulations!=null)
					for(Cluster c: groups.get(s).subPopulations)
						for(Integer i: c.getMembership()){
							SVEntry sv = groups.get(s).getSVs().get(i);
							if(sv.getId()!=null)
								notMissing[Integer.parseInt(sv.getId())-1]=true;
						}
			}
			for(int x = 0; x<notMissing.length; x++){
				if(!notMissing[x])
					System.out.println(x-1+" missing!");
			}
			System.out.println();
		}
}
