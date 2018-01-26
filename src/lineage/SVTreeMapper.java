package lineage;


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Set;



import lineage.SVReClusterer.Cluster;
import cc.mallet.types.Dirichlet;

/** 
 * Handles the assignment of SVs to tree nodes
 */
public class SVTreeMapper {
	private int numSamples;
	protected LineageTree tree;
	protected ArrayList<SVEntry> svs;
	//protected ArrayList<SVEntry> destructibleSVS;//Used for efficiency purposes. Is a copy of svs that can be destroyed without modifying svs;
	protected ArrayList<SVEntry> assignedSVs;
	protected ArrayList<SVEntry> noiseSVs;
	protected ArrayList<SVEntry> unassignedSVs;
	protected ArrayList<TreeNode> potentialNewNodes;
	protected ArrayList<SVEntry> oldSSNVs;
	protected SampleStats[] ssnvStats;//experimental
	protected SampleStats[] svStats;//experimental
	final static double SIMILARITY_CUTOFF = 0.5;
	final static double NOISE_LIKELIHOOD_CUTOFF = 0.5;
	protected int filterCount0s;
	protected int filterCountsSim;
	protected int filterCountNodes;
	HashMap <String, SVGroup> greturn;
	
	public SVTreeMapper(ArrayList<SVEntry> inputSVs, LineageTree inputTree, int numInputSamples, ArrayList<SVEntry> ssnvs)throws IOException {
		numSamples = numInputSamples;
		svs = inputSVs;
		//destructibleSVS = new ArrayList<SVEntry>(svs);
		tree = inputTree;
		assignedSVs = new ArrayList<SVEntry>();
		unassignedSVs = new ArrayList<SVEntry>();
		noiseSVs = new ArrayList<SVEntry>();
		potentialNewNodes = new ArrayList<TreeNode>();
		oldSSNVs = ssnvs;//experimental
		ssnvStats = new SampleStats[numInputSamples];//experimental
		svStats = new SampleStats[numInputSamples];//experimental

		for(int x = 0; x<ssnvStats.length; x++){
			ssnvStats[x] = new SampleStats(oldSSNVs, x);//experimental
			svStats[x] = new SampleStats(inputSVs, x);//experimental
		}
		FileWriter fwV = new FileWriter("testdata/translatedEstimatedVAFs.txt");
		for(int x = 0; x<inputSVs.size(); x++){
			SVEntry sv = inputSVs.get(x);
			for(int s = 0; s<svStats.length; s++){
				if(sv.approxVAF[s]!=0){
					sv.approxVAF[s]=ssnvStats[s].convertVAF(sv.approxVAF[s], svStats[s]);
					if(sv.approxVAF[s]<0)
						sv.approxVAF[s]=0;
					//else if(approxVAF[x]>.5)
						//approxVAF[x]=.5;
				}
			}
			fwV.write(sv.id+"\t"+Arrays.toString(sv.approxVAF)+"\n");
			//System.out.println(sv);
		}
		fwV.close();
		filterCount0s = 0;
		filterCountNodes = 0;
		filterCountsSim = 0;
		// any filtering...
		
		
		///TESTVAFS!!!
		//makeFakeTestVAFS2(svs);
		
	}
	
	public void makeFakeTestVAFS(ArrayList<SVEntry> realSVs){
		int fakeCount = 1000;
		for(int x =0; x<40; x++){
			for(Integer i: tree.nodes.keySet()){
				TreeNode node = tree.nodes.get(i);
				double[] fakeVAFs = new double[node.bProfile.length()];
				for(int vafInd = 0; vafInd<node.bProfile.length(); vafInd++){
					if(node.bProfile.charAt(vafInd)=='1')
						fakeVAFs[vafInd]=Math.random()*(.6-Parameters.MIN_VAF_PRESENT)+Parameters.MIN_VAF_PRESENT;
					else
						fakeVAFs[vafInd]=Math.random()*Parameters.MAX_VAF_ABSENT;
						
				}
				SVEntry fakeSV = new SVEntry(numSamples, ""+fakeCount);
				fakeCount++;
				fakeSV.approxVAF=fakeVAFs;
				fakeSV.description = "fake";
				fakeSV.makeProfile(Parameters.MAX_VAF_ABSENT, Parameters.MIN_VAF_PRESENT);
				realSVs.add(fakeSV);
			}
		}
	}
	public void makeFakeTestVAFS2(ArrayList<SVEntry> realSVs){
		int fakeCount = 1000;
		for(int x =0; x<40; x++){
			for(Integer i: tree.nodes.keySet()){
				TreeNode node = tree.nodes.get(i);
				double[] fakeVAFs = new double[node.bProfile.length()];
				double avVAF = Parameters.MIN_VAF_PRESENT/2.0+Parameters.MAX_VAF_ABSENT/2.0;
				for(int vafInd = 0; vafInd<node.bProfile.length(); vafInd++){
					if(node.bProfile.charAt(vafInd)=='1')
						fakeVAFs[vafInd]=Math.random()*(.6-avVAF)+avVAF;
					else
						fakeVAFs[vafInd]=Math.random()*avVAF;
						
				}
				SVEntry fakeSV = new SVEntry(numSamples, ""+fakeCount);
				fakeCount++;
				fakeSV.approxVAF=fakeVAFs;
				fakeSV.description = "fake";
				fakeSV.makeProfile(Parameters.MAX_VAF_ABSENT, Parameters.MIN_VAF_PRESENT);
				realSVs.add(fakeSV);
			}
		}
	}
	
	
	public void assignSVs() {
		
		for(SVEntry sv : svs) {
			double bestScore = 0;
			// compute the distance to each tree node
			for(TreeNode node : tree.nodes.values()) {
				double similarity = getSimilarityScore(sv, node);
				if(similarity > bestScore) {
					bestScore = similarity;
					sv.assignedNodeId = node.getNodeId();
				}
			}
			// check if we can assign this node (i.e. its score is within the allowed similarity threshold)
			if(bestScore > SIMILARITY_CUTOFF) {
				assignedSVs.add(sv);
			} else if(getNoiseLikelihood(sv) > NOISE_LIKELIHOOD_CUTOFF) {
				noiseSVs.add(sv);				
			} else {
				unassignedSVs.add(sv);
			}
		}
	}
	
	public void checkMissing(int count,	HashMap<String, SVGroup> groups){
		boolean[] notMissing = new boolean[count];
		for(String s: groups.keySet()){
			for(SVEntry sv: groups.get(s).getSVs()){
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
	
	//Cluster-together-based node assignment algorithm. Can be substituted for the assignSVs function without additional modification.
	public void assignSVsViaClustering(){
		
		ArrayList<SVEntry> formerSSNVs = new ArrayList<SVEntry>();
		
		for(Integer i: tree.nodes.keySet()){
			TreeNode tNode = tree.nodes.get(i);
			//tNode.translate(ssnvStats,svStats);
			SVEntry exportedEntry = tNode.exportEntryNonCompressed();
			exportedEntry.description="snv";
			//exportedEntry.translate(ssnvStats, svStats);
			//exportedEntry.makeProfile(.2,.4);
			formerSSNVs.add(exportedEntry);
			tNode.compress();
			
		}
		
		SVReClusterer reClusterer = new SVReClusterer(svs, formerSSNVs, numSamples);
		
		HashMap<String, SVGroup> groups = reClusterer.createClusters();
		//checkMissing(40,groups);
		HashMap<GCKey, ArrayList<SVEntry>> comboClusters = new HashMap<GCKey, ArrayList<SVEntry>>();
		String all0s = "";
		for (int i = 0; i < numSamples; i++) { 
			//all1s += "1";
			all0s += "0";
		}
		for(SVEntry ssnv: formerSSNVs)
		{
			if(ssnv.bProfile.equals(all0s))
				continue;
			GCKey key = new GCKey(ssnv.bProfile, ssnv.clusterID);
			if(comboClusters.containsKey(key))
				comboClusters.get(key).add(ssnv);//Denotes that this cluster has more than one former SSNV node in it.
			else{
				comboClusters.put(key, new ArrayList<SVEntry>());
				comboClusters.get(key).add(ssnv);
			}
					
		}
		for(GCKey key: comboClusters.keySet())
		{
			String groupName = key.getGroupName();
			int clusterIndex = key.getClusterIndex();
			SVGroup group = groups.get(groupName);
			if(group==null)
				continue;
			//if(group.subPopulations.length==0)
				//System.out.println("?");
			Cluster cluster = group.subPopulations[clusterIndex];
			
			ArrayList<SVEntry> toDeal = comboClusters.get(key);
			
			if(toDeal.size()==1){
				assignAll(tree.nodes.get(toDeal.get(0).assignedNodeId), cluster, group);
			}
			else{
				assignDivy(toDeal, cluster, group);
			}
		}
		//checkMissing(40,groups);
		checkAssignment();
		buildNodes(comboClusters, groups);
		ArrayList<SVGroup> groupOutput = new ArrayList<SVGroup>();
		for(String i: groups.keySet())
			groupOutput.add(groups.get(i));
		greturn = groups;
//		HashMap<Integer, ArrayList<Integer>> orEdges = new HashMap<Integer, ArrayList<Integer>>();
		
		//makePhy(groupOutput,numSamples, tree.edges);
		
	}
	
	public void makePhy(ArrayList<SVGroup> groups, int numSamples, HashMap<Integer, ArrayList<Integer>> orEdges){
		Args args = new Args();
		ArrayList<String> sampleNames = new ArrayList<String>();
		for(int x = 0; x<numSamples; x++)
			sampleNames.add(""+x);
		
		
		PHYNetwork constrNetwork = new PHYNetwork(groups, numSamples, orEdges);
		//logger.fine(constrNetwork.toString());
		
		// 5. find all the lineage trees that pass the VAF constraints
		ArrayList<PHYTree> spanningTrees = constrNetwork.getLineageTrees();  
		//logger.info("Found " + spanningTrees.size() + " valid tree(s)");
		
		
		//testPrintTrees(spanningTrees);
		
		if(spanningTrees.size() == 0) {
		//	logger.info("Adjusting the network...");	
			// if no valid trees were found, fix the network 
			// remove group nodes that are not robust, complete edges
			int delta = 0;
			do {
				int numNodes = constrNetwork.numNodes;
				constrNetwork = constrNetwork.fixNetwork();
				spanningTrees = constrNetwork.getLineageTrees();  
				delta = numNodes - constrNetwork.numNodes; 
			} while((delta != 0) && (spanningTrees.size() <= 0));
			if(spanningTrees.size() <= 0) {
				Parameters.ALL_EDGES = true;
				constrNetwork = new PHYNetwork(groups, numSamples, orEdges);
				spanningTrees = constrNetwork.getLineageTrees();
			}	
			//logger.info("Found " + spanningTrees.size() + " valid trees after network adjustments");	
		}
		
		// 6. evaluate/rank the trees
		if(spanningTrees.size() > 0) {
			constrNetwork.evaluateLineageTrees();
			//logger.fine("Top tree\nError score: " + spanningTrees.get(0).getErrorScore());	
			//logger.fine(spanningTrees.get(0).toString());
		} 
		
		// 7. result visualization
		if(args.showNetwork) {
			constrNetwork.displayNetwork();
		}
		if(spanningTrees.size() > 0) {
			for(int i = 0; i < spanningTrees.size(); i++) {
				if(i >= args.numShow && i > 0) break;
				LineageDisplayConfig display = new LineageDisplayConfig(constrNetwork, spanningTrees.get(i), sampleNames, args.color);
				if(i < args.numShow) {
					Visualizer.showLineageTree(display);
				}
				if(i == 0 && args.outputDOTFileName != null) { // top tree only for simplicity
					String parentDir = new File(args.outputDOTFileName).getParent();
					writeTreeToDOTFile(display.toDOT(parentDir), sampleNames, args);
				}
			}
			// 8. persistent storage	
			if(args.numSave > 0) {
				writeTreesToTxtFile(constrNetwork, spanningTrees, sampleNames, args);
			}	
		} 
	}
	
	public void testPrintTrees(ArrayList<PHYTree> spanningTrees){
		for(PHYTree ph: spanningTrees){
			for(PHYNode tNode: ph.treeNodes){
				System.out.println("Node ID:"+tNode.getNodeId());
				System.out.println(tNode.getLabel());
				if(tNode.snvGroup!=null)
					for(SVEntry sv: tNode.getSVs(tNode.snvGroup.getSVs())){
						System.out.print(sv.id+" ");
					}
				System.out.print("\n");
			}
			for(PHYNode tFrom: ph.treeEdges.keySet()){
				
				System.out.println(tFrom.getNodeId());
				for(PHYNode tTo: ph.treeEdges.get(tFrom)){
					System.out.println("\t->"+tTo.getNodeId());
				}
			}
		}
	}
	
	public void assignAll(TreeNode t, Cluster c, SVGroup g){
		for(Integer index: c.getMembership()){
			SVEntry sv = g.getSVs().get(index);
			if(sv.assignedNodeId==-1){
				
				double similarity = getSimilarityScore(sv, t);
				if(similarity > SIMILARITY_CUTOFF) {
					assignedSVs.add(sv);
					t.receiveSV(sv);
					sv.assignedNodeId=t.getNodeId();
				} else if(getNoiseLikelihood(sv) > NOISE_LIKELIHOOD_CUTOFF) {
					noiseSVs.add(sv);
					filterCountsSim++;
					sv.assignedNodeId=-2;
				} //else
					//unassignedSVs.add(sv);
					
			}
		}
	}
	
	//Used when multiple SSNVs are found in the same cluster. Very similar to the orignal assignSVs function
	public void assignDivy(ArrayList<SVEntry> oldEntries, Cluster c, SVGroup g){
		for(Integer index: c.getMembership()){
			SVEntry sv = g.getSVs().get(index);
			
			double bestScore = 0;
			SVEntry bestEntry = null;
			for(SVEntry ssnv: oldEntries){
				double similarity = getSimilarityScore(sv, tree.nodes.get(ssnv.assignedNodeId));
				if(similarity>=bestScore){
					bestScore = similarity;
					bestEntry = ssnv;
				}
			}
			
			TreeNode t = tree.nodes.get(bestEntry.assignedNodeId);
			
			if(sv.assignedNodeId==-1){
				
				double similarity = getSimilarityScore(sv, t);
				if(similarity > SIMILARITY_CUTOFF) {
					assignedSVs.add(sv);
					t.receiveSV(sv);
					sv.assignedNodeId=t.getNodeId();
				} else if(getNoiseLikelihood(sv) > NOISE_LIKELIHOOD_CUTOFF) {
					noiseSVs.add(sv);
					filterCountsSim++;
					sv.assignedNodeId=-2;
				} //else
					//unassignedSVs.add(sv);
				
			}
		}
		
		
		
	}
	
	
	//Finds the unassigned nodes. There will probably be a better way to do this later.
	public void checkAssignment(){
		for(SVEntry s: svs){
			if(s.assignedNodeId==-1){
				unassignedSVs.add(s);
			}
		}
	}
	
	// remaining SVs had no suitable matches and potentially represent true branches
	// create new nodes for SVs that cannot be assigned
	public void handleUnassigned() {
		attachNewNodes(potentialNewNodes,tree);
	}
	
	
	//No longer used. Originally a simplified version of our node-addition algorithm.
	public void attachAuxiliaryNodes(HashMap<TreeNode, ArrayList<SVEntry>> auxNodes, LineageTree tree){
		
		///Does not use binary profiles to reduce checks at this time. Checks all nodes instead.
		PriorityQueue<WeightedEdge> potentialEdges = new PriorityQueue<WeightedEdge>();
		for(int nodeID: tree.nodes.keySet()){
			TreeNode parentNode = tree.nodes.get(nodeID);
			for(TreeNode child: auxNodes.keySet()){
				potentialEdges.add(new WeightedEdge(parentNode, child, getMinMaxScore2Trees(parentNode, child)));
			}
		}
		
		///NOTE! THIS KEEPS TRACK OF ADDED NODES BY DESTROYING THE auxNodes HashMap!
		while(!potentialEdges.isEmpty()){
			WeightedEdge e = potentialEdges.remove();
			if(auxNodes.containsKey(e.getChild())||!meetsConstraints(e.getParent(), e.getChild(), tree))
				continue;
			attachAuxNode(e, tree);
		}
		
	}
	
	//The function to call to add all the unassigned nodes;
	public void attachNewNodes(ArrayList<TreeNode> auxNodes, LineageTree tree){
		//ArrayList<TreeNode> toBeAdded = new ArrayList<TreeNode>(auxNodes.keySet());//Originally I used a priorityqueue, but I modify it so many times there wasn't any point.
		ArrayList<TreeNode> toBeAdded = auxNodes;
		
		ArrayList<WeightedEdgeSystem> newSystems = generateAllSystems(tree, toBeAdded);//Makes all possible WeightedEdgeSystems that meet constraints 1 and 2
		Collections.sort(newSystems);//Sorts them by net weight gain. weight function is the similarity score at this time.
		HashSet<Integer> addedNodes = new HashSet<Integer>();//Keeps track of which nodes have been used already without removing them from previous lists.
		
		//Runs until there are no more possible weightedEdgeSystems, which also implies there are no more possible nodes to add.
		while(!newSystems.isEmpty()){
			WeightedEdgeSystem newSystem = newSystems.remove(0);//Takes the highest weighted WeightedEdgeSystem out of the list.
			if(meetsConstraint3(newSystem.getParentNode(), newSystem.getNewNode(), tree)){
				if(applySystem(newSystem, tree)){//Makes the changes listed in the WeightedEdgeSystem
					addedNodes.add(newSystem.getNewNode().getNodeId());//Adds the newly added node to the tree's hashmap.
					addNewSystems(newSystem.getNewNode(), tree, toBeAdded, newSystems, addedNodes);//Makes new WeightedEdgeSystems based around the new node
					deleteOldSystems(newSystem, newSystems, addedNodes);//Removes old WeightedEdgeSystems that require an edge that was broken during the most recent change.
					Collections.sort(newSystems);//Resorts the weightedEdgeSystems.
				}
			}
			
		}
		
		
	}
	
	//Creates all possible weightedEdgeSystems involving unassigned nodes and the former lineage tree.
	public ArrayList<WeightedEdgeSystem> generateAllSystems(LineageTree tree, ArrayList<TreeNode> nodes){
		ArrayList<WeightedEdgeSystem> newSystems = new ArrayList<WeightedEdgeSystem>();
		
		Set<Integer> treeNodeSt = tree.nodes.keySet();
		double[] rootVAFs = new double[tree.numSamples];
		int[] rootInds = new int[tree.numSamples];
		double[] rootVars = new double[tree.numSamples];
		for(int x = 0; x<tree.numSamples; x++){
			rootVAFs[x]=.5;
			rootInds[x]=x;
			rootVars[x]=Math.pow(svStats[x].std,2);
		}
		NodeStats rootStats = new NodeStats(rootVAFs, rootVars, numSamples);
		TreeNode rootNode = new TreeNode(0, rootVAFs, rootInds, numSamples);
		rootNode.nodeStats = rootStats;
		tree.nodes.put(0, rootNode);
		for(Integer i: tree.nodes.keySet()){
			TreeNode oldNode = tree.nodes.get(i);
			
			ArrayList<TreeNode> candidates = new ArrayList<TreeNode>();
			for(TreeNode newNode: nodes)
				if(meetsConstraints1And2(oldNode, newNode, tree))
					candidates.add(newNode);
			if(candidates.isEmpty())
				continue;
			
			
			 
			
			for(TreeNode newNode: candidates){
				ArrayList<TreeNode> children = new ArrayList<TreeNode>();
				if(tree.edges.containsKey(i)){
					
					for(Integer j: tree.edges.get(i)){
						TreeNode childNode = tree.nodes.get(j);
						if(meetsConstraints1And2(newNode, childNode, tree))
							children.add(tree.nodes.get(j));
					}
					ArrayList<ArrayList<TreeNode>> childLists = powerList(children);
					
					ArrayList<WeightedEdge> toMake = new ArrayList<WeightedEdge>();
					toMake.add(new WeightedEdge(oldNode, newNode, oldNode.nodeStats.probGenerate(newNode.nodeStats)));
					ArrayList<WeightedEdge> toBreak = new ArrayList<WeightedEdge>();
					for(ArrayList<TreeNode> toChange: childLists){
						for(TreeNode t: toChange){
							toMake.add(new WeightedEdge(newNode, t, t.nodeStats.probGenerate(newNode.nodeStats)));
							toBreak.add(new WeightedEdge(oldNode, t, t.nodeStats.probGenerate(newNode.nodeStats)));
						}
					}
					
					
					WeightedEdgeSystem newSystem = new WeightedEdgeSystem(newNode, oldNode, toBreak, toMake);
					newSystems.add(newSystem);
				}
				else{
					ArrayList<WeightedEdge> toMake = new ArrayList<WeightedEdge>();
					toMake.add(new WeightedEdge(oldNode, newNode, newNode.nodeStats.probGenerate(oldNode.nodeStats)));
					WeightedEdgeSystem newSystem = new WeightedEdgeSystem(newNode, oldNode, new ArrayList<WeightedEdge>(), toMake);
					newSystems.add(newSystem);
				}
			}
		}
		
		
		
		return newSystems;
	}
	
	
	//Adds new WeightedEdgeSystems to the list of possible tree modifications based around the new node that was just added to the tree in the last modification.
	public void addNewSystems(TreeNode addedNode, LineageTree tree, ArrayList<TreeNode> nodes, ArrayList<WeightedEdgeSystem> systems, HashSet<Integer> done){
		
		TreeNode oldNode = addedNode;
		
		ArrayList<TreeNode> candidates = new ArrayList<TreeNode>();
		for(TreeNode newNode: nodes)
			if(!done.contains(newNode) && meetsConstraints1And2(oldNode, newNode, tree))
				candidates.add(newNode);
		if(candidates.isEmpty())
			return;
		
		
		
		
		for(TreeNode newNode: candidates){
			ArrayList<TreeNode> children = new ArrayList<TreeNode>();
			if(tree.edges.containsKey(addedNode.getNodeId()))
				for(Integer j: tree.edges.get(addedNode.getNodeId())){
					TreeNode childNode = tree.nodes.get(j);
					if(meetsConstraints1And2(newNode, childNode, tree))
						children.add(tree.nodes.get(j));
				}
			ArrayList<ArrayList<TreeNode>> childLists = powerList(children);
			
			ArrayList<WeightedEdge> toMake = new ArrayList<WeightedEdge>();
			toMake.add(new WeightedEdge(newNode, oldNode, newNode.nodeStats.probGenerate(oldNode.nodeStats)));
			ArrayList<WeightedEdge> toBreak = new ArrayList<WeightedEdge>();
			for(ArrayList<TreeNode> toChange: childLists){
				for(TreeNode t: toChange){
					toMake.add(new WeightedEdge(newNode, t, getMinMaxScore2Trees(t, newNode)));
					toBreak.add(new WeightedEdge(oldNode, t, getMinMaxScore2Trees(t, newNode)));
				}
			}
			
			
			WeightedEdgeSystem newSystem = new WeightedEdgeSystem(newNode, oldNode, toBreak, toMake);
			systems.add(newSystem);
			
		}
		
	}
	
	//Searches the possible modifications to the tree for WeightedEdgeSystems that are rendered invalid by the last change to the tree and removes them.
	public void deleteOldSystems(WeightedEdgeSystem modifications, ArrayList<WeightedEdgeSystem> systems, HashSet<Integer> done){
		ArrayList<WeightedEdge> wasBroken = modifications.toBreak;
		for(int x = 0; x<systems.size(); x++){
			if(done.contains(systems.get(x).getNewNode().getNodeId())){
				systems.remove(x);
				x--;
				continue;
			}
			
			
			if(systems.get(x).getParentNode().getNodeId()==modifications.getParentNode().getNodeId())
				for(WeightedEdge e: wasBroken)///May want to replace toMake and toBreak with HashSets.
					if(systems.get(x).toMake.contains(e))
					{
						systems.remove(x);
						x--;
						break;
					}
		}
	}
	
	//Makes changes to the lineage tree based on the WeightedEdgeSystem's specifications
	public boolean applySystem(WeightedEdgeSystem system, LineageTree tree){
		if(tree.nodes.containsKey(system.getNewNode()))
			return false;
		for(WeightedEdge edge: system.toMake){
			if(!meetsConstraints(edge.getParent(), edge.getChild(), tree)){
				return false;
			}
		}
		for(WeightedEdge edge: system.toBreak){
			if(tree.edges.get(edge.getParentId())==null || !tree.edges.get(edge.getParentId()).contains(edge.getChildId())){
				return false;
			}
		}
		for(WeightedEdge edge: system.toMake)
			addEdge(edge, tree);
		for(WeightedEdge edge: system.toBreak)
			breakEdge(edge, tree);
		tree.nodes.put(system.getNewNode().getNodeId(), system.getNewNode());
		
		return true;
	}
	
	
	//Dead end function. Does nothing at this time.
	public void attachAllNodeTypesExperimental(HashMap<TreeNode, ArrayList<SVEntry>> auxNodes, LineageTree tree){
		ArrayList<WeightedEdgeSystem> toCheck = new ArrayList<WeightedEdgeSystem>();
	}
	
	//Adds both an edge and a new node to the tree
	public void attachAuxNode(WeightedEdge e, LineageTree tree){
		tree.nodes.put(e.getChildId(), e.getChild());
		tree.edges.get(e.getParentId()).add(e.getChildId());
	}
	
	//Adds the corresponding edge to the tree. DOES NOT ADD EITHER NODE IN THE EDGE TO THE TREE.
	public void addEdge(WeightedEdge e, LineageTree tree){
		//if(tree.edges.containsKey(e.getParentId()))
		//	tree.edges.get(e.getParentId()).add(e.getChildId());
		//else
		//	tree.addEdge
		tree.addEdge(e.getParentId(), e.getChildId());
	}
	
	//Places a new node into the lineage tree.
	public void addNode(TreeNode t, LineageTree tree){
		tree.nodes.put(t.getNodeId(), t);
	}
	
	//Takes the corresponding edge e, and removes it from the tree.
	public void breakEdge(WeightedEdge e, LineageTree tree){
		tree.edges.get(e.getParentId()).remove((Integer)e.getChildId());///I might want to make edges a hashmap of integers to hashsets instead of arraylists.
	}
	
	public double getSimilarityScore(SVEntry sv, TreeNode treeNode) {
		//return getMinMaxScore(sv, treeNode);
		return treeNode.nodeStats.probGenerateSV(sv);
	}
	
	public double getNoiseLikelihood(SVEntry sv) {
		return 0;
	}
	
	private double getMinMaxScore(SVEntry sv, TreeNode treeNode) {
		double score = 0;
		for(int i = 0; i < numSamples; i++) {
			int tIndex = treeNode.getSampleIndex(i);
			if(tIndex==-1)
				continue;
			double min = treeNode.getVAF(tIndex) < sv.getVAF(i) ? treeNode.getVAF(tIndex) : sv.getVAF(i);
			double max = treeNode.getVAF(tIndex) > sv.getVAF(i) ? treeNode.getVAF(tIndex) : sv.getVAF(i);
			if(min == 0) {
				min = 0.0001;
			}
			if(max!=0)
				score += min/max;
			else
				score+=1;
		}
		return score;
	}
	
	//Just a duplicate function of getMinMaxScore only with Tree inputs instead for convenience
	private double getMinMaxScore2Trees(TreeNode sv, TreeNode treeNode) {
		double score = 0;
		for(int i = 0; i < numSamples; i++) {
			int sIndex = sv.getSampleIndex(i);
			int tIndex = treeNode.getSampleIndex(i);
			if(sIndex==-1 || tIndex==-1)
				continue;
			
			double min = treeNode.getVAF(tIndex) < sv.getVAF(sIndex) ? treeNode.getVAF(tIndex) : sv.getVAF(sIndex);
			double max = treeNode.getVAF(tIndex) > sv.getVAF(sIndex) ? treeNode.getVAF(tIndex) : sv.getVAF(sIndex);
			if(min == 0) {
				min = 0.0001;
			}
			if(max!=0)
				score += min/max;
			else
				score+=1;
		}
		return score;
	}
	
	/**
	 * Finds the minimum number of nodes that can incorporate the unassigned input SVs
	 * applying the greedy set cover algorithm
	 */
	protected HashMap<TreeNode, ArrayList<SVEntry>> minSetCover() {
		HashMap<TreeNode, ArrayList<SVEntry>> assignments = new HashMap<TreeNode, ArrayList<SVEntry>>();
		HashMap<TreeNode, ArrayList<SVEntry>> graph = new HashMap<TreeNode, ArrayList<SVEntry>>();
		
		ArrayList<TreeNode> targets = generateAllPossibleTargets(unassignedSVs);
		for(TreeNode t : targets) {
			if(!graph.keySet().contains(t)) {//Wouldn't "containsKey" be more efficient here? (Daniel Comment)
				graph.put(t, new ArrayList<SVEntry>());
			}
		}	
		for(SVEntry entry : unassignedSVs) {
			for(TreeNode target : graph.keySet()) {
				if(canAssign(entry, target)) {
					graph.get(target).add(entry);
				}
			}
		}
		while(unassignedSVs.size() > 0) {
			// find the largest set
			int maxSize = 0;
			TreeNode maxSet = null;
			for(TreeNode target : graph.keySet()) {
				if(graph.get(target).size() > maxSize) {
					maxSize = graph.get(target).size();
					maxSet = target;
				}
			}
			assignments.put(maxSet, graph.get(maxSet));
			for(SVEntry entry : graph.get(maxSet)) {
				entry.assignedNodeId = maxSet.getNodeId();
				for(TreeNode target : graph.keySet()) {
					if(target.equals(maxSet)) continue;
					ArrayList<SVEntry> l = graph.get(target);
					if(l.contains(entry)) {
						l.remove(entry);
					}
				}
				unassignedSVs.remove(entry);
			}
			graph.remove(maxSet);
		}
		return assignments;
	}
	
	
	//This is Camir's. Not sure what it's for.
	private ArrayList<TreeNode> generateAllPossibleTargets(ArrayList<SVEntry> svs) {
		ArrayList<TreeNode> targets = new ArrayList<TreeNode>();
		return targets;
	}
	
	
	//This is Camir's. Not sure what it's for.
	private boolean canAssign(SVEntry sv, TreeNode target) {
		return true;
	}	
	
	private boolean meetsConstraints(TreeNode parent, TreeNode child, LineageTree tree){
		/**boolean parentVAFGreater = true;
		boolean bothZero = true;
		boolean parentSumGreater = true;
		
		
		
		for(int vInd = 0; vInd< parent.getVAFLength(); vInd++){
			if(parent.getVAF(vInd)<child.getVAF(vInd)-noiseErrorMargin())
				parentVAFGreater = false;
			if(parent.getVAF(vInd)==0 && child.getVAF(vInd)!=0)
				bothZero = false;
			
			double parentSum = parent.getVAF(vInd) + noiseErrorMargin();
			double childSum = child.getVAF(vInd);
			
			for(int tID: tree.edges.get(parent.getNodeId())){
				TreeNode altChild = tree.nodes.get(tID);
				childSum+=altChild.getVAF(vInd);
			}
			if(parentSum<childSum)
				parentSumGreater = false;
		}
		
		
		return parentVAFGreater && bothZero && parentSumGreater;**/
		return meetsConstraints1And2(parent, child, tree) && meetsConstraint3(parent, child, tree);
	}
	
	
	//Checks that all the parent's VAFs are greater than its child's, and that if a Parent's VAFs are 0, so are the childs.
	private boolean meetsConstraints1And2_oldVersion(TreeNode parent, TreeNode child, LineageTree tree){
		boolean parentVAFGreater = true;
		boolean bothZero = true;
		//if(debugconstraint1test1(parent, child)){
		//	System.out.println("Trying");
		//}
		for(int vInd = 0; vInd< parent.getFullSampleLength(); vInd++){//MIGHT NEED TO CHANGE getVAFLength() to numSamples!
			
			int pInd = parent.getSampleIndex(vInd);
			int cInd = child.getSampleIndex(vInd);
			if(cInd==-1){
				if(pInd!=-1 && parent.getVAF(pInd)!=0)
					return false;
				else
					continue;
			}
			else if(child.getVAF(cInd)==0){//I may fix this later so that it is unnecessary
				if(pInd==-1 || parent.getVAF(pInd)==0)
					continue;
			}
			else if(pInd==-1||parent.getVAF(pInd)==0)
				continue;
			
			double pVAF = parent.getVAF(parent.getSampleIndex(vInd));
			double cVAF = child.getVAF(child.getSampleIndex(vInd));
			if(pVAF<cVAF-noiseErrorMargin())
				return false;
			if(pVAF==0 && cVAF!=0)
				return false;
		}
		
		
		return true;
	}
	private boolean meetsConstraints1And2(TreeNode parent, TreeNode child, LineageTree tree){
		//if(parent.getNodeId()==0)
			//System.out.println("CheckThis");
		//if(child.getNodeId()==1)
			//System.out.println("WHAT.");
		
		for(int vInd = 0; vInd< parent.getFullSampleLength(); vInd++){//MIGHT NEED TO CHANGE getVAFLength() to numSamples!
			
			int pInd = parent.getSampleIndex(vInd);
			int cInd = child.getSampleIndex(vInd);
			if(pInd==-1 || parent.getVAF(parent.getSampleIndex(vInd))==0){
				if(cInd==-1)
					continue;
				else if(child.getVAF(child.getSampleIndex(vInd))==0)
					continue;
				else
					return false;
			}
			if(cInd==-1)
				continue;
			double pVAF = parent.getVAF(parent.getSampleIndex(vInd));
			double cVAF = child.getVAF(child.getSampleIndex(vInd));
			if(pVAF<cVAF-noiseErrorMargin())
				return false;
			if(pVAF==0 && cVAF!=0)
				return false;
		}
		
		
		return true;
	}
	
	public boolean debugconstraint1test1(TreeNode Parent, TreeNode Child){
		int count1 = 0;
		int count2 = 0;
		for(int vf: Parent.vafIndices){
			if(vf!=0)
				count1++;
			for(int vf2: Child.vafIndices)
				if(vf!=0 && vf==vf2)
					count2++;
		}
		return count2==count1;
	}
	
	
	//This checks if the parent node has VAFs greater than the sum of its child VAFs.
	//Need to be separate rom the other constraints because it changes during the tree modification process, and checking it simultaneously with them would be useless.
	private boolean meetsConstraint3(TreeNode parent, TreeNode child, LineageTree tree){
		boolean parentSumGreater = true;
		
		//if(parent.getNodeId()==0)
		//	System.out.println("RootCheck");
			
		for(int vInd = 0; vInd< parent.getFullSampleLength(); vInd++){//MIGHT NEED TO CHANGE getVAFLength() to numSamples!
			int pInd = parent.getSampleIndex(vInd);
			
			int cInd = child.getSampleIndex(vInd);
			//if(pInd==-1 && (cInd!=-1||child.getVAF(cInd)!=0))//Hopefully redundant. Will experiment with removing it later.
			//	return false;
			
			if(cInd==-1)
				continue;
			
			double parentSum =  noiseErrorMargin();
			double childSum = 0;
			
			if (cInd!=-1)
				childSum += child.getVAF(cInd);
			if(pInd!=-1)
				parentSum += parent.getVAF(pInd);
			 
			if(tree.edges.containsKey(parent.getNodeId())){
				ArrayList<Integer> testChildren = tree.edges.get(parent.getNodeId());
				//if(testChildren.size()>4)
					//System.out.println("TESTSIZE!");
				for(int tID: tree.edges.get(parent.getNodeId())){
					
					TreeNode altChild = tree.nodes.get(tID);
					
					cInd = altChild.getSampleIndex(vInd);
					if (cInd!=-1)
						childSum += altChild.getVAF(cInd);
					
					//childSum+=altChild.getVAF(vInd);
				}
			}
			if(parentSum<childSum)
				return false;//parentSumGreater = false;
		}
		
		return true;//parentSumGreater;
	}
	
	//May need to add some filtering for low significance nodes later.
	public void buildNodes(HashMap<GCKey, ArrayList<SVEntry>> usedClusters, HashMap<String, SVGroup> groups){
		
		int count = tree.nodes.size()+1;
		for(String bProfile: groups.keySet()){
			SVGroup group = groups.get(bProfile);
			if(group.subPopulations!=null)
				for(int clusterIndex = 0; clusterIndex<group.subPopulations.length; clusterIndex++){
					Cluster cluster= group.getSubPopulations()[clusterIndex];
					if(cluster.getMembership().size()==0)
						continue;//Filters empty clusters, which apparently can happen from weka EM.
					if(!usedClusters.containsKey(new GCKey(bProfile, clusterIndex))){
						TreeNode newNode = new TreeNode(count, cluster.getCentroid(), group.getSampleIds(), group.getNumSamples());
						potentialNewNodes.add(newNode);
						for(Integer i: cluster.getMembership()){
							newNode.receiveSV(group.getSVs().get(i));//May not be necessary in the long run, but stores references to SVs in the node.
							group.getSVs().get(i).setNodeId(count);
						}
						newNode.decompress();
						newNode.makeNodeStats();
						newNode.bProfile=bProfile;
						count++;
						//if(count==59)
							//System.out.println("test");
					}
				}
		}
		//for(Integer i: tree.nodes.keySet())
			//System.out.println(tree.nodes.get(i));
		//System.out.println();
		//for(TreeNode treeN: potentialNewNodes)
			//System.out.println(treeN);
			
	}
	
	
	
	
	private double noiseErrorMargin(){
		return Parameters.VAF_ERROR_MARGIN;//Need to ask about how to get this. (Daniel Comment)
	}
	
	
	
	
	//For the purpose of hashing on both the ClusterID and the GroupID. Used for retrieval of clusters which share multiple SSNV nodes from the former tree.
	public class GCKey {

	    private final String bProfile;
	    private final int cIndex;
	    private String fullKey;
	    
	    
	    public GCKey(String x, int y) {
	        bProfile = x;
	        cIndex = y;
	        fullKey = bProfile+cIndex;
	    }

	    @Override
	    public boolean equals(Object o) {
	        if (this == o) return true;
	        if (!(o instanceof GCKey)) return false;
	        GCKey key = (GCKey) o;
	        return fullKey.equals(key.fullKey);
	    }

	    @Override
	    public int hashCode() {
	        String result = fullKey;
	        return result.hashCode();
	    }
	    
	    public String getGroupName(){
	    	return bProfile;
	    }
	    
	    public int getClusterIndex(){
	    	return cIndex;
	    }

	}
	
	//Helper function for powerList
	public void powerListRecursion(ArrayList<TreeNode> input, ArrayList<TreeNode> extend, ArrayList<ArrayList<TreeNode>> lists, int integerIndex){
		ArrayList<TreeNode> extend2 = new ArrayList<TreeNode>(extend);
		extend.add(input.get(integerIndex));
		integerIndex++;
		if(integerIndex==input.size()){
			lists.add(extend);
			lists.add(extend2);
		}
		else{
			powerListRecursion(input, extend, lists, integerIndex);
			powerListRecursion(input, extend2, lists, integerIndex);
		}
	}
	
	//Used for finding all subsets of possible children to shift around for finding hypothetical weighted edge systems to modify the existing lineage tree with.
	public ArrayList<ArrayList<TreeNode>> powerList(ArrayList<TreeNode> input){
		ArrayList<ArrayList<TreeNode>> output = new ArrayList<ArrayList<TreeNode>>();
		if(input.size()!=0)
			powerListRecursion(input, new ArrayList<TreeNode>(), output, 0);
		return output;
	}
	
	
	public HashMap<Integer, TreeNode> getNodes(){
		return tree.nodes;
	}
	
	public ArrayList<TreeNode> getUnaddedNodes(){
		return potentialNewNodes;
	}
	
	//Test Technique for attaching nodes.
	public void newNodeTechnique1(){
		attachNewNodes(potentialNewNodes,tree);
	}
	
	protected static class Args {
		// --- 'build' command ---
		String inputFileName;
		String outputFileName;
		String outputDOTFileName;	
		String clustersFileName;
		int normalSampleId = 0;
		String cnvFileName;
		String annFileName;
		String cosmicFileName;
		String tcgaFileName;
		
		// --- 'show' command ---
		String showFileNamePrefix;
		int numShow = 0;
		int numSave = 10;
		
		// flags
		boolean showNetwork = false;
		boolean verbose = false;
		boolean color = false;
		protected Args(){
			showFileNamePrefix="showFile";
			outputFileName="newTreesFromNonGreedy/Trees.txt";
			outputDOTFileName = "newTreesFromNonGreedy/DotTrees.dot";
			showNetwork = false;
			
		}
	}
	
	private static void writeTreesToTxtFile(PHYNetwork net, ArrayList<PHYTree> trees, ArrayList<String> sampleNames, Args args) {
		String treeFileName = args.outputFileName;
		try {
			FileWriter fw = new FileWriter(treeFileName);
			fw.write("Nodes:\n" + net.getNodesWithMembersAsString() + "\n");
			for(int i = 0; i < args.numSave; i++) {
				if(trees.size() > i) {
					fw.write("****Tree " + i + "****\n");
					String edges = trees.get(i).toString();
					fw.write(edges);
					fw.write("Error score: " + trees.get(i).getErrorScore()+"\n\n");	
					fw.write("Sample decomposition: \n");
					String lineage = "";
					for(int j = 0; j < sampleNames.size(); j++) {
						lineage += trees.get(i).getLineage(j, sampleNames.get(j));
						lineage += "\n";
					}
					fw.write(lineage);
					fw.write("\n");
				}
			}
			fw.write("SNV info:\n" + net.getNodeMembersOnlyAsString() + "\n");
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + treeFileName);
			System.exit(-1);
		}
	}
	
	private static void writeTreeToDOTFile(String dotTree, ArrayList<String> sampleNames, Args args) {
		String treeFileName = args.outputDOTFileName;
		try {
			FileWriter fw = new FileWriter(treeFileName);
			fw.write(dotTree);
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + treeFileName);
			System.exit(-1);
		}
	}
}
