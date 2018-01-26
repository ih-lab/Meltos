package lineage;

import java.util.ArrayList;
import java.util.Arrays;


public class TreeNode implements Comparable<TreeNode> {
	private boolean isRoot;
	private int nodeId;
	private double[] vafs;
	protected int vafIndices[];
	private ArrayList<SVEntry> svs;
	private int numSamples;
	NodeStats nodeStats;
	static String bProfile;
	
	public TreeNode(int uniqueId, double[] centroidVafs, int[] centIndices, int numS) {
		nodeId = uniqueId;
		vafs = centroidVafs;
		svs = new ArrayList<SVEntry>();
		vafIndices = centIndices;
		//System.out.println(Arrays.toString(centIndices));
		numSamples = numS;
	}
	
	public void setProfile(String bp){
		bProfile = bp;
	}
	
	public void makeStats(int numS){
		nodeStats = new NodeStats(svs, numS);
	}
	
	public TreeNode(int uniqueId) { // Root constructor
		isRoot = true;
		nodeId = uniqueId;
	}

	public boolean isRoot() {
		return isRoot;
	}
	
	public int getNodeId() {
		return nodeId;
	}
	
	/**
	 * Returns the index of this sample in the centroid/AAF data of the group
	 * @return -1 if this sample is not represented in the group
	 */
	public int getSampleIndex(int sampleId) {
		for(int i = 0; i < numSamples; i++) {
			if(vafIndices[i] == sampleId) {
				return i;
			}
		}
		return -1;
	}
	
	public double getVAF(int sampleId) {
		
		if(isRoot) {
			return 0.5;
		} 
		return vafs[sampleId];
	}
	
	public int getVAFLength(){
		return vafs.length;
	}
	
	public void setVAF(double VAF, int index){
		vafs[index]=VAF;
	}
	
	public int getFullSampleLength(){
		return vafIndices.length;
	}
	
	public boolean equals(Object o) {
		if(!(o instanceof TreeNode)) {
			return false;
		}
		TreeNode n = (TreeNode) o;
		if(this.nodeId == n.nodeId) {
			return true;
		} else {
			return false;
		}
	}
	
	public int hashCode() {
		return nodeId;
	}

	@Override
	public int compareTo(TreeNode arg0) {
		if(arg0.getNodeId() < this.nodeId) {
			return -1;
		} else if(arg0.getNodeId() > this.nodeId) {
			return 1;
		}
		return 0;
	}
	
	public SVEntry exportEntry(){
		
		SVEntry asEntry = new SVEntry(vafs.length);//TO MODIFY
		//asEntry.setTemporaryVAFs(vafs, vafIndices);
		asEntry.setTemporaryVAFsFromCondensed(vafs, vafIndices);
		if(bProfile==null)
			asEntry.makeProfile(Parameters.MAX_VAF_ABSENT,Parameters.MIN_VAF_PRESENT);
		else
			asEntry.bProfile=bProfile;
		//asEntry.makeProfileFromCondensed(.1, .1, vafIndices);
		asEntry.assignedNodeId = nodeId;
		asEntry.description = "Centroid";
		return asEntry;
	}
	
	public SVEntry exportEntryNonCompressed(){
		
		SVEntry asEntry = new SVEntry(vafs.length);//TO MODIFY
		//asEntry.setTemporaryVAFs(vafs, vafIndices);
		asEntry.setTemporaryVAFs(vafs);
		if(bProfile==null)
			asEntry.makeProfile(Parameters.MAX_VAF_ABSENT,Parameters.MIN_VAF_PRESENT);
		else
			asEntry.bProfile=bProfile;
		//asEntry.makeProfileFromCondensed(.1, .1, vafIndices);
		asEntry.assignedNodeId = nodeId;
		asEntry.description = "Centroid";
		return asEntry;
	}
	
	public void receiveSV(SVEntry sv){
		svs.add(sv);
	}
	
	public ArrayList<SVEntry> getSVList(){
		return svs;
	}
	
	public void translate(SampleStats[] ssnvStats, SampleStats[] svStats){
		for(int x = 0; x<numSamples; x++){
			if(vafs[x]!=0)
				vafs[x] = svStats[x].convertVAF(vafs[x], ssnvStats[x]);
		}
	}
	
	public void recount(){
		int count = 0;
		for(int x = 0; x<numSamples; x++){
			if(vafs[x]!=0)
				count++;
		}
		numSamples = count;
	}
	
	public void compress(){
		recount();
		double[] newVafs = new double[numSamples];
		int count = 0;
		for(int x = 0; x<vafs.length; x++)
			if(vafs[x]!=0){
				newVafs[count]=vafs[x];
				count++;
			}
		vafs = newVafs;
	}
	public void decompress(){
		if(vafIndices.length==vafs.length)
			return;
		double[] newVafs = new double[vafIndices.length];
		int[] newVafIndices = new int[vafIndices.length];
		for(int x = 0; x<vafs.length; x++){
			newVafs[vafIndices[x]]=vafs[x];
			
		}
		for(int x =0; x<vafIndices.length; x++)
			newVafIndices[x]=x;
		vafs=newVafs;
		vafIndices = newVafIndices;
		numSamples = newVafs.length;
	}
	
	public double[] getVAFs(){
		return vafs;
	}
	
	public void makeNodeStats(){
		nodeStats = new NodeStats(svs, numSamples);
	}
	public String toString(){
		return Arrays.toString(vafIndices);
	}
	
	
}
