package lineage;

import java.util.ArrayList;
import java.util.Arrays;


public class SVEntry {
	
	
	protected String id;
	protected String description;
	protected int[] readSupport;
	protected double[] readDepth;
	protected double[] approxVAF;
	protected int assignedNodeId;
	protected int numSamples;
	protected String bProfile;
	String originalProfile;
	protected int clusterID;
	protected boolean robust;
	
	public SVEntry(int numSamples) {
		this.numSamples = numSamples;
		readSupport = new int[numSamples];
		readDepth = new double[numSamples];
		approxVAF = new double[numSamples];
		assignedNodeId=-1;
		robust = true;
		
	}
	
	public SVEntry(int numSamples, String testID) {
		this.numSamples = numSamples;
		readSupport = new int[numSamples];
		readDepth = new double[numSamples];
		approxVAF = new double[numSamples];
		assignedNodeId=-1;
		id = testID;
	}
	
	public void approximateVAFs() {}
	
	public int getReadSupport(int sampleId) {
		return readSupport[sampleId];
	}
	
	public double getReadDepth(int sampleId) {
		return readDepth[sampleId];
	}
	
	public double getVAF(int sampleId) {
		return approxVAF[sampleId];
	}
		
	public String toString() {
		return id + " " + description;
	}
	
	public void makeProfile(double low, double high){
		low = low;
		high = high;
		int ambiCount = 0;
		String s = "";
		String orS = "";
		for(double x: approxVAF){
			if(x>0)
				orS+="1";
			else
				orS+="0";
			
			if(x>high)
				s+="1";
			else if(x<low)
				s+="0";
			else{
				robust = false;
				String[] stuff = {"1","3","10","12","14","15","17","18","20","22","23","24","26","27"};
				ArrayList<String> missingList = new ArrayList<String>(Arrays.asList(stuff));
				//if(missingList.contains(id))
				//	System.out.println("AHAH!: "+id);
				s+="*";
				ambiCount++;
			}
		}
		bProfile = s;
		originalProfile = orS;
		//if(ambiCount>5)
			//System.out.println("Ambiguities in: "+id+": "+s);
	}
	
	public boolean isRobust(){
		return robust;//TO DO? UNCERTAIN IF THIS ACTUALLY EXISTS IN SNVENTRIES
	}
	
	public String getProfile(){
		return bProfile;
	}
	
	public String getId(){
		return id;
	}
	
	public String getChromosome(){
		return "notImplemented";
	}
	
	public String getPosition(){
		return "notImplemented";
	}
	public String getDescription(){
		return description;
	}
	
	public void updateGroup(String newGroup){
		bProfile = newGroup;
	}
	
	public boolean evidenceOfPresence(int sample){
		return (approxVAF[sample] > Parameters.MAX_VAF_ABSENT );
	}
	
	//stores information to help identify the location of this entry in the group/cluster system.
	public void setCluster(int index){
		clusterID = index;
	}
	
	
	public int getCluster(){
		return clusterID;
	}
	
	///for testing purposes:
	public void setTemporaryVAFs(double[] temps){
		approxVAF = temps;
		//makeProfileFromCondensed(.1, .1, inds);
		makeProfile(Parameters.MAX_VAF_ABSENT,Parameters.MIN_VAF_PRESENT);
	}
	
	public void setNodeId(int id){
		assignedNodeId=id;
	}
	
	public String getEntryId(){
		return id;
	}
	
	//for testing purposes:
	public int getNodeId(){
		return assignedNodeId;
	}
	
	
	//Creates profile when given vafs that don't include 0s.
	public void makeProfileFromCondensed(double low, double high, int[] indices){
		low = low;
		high = high;
		int y = 0;
		String s = "";
		for(int x = 0; x<indices.length; x++){
			if(indices[y]==x){
				if(approxVAF[y]>high)
					s+="1";
				else if(approxVAF[y]<low)
					s+="0";
				else{
					s+="*";
					robust = false;
				}
				y++;
			}
			else
				s+="0";
		}
		bProfile = s;
			
	}
	
	public void setTemporaryVAFsFromCondensed(double[] temps, int[] indices){
		int y = 0;
		approxVAF = new double[indices.length];
		for(int x = 0; x<indices.length; x++){
			if(indices[y]==x){
				approxVAF[x]=temps[y];
				y++;
			}

		}

	}
	
	public void translate(SampleStats[] ssnvStats, SampleStats[] svStats){
		for(int x = 0; x<numSamples; x++){
			if(approxVAF[x]!=0)
				approxVAF[x] = svStats[x].convertVAF(approxVAF[x], ssnvStats[x]);
			if(approxVAF[x]<0)
				approxVAF[x]=0;
			else if(approxVAF[x]>.5)
				approxVAF[x]=.5;
		}
	}
	
	 
}
