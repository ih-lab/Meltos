package lineage;


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.logging.ConsoleHandler;
import java.util.logging.Formatter;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

//import cc.mallet.util.FileUtils;


public class LapsiEngine {
	static String[] originalProfiles;
	static String[] newProfiles;
	
	
	protected static final Logger logger = Logger.getLogger("lapsi.engine");
	
	/**
	 * The main pipeline for calling SVs using a known lineage tree.
	 */
	public static void callSVs(Args args) throws IOException{
		//System.out.println("TEST!");
		ArrayList<SVEntry> ssnvs = SSNVIO.loadSVFile(args.snvFilePath, args.numSamples);
		//System.out.println("TEST!");
		LineageTree tree = LineageTreeIO.loadTree(args.treeFilePath, args.numSamples, ssnvs);
		//System.out.println("TEST");
		ArrayList<SVEntry> svs = SVIO.loadSVFile(args.svFilePath, args.numSamples, args.calcVAF);
		String[] originals= new String[svs.size()];
		for(int x = 0; x<svs.size(); x++){
			originals[x]=svs.get(x).originalProfile;
			//System.out.println(originals[x]);
		}
		
		String[] news= new String[svs.size()];
		//System.out.println("SVS size "+svs.size());
		
		
		
		// Place the SVs onto the tree.
		SVTreeMapper mapper = new SVTreeMapper(svs, tree, args.numSamples, ssnvs);
		mapper.assignSVsViaClustering();//Testing Clustering technique
		mapper.handleUnassigned();
		//System.out.println("");
		for(String s: mapper.greturn.keySet()){
			for(SVEntry sv: mapper.greturn.get(s).getSVs())
				if(sv.description!="snv"){
					int index = Integer.parseInt(sv.id)-1;
					news[index]=s;
					//System.out.println(news[index]);
					if(!news[index].equals(originals[index]))
						System.out.println(index+" "+originals[index]+" "+news[index]);
				}
		}
		originalProfiles=originals;
		newProfiles=news;
		// Call SVs in each sample.
		
		// Report the results.
		outputSVsToFile(args, mapper, tree);
	}
	
	private static void outputSVsToFile(Args args, SVTreeMapper mapper, LineageTree t) {
		try {
			//FileWriter fw0 = new FileWriter("variedOutput/SV_VAF_Output"+"ERMAR"+Parameters.VAF_ERROR_MARGIN+"_"+"MINV"+Parameters.MAX_VAF_ABSENT+"_MAXV"+Parameters.MIN_VAF_PRESENT+"_COLAPS"+Parameters.MAX_COLLAPSE_CLUSTER_DIFF+".txt");
			//vafDistributionWriter(fw0,mapper);
			//fw0.close();
			
			int goodCount = SVCounterGoodOnly(t,args);
			FileWriter fw = new FileWriter(args.outputFilePath+"Good"+goodCount+""+"_ERMAR"+Parameters.VAF_ERROR_MARGIN+"_"+"MINV"+Parameters.MAX_VAF_ABSENT+"_MAXV"+Parameters.MIN_VAF_PRESENT+"_COLAPS"+Parameters.MAX_COLLAPSE_CLUSTER_DIFF+".txt");
			fw.write("#VAF ERROR:" + Parameters.VAF_ERROR_MARGIN);
			fw.write("#VAF LOWER THRESHOLD:" + Parameters.MAX_VAF_ABSENT);
			fw.write("#VAF UPPER THRESHOLD:" + Parameters.MIN_VAF_PRESENT);
			fw.write("#THRESHOLD FOR MERGING CLUSTERS:"+ Parameters.MAX_COLLAPSE_CLUSTER_DIFF);
			snvNodeWriterDebugger(fw, mapper.tree.nodes);
			svNodeWriterDebugger(fw, mapper.potentialNewNodes, mapper.tree.nodes);
			svInviableNodeWriterDebugger(fw, mapper.potentialNewNodes, mapper.tree.nodes);
			unassignedSVWriter(fw, mapper.potentialNewNodes);
			filteredSVWriter(fw, mapper.unassignedSVs);
			dataCollector(fw, mapper, mapper.tree);
			treePrinter(fw, t);
			fw.close();
			
			if(args.usingThreshold){
				int badandGoodCount = SVCounterBadAndGood(t,args);
				fw = new FileWriter(args.outputFilePath+"BadAndGood"+badandGoodCount+""+"_ERMAR"+Parameters.VAF_ERROR_MARGIN+"_"+"MINV"+Parameters.MAX_VAF_ABSENT+"_MAXV"+Parameters.MIN_VAF_PRESENT+"_COLAPS"+Parameters.MAX_COLLAPSE_CLUSTER_DIFF+".txt");
				fw.write("#VAF ERROR:" + Parameters.VAF_ERROR_MARGIN);
				fw.write("#VAF LOWER THRESHOLD:" + Parameters.MAX_VAF_ABSENT);
				fw.write("#VAF UPPER THRESHOLD:" + Parameters.MIN_VAF_PRESENT);
				fw.write("#THRESHOLD FOR MERGING CLUSTERS:"+ Parameters.MAX_COLLAPSE_CLUSTER_DIFF);
				snvNodeWriterDebugger(fw, mapper.tree.nodes);
				svNodeWriterDebugger(fw, mapper.potentialNewNodes, mapper.tree.nodes);
				svInviableNodeWriterDebugger(fw, mapper.potentialNewNodes, mapper.tree.nodes);
				unassignedSVWriter(fw, mapper.potentialNewNodes);
				filteredSVWriter(fw, mapper.unassignedSVs);
				dataCollector(fw, mapper, mapper.tree);
				profileDifferenceWriter(fw, originalProfiles, newProfiles);
				treePrinter(fw, t);
				fw.close();
			}
		
		
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + args.outputFilePath);
			System.exit(-1);
		}
		
	}
	
	
	
	private static void vafDistributionWriter(FileWriter fw, SVTreeMapper mapper)throws IOException{
		fw.write("#SV VAF distributions for each sample\n");
		SampleStats[] svStats = mapper.svStats;
		for(SampleStats s: svStats){
			String meanVal = ""+s.mean;
			String varVal = ""+Math.pow(s.std,2);
			meanVal = meanVal.substring(0,Math.min(5,meanVal.length()));
			varVal = varVal.substring(0,Math.min(5,varVal.length()));
			
			fw.write("Mean: "+meanVal+" Variance: "+varVal+"|");
		}
		fw.write("\n#SV VAFs\n");
		for(SVEntry s: mapper.svs){
			for(double numVal: s.approxVAF){
				String val = ""+numVal;
				val = val.substring(0, Math.min(5,val.length()));
				fw.write(val+"\t");
			}
			fw.write("\n");
		}
		
		
		
	}
	
	private static void snvNodeWriter(FileWriter fw, HashMap<Integer, TreeNode> nodes)throws IOException{
		fw.write("#Lineage Tree Node Assignments:\n");
		fw.write("Node_ID\tSV_IDs\n");
		for(Integer i: nodes.keySet()){
			
			TreeNode t = nodes.get(i);
			String s = ""+t.getNodeId()+"\t";
			for(SVEntry sv: t.getSVList()){
				if(!sv.description.equals("Centroid"))
					s+= sv.getEntryId()+",";
				
			}
			fw.write(s+"\n");
		}
		
	}
	
	private static void svNodeWriter(FileWriter fw, ArrayList<TreeNode> potentialNewNodes, HashMap<Integer, TreeNode> nodes)throws IOException{
		fw.write("#SVs Assigned to Added New SV Nodes:\n");
		fw.write("Node_ID\tSV_IDs\n");
		for(TreeNode t: potentialNewNodes){
			if(!nodes.containsKey(t.getNodeId()))
					continue;
			String s = ""+t.getNodeId()+"\t";
			for(SVEntry sv: t.getSVList()){
				if(!sv.description.equals("Centroid"))
					s+= sv.getEntryId()+",";
				
			}
			fw.write(s+"\n");
		}
	}
	
	private static void svInviableNodeWriter(FileWriter fw, ArrayList<TreeNode> potentialNewNodes, HashMap<Integer, TreeNode> nodes)throws IOException{
		fw.write("#SVs Assigned to Inviable New SV Nodes:\n");
		fw.write("Node_ID\tSV_IDs\n");
		for(TreeNode t: potentialNewNodes){
			if(nodes.containsKey(t.getNodeId()))
					continue;
			String s = ""+t.getNodeId()+"\t";
			for(SVEntry sv: t.getSVList()){
				if(!sv.description.equals("Centroid"))
					s+= sv.getEntryId()+",";
				
			}
			fw.write(s+"\n");
		}
	}
	
	private static void unassignedSVWriter(FileWriter fw, ArrayList<TreeNode> potentialNewNodes)throws IOException{
		fw.write("#Unassigned SVs:\n");
		
		fw.write("Not implemented yet\n");
	}
	
	private static void filteredSVWriter(FileWriter fw, ArrayList<SVEntry> unassignedSVs)throws IOException{
		fw.write("#SVs filtered due to VAF insignificance:\n");
		fw.write("Not implemented yet\n");
	}
	
	
	private static void profileDifferenceWriter(FileWriter fw, String[] orgs, String[] news)throws IOException{
		fw.write("Differences in profiles during clustering\n");
		fw.write("Node_ID\toriginal\tmodified");
		for(int x = 0; x<orgs.length; x++){
			if(!orgs[x].equals(news[x]))
				fw.write((x+1)+"\t"+orgs[x]+"\t"+news[x]+"\n");;
				
			}
			
		
	}
	
	// ---- LAUNCH ----
	public static void main(String[] args)throws IOException {
		Options options = new Options(); 
		
		// Input/Output
		options.addOption("treeFile", "t", true, "Input lineage tree file path [required]");
		options.addOption("svFile", "sv", true, "Input SV file path [required]");
		options.addOption("numSamples", "n", true, "Number of input samples (must be consistent with both input files) [required]");
		options.addOption("outputFile", "o", true, "Output file path (default: prefix is the input SV file path)");
                
                //Calculate VAF from Genomic Counts
                options.addOption("calcVAF", "vaf", false, "Calculate VAF from counts provided in file");
                
		// Tree placement (TODO: similarity)
		options.addOption("allowSVOnlyNodes", true, "Allows the placement of SVs to new tree nodes, as long as they satisfy the PP constraints.");
	
		// SV filtering
		//options.addOption("minDepth", true, "Minimum number of reads supporting the event.");
				
		options.addOption("v", "verbose", false, "Verbose mode");
		options.addOption("h", "help", false, "Print usage");
		
		options.addOption("ssnvFile", "snv", true, "Input the file containing the SVs from the lineageTree");//experimental
		options.addOption("negativeCutoff", "nc", true, "Give a threshold for false SVs for testing purpose. SVs indexed after that threshold are considered fake.");//experimental
		
		// Set the display order.
		ArrayList<Option> optionsList = new ArrayList<Option>();
		optionsList.add(options.getOption("treeFile"));
		optionsList.add(options.getOption("svFile"));
		optionsList.add(options.getOption("numSamples"));
		optionsList.add(options.getOption("outputFile"));
                optionsList.add(options.getOption("calcVAF"));
		optionsList.add(options.getOption("allowSVOnlyNodes"));
		optionsList.add(options.getOption("v"));
		optionsList.add(options.getOption("h"));
		optionsList.add(options.getOption("ssnvFile"));
		optionsList.add(options.getOption("negativeCutoff"));
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmdLine = null;
		HelpFormatter hf = new HelpFormatter();
		hf.setOptionComparator(new OptionComarator<Option>(optionsList));
		try {
			cmdLine = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			hf.printHelp("lapsi", options);
			System.exit(-1);
		}
		
		// Set-up the input parameters.
		Args params = new Args();	
		if(cmdLine.hasOption("treeFile")) {
			params.treeFilePath = cmdLine.getOptionValue("treeFile");
		} else {
			System.out.println("Required parameter: input lineage tree file path [-treeFile]");
			hf.printHelp("lapsi", options);
			System.exit(-1);
		}
		if(cmdLine.hasOption("svFile")) {
			params.svFilePath = cmdLine.getOptionValue("svFile");
		} else {
			System.out.println("Required parameter: input SV file path [-svFile]");
			hf.printHelp("lapsi", options);
			System.exit(-1);
		}
		if(cmdLine.hasOption("numSamples")) {
			params.numSamples = Integer.parseInt(cmdLine.getOptionValue("numSamples"));
		} else {
			System.out.println("Required parameter: number of input samples [-numSamples]");
			hf.printHelp("lapse", options);
			System.exit(-1);
		}	
		if(cmdLine.hasOption("outputFile")) {
			params.outputFilePath = cmdLine.getOptionValue("outputFile");	
		} else {
			params.outputFilePath = params.svFilePath + ".lapsi";
		}
                
                if(cmdLine.hasOption("calcVAF")) {
			params.calcVAF = true;
		}
		
		if(cmdLine.hasOption("allowSVOnlyNodes")) {
			params.allowSVOnlyNodes = true;
		}	
		
		if(cmdLine.hasOption("h")) {
			new HelpFormatter().printHelp(" ", options);
		}
		
		if(cmdLine.hasOption("negativeCutoff")) {
			params.negativeThreshold=Integer.parseInt(cmdLine.getOptionValue("negativeCutoff"));
			params.usingThreshold=true;
		}
		// Setup the logger verbosity level.
		ConsoleHandler h = new ConsoleHandler();
		h.setFormatter(new LogFormatter());
		h.setLevel(Level.INFO);
		logger.setLevel(Level.INFO);
		if(cmdLine.hasOption("v")) {
			h.setLevel(Level.FINEST);
			logger.setLevel(Level.FINEST);
		}
		
		//Experimental
		if(cmdLine.hasOption("ssnvFile")){
			params.snvFilePath = cmdLine.getOptionValue("ssnvFile");
		}
		else{
			System.out.println("Required parameter: snvFile");
			hf.printHelp("lapse", options);
			System.exit(-1);
		}
		
		logger.addHandler(h);
		logger.setUseParentHandlers(false);
		File fileFolder = new File(params.outputFilePath);
		fileFolder.mkdir();
		File[] files = fileFolder.listFiles();
	    if(files!=null) { //some JVMs return null for empty dirs
	        for(File f: files) {
	                f.delete();
	        }
	    }
	    
	    //VAF ERROR:{0.02};//
	    //VAF LOWER THRESHOLD:{0.05};//
	    //VAF UPPER THRESHOLD:{0.3};//
	    //THRESHOLD FOR MERGING CLUSTERS:{0.2};//

		// Call into the main pipeline.#VAF ERROR:0.08#VAF LOWER THRESHOLD:0.3#VAF UPPER THRESHOLD:0.4#THRESHOLD FOR MERGING CLUSTERS:0.3
		double vafMax[] = {.05};///{.05,.1,.2,.3,.4};///
		double vafMin[] = {.05};////{.05,.1,.2,.3,.4};//
		double vafEr[] = {.02};////{.02,.04,.06,.08,.1};//
		double vafClusDif[] ={.05};////{.05,.1,.2,.3,.4};//
		for(double vM: vafMax)
			for(double vMin: vafMin)
				for(double vE: vafEr)
					for(double vCD: vafClusDif)
						if(vM>=vMin)
						{
							Parameters.MAX_VAF_ABSENT=vMin;
							Parameters.MIN_VAF_PRESENT=vM;
							Parameters.VAF_ERROR_MARGIN=vE;
							Parameters.MAX_COLLAPSE_CLUSTER_DIFF=vCD;
							
							callSVs(params);
						}
	}
	
	protected static class Args {
		String treeFilePath;
		String svFilePath;
		String snvFilePath;//Experimental
		int numSamples;
		String outputFilePath;
		int negativeThreshold;
		boolean usingThreshold = false;
		boolean allowSVOnlyNodes = false;
        boolean calcVAF = false;
	}

	protected static class LogFormatter extends Formatter {
		public String format(LogRecord rec) {
			return rec.getMessage() + "\r\n";
		}
	}
	
	protected static class OptionComarator<T extends Option> implements Comparator<T> {
	    protected ArrayList<Option> orderedOptions;
	    public OptionComarator(ArrayList<Option> options) {
	    	orderedOptions = options;
	    }
	    public int compare(T o1, T o2) {
	        return orderedOptions.indexOf(o1) - orderedOptions.indexOf(o2);
	    }
	}
	
	
	
	private static void snvNodeWriterDebugger(FileWriter fw, HashMap<Integer, TreeNode> nodes)throws IOException{
		fw.write("#Lineage Tree Node Assignments:\n");
		fw.write("Node_ID\tSV_IDs\n");
		for(Integer i: nodes.keySet()){
			
			TreeNode t = nodes.get(i);
			String s = ""+t.getNodeId()+"\t";
			int goodcount = 0;
			int badcount = 0;
			
			for(SVEntry sv: t.getSVList())
				if(!sv.description.equals("Centroid") && !sv.description.equals("snv"))
					if((Integer.parseInt(sv.getEntryId())>1000)||(Integer.parseInt(sv.getEntryId())<34))
						goodcount++;
					else if(Integer.parseInt(sv.getEntryId())>=34)
						badcount++;
			fw.write("ValidSV"+goodcount+"\t");
			fw.write("invalidSV"+badcount+"\t");
			if(t.getSVList().size()!=0)
				s+= t.bProfile+"\t";//t.getSVList().get(0).bProfile+"\t";
			s+= Arrays.toString(t.getVAFs())+"\t";
			for(SVEntry sv: t.getSVList()){
				if(!sv.description.equals("Centroid") && !sv.description.equals("snv"))
					s+= sv.getEntryId()+",";
				
			}
			fw.write(s+"\n");
		}
		
	}
	
	private static void svNodeWriterDebugger(FileWriter fw, ArrayList<TreeNode> potentialNewNodes, HashMap<Integer, TreeNode> nodes)throws IOException{
		fw.write("#SVs Assigned to Added New SV Nodes:\n");
		fw.write("Node_ID\tSV_IDs\n");
		for(TreeNode t: potentialNewNodes){
			if(!nodes.containsKey(t.getNodeId()))
					continue;
			String s = ""+t.getNodeId()+"\t";
			s+= t.bProfile+"\t";//t.getSVList().get(0).bProfile+"\t";
			s+= Arrays.toString(t.getVAFs())+"\t";
			int goodcount = 0;
			int badcount = 0;;
			for(SVEntry sv: t.getSVList())
				if(!sv.description.equals("Centroid") && !sv.description.equals("snv"))
					if((Integer.parseInt(sv.getEntryId())>1000)||(Integer.parseInt(sv.getEntryId())<34))
						goodcount++;
					else if(Integer.parseInt(sv.getEntryId())>=34)
						badcount++;
			//System.out.println("ValidSV"+goodcount+"\t");
			//System.out.println("invalidSV"+badcount+"\t");
			for(SVEntry sv: t.getSVList()){
				if(!sv.description.equals("Centroid") && !sv.description.equals("snv"))
					s+= sv.getEntryId()+",";
				
			}
			fw.write(s+"\n");
		}
	}
	
	private static void svInviableNodeWriterDebugger(FileWriter fw, ArrayList<TreeNode> potentialNewNodes, HashMap<Integer, TreeNode> nodes)throws IOException{
		fw.write("#SVs Assigned to Inviable New SV Nodes:\n");
		fw.write("Node_ID\tSV_IDs\n");
		for(TreeNode t: potentialNewNodes){
			if(nodes.containsKey(t.getNodeId()))
					continue;
			String s = ""+t.getNodeId()+"\t";
			//if(t.getSVList().size()==0)
				//System.out.println("test");
			s+= t.bProfile+"\t";//s+= t.getSVList().get(0).bProfile+"\t";
			s+= Arrays.toString(t.getVAFs())+"\t";
			for(SVEntry sv: t.getSVList()){
				if(!sv.description.equals("Centroid") && !sv.description.equals("snv"))
					s+= sv.getEntryId()+",";
				
			}
			fw.write(s+"\n");
		}
	}
	
	private static void treePrinter(FileWriter fw, LineageTree t)throws IOException{
		fw.write(t.toString());
	}
	
	private static void dataCollector(FileWriter fw,SVTreeMapper mapper, LineageTree tree)throws IOException{
		int addedNodeSVCount = 0;
		int unaddedNodeSVCount = 0;
		int insignificantNodeSVCount = 0;
		int filteredSVs = 0;
		
		for(SVEntry s: mapper.svs){
			if(s.assignedNodeId==-2)
				filteredSVs++;
		}
		for(TreeNode t: mapper.potentialNewNodes){
			if(t.getSVList().size()>1)
				unaddedNodeSVCount+=t.getSVList().size();
			else
				insignificantNodeSVCount+=t.getSVList().size();
		}
		for(Integer i: tree.nodes.keySet()){
			
			TreeNode t = tree.nodes.get(i);
			for(SVEntry sv: t.getSVList())
				if(!sv.description.equals("Centroid"))
					addedNodeSVCount++;
			//addedNodeSVCount+=t.getSVList().size();
			
		}
		
		fw.write("addedNodeSVCount"+addedNodeSVCount+"\n");
		fw.write("unaddedNodeSVCount"+unaddedNodeSVCount+"\n");
		fw.write("insignificantNodeSVCount"+insignificantNodeSVCount+"\n");
		fw.write("filteredSVs"+filteredSVs+"\n");
		fw.write("(addedNodeSVCount+unaddedNodeSVCount+insignificantNodeSVCount+filteredSVs)"+(addedNodeSVCount+unaddedNodeSVCount+insignificantNodeSVCount+filteredSVs)+"\n");
		
		//System.out.println(unaddedNodeSVCount);
		//System.out.println(insignificantNodeSVCount);
		//System.out.println(filteredSVs);
		//System.out.println("Sum "+ (addedNodeSVCount+unaddedNodeSVCount+insignificantNodeSVCount+filteredSVs));
		
		int negCount = 0;
		//for(SVEntry s: mapper.svs){
			//System.out.println(s.id+" "+s.assignedNodeId);
			//if(s.assignedNodeId==-1||s.assignedNodeId==-2){
			//	System.out.println(s.bProfile);
			//}
			
		
		//}
		//System.out.println(negCount);
	}
	
	private static int SVCounterGoodOnly(LineageTree t, Args args){
		int count=0;
		if(args.usingThreshold){
			for(int nodeInd: t.nodes.keySet()){
				TreeNode tNode = t.nodes.get(nodeInd);
				for(SVEntry sv: tNode.getSVList())
					if(!sv.description.equals("snv"))
						if((Integer.parseInt(sv.getEntryId())>1000&&tNode.getNodeId()<5)||(!sv.description.equals("snv")&&Integer.parseInt(sv.getEntryId())<args.negativeThreshold))
							count++;
			}
		}
		else{
			for(int nodeInd: t.nodes.keySet()){
				TreeNode tNode = t.nodes.get(nodeInd);
				for(SVEntry sv: tNode.getSVList())
					if(!sv.description.equals("snv"))
						if((Integer.parseInt(sv.getEntryId())>1000&&tNode.getNodeId()<5)||(!sv.description.equals("snv")))//&&Integer.parseInt(sv.getEntryId())<args.negativeThreshold))
							count++;
			}
		}
		return count;
		
	}
	private static int SVCounterBadAndGood(LineageTree t, Args args){
		int count=0;
		for(int nodeInd: t.nodes.keySet()){
			TreeNode tNode = t.nodes.get(nodeInd);
			for(SVEntry sv: tNode.getSVList()){
				if(!sv.description.equals("snv")){
					if((Integer.parseInt(sv.getEntryId())>1000&&tNode.getNodeId()<5)||(!sv.description.equals("snv")&&Integer.parseInt(sv.getEntryId())<args.negativeThreshold))
						count++;
					else if(!sv.description.equals("snv")&&Integer.parseInt(sv.getEntryId())>=args.negativeThreshold)
						count--;
				}
					
			}
		}
		return count;
	}
	
	
}
