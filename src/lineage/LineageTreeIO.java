package lineage;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class LineageTreeIO {
	public static LineageTree loadTree(String filename, int numSamples, ArrayList<SVEntry> ssnvs) {
		LineageTree tree = new LineageTree(numSamples);
		try {
			BufferedReader rd = new BufferedReader(new FileReader(filename));
			String currLine = rd.readLine().trim();
			if(currLine == null) returnInputFileFormatError("Empty file", null); 
                        
                        if(currLine.equals("Nodes:")) {
                            currLine = rd.readLine();

                            while (!currLine.equals("")) { // load nodes  
                                tree.addNode(parseTreeNode(currLine, numSamples, ssnvs));				
                                currLine = rd.readLine();
                            }            
                        }
                        
                        currLine = rd.readLine();
                        
			if(!currLine.contains("****")) {
                            returnInputFileFormatError("Did not find the node/edge list separator line #", null);
			}
                        
                        currLine = rd.readLine();
                        
			while (!currLine.contains("Error")) { // load edges
				int[] edge = parseTreeEdge(currLine);
				tree.addEdge(edge[0], edge[1]);				
				currLine = rd.readLine();
			}
			rd.close();
			System.out.println("Loaded " + tree.nodes.size() + " nodes from the input file.\n");
		} catch (IOException e){
			returnInputFileFormatError("Could not read file: " + filename, null);
		}
		return tree;
	}
	
	final static int NUM_NODE_METADATA_FIELDS = 2;
	private static TreeNode parseTreeNode(String line, int numSamples, ArrayList<SVEntry> ssnvs) {	
		String[] entryParts = line.split("\t");
		System.out.println("RAMPAGING FISH!");        
		
		System.out.println("RAMPAGING FISH2!");
		Integer id = Integer.parseInt(entryParts[0].trim());
		String profile = entryParts[1].trim();
		if(profile.length() != numSamples) {
			returnInputFileFormatError("Cluster presence profile " + profile + " length does not match the number of input samples", line);
		}
		/*
		double[] centroid = new double[numSamples];
		for(int i = NUM_NODE_METADATA_FIELDS; i < entryParts.length; i++) {
			try {
                                System.out.println(i);
				if(profile.charAt(i) == '1') {
                                    centroid[i - NUM_NODE_METADATA_FIELDS] = Double.parseDouble(entryParts[i]);
				}
			} catch (NumberFormatException e) {
				returnInputFileFormatError("Cluster centroid VAF value " + entryParts[i], line);
			}
		}
		return new TreeNode(id, centroid);
                */
		double[] centroid = new double[numSamples];
                String vals = entryParts[2].replaceAll("\\[|\\]","");
                vals = vals.trim();
                String[] centroid_vals = vals.split("\\s+");
                int y = 0;
        int[] vafIndices = new int[numSamples];
        int count = 0;
		for(int i = 0; i < numSamples; i++) {
                    
			try {
				if(profile.charAt(i) == '1') {
                                    /*System.out.println(i);
                                    System.out.println(centroid_vals[y]);*/
                                    centroid[i] = Double.parseDouble(centroid_vals[y]);
                                    vafIndices[count]=i;//Part of a workaround.
                                    count++;//workaround continued
                                    /*centroid[i - NUM_NODE_METADATA_FIELDS] = Double.parseDouble(entryParts[i]);*/
                                    y++;
				}
			} catch (NumberFormatException e) {
				returnInputFileFormatError("Cluster centroid VAF value " + entryParts[i], line);
			}
		}
		
		TreeNode tNode = new TreeNode(id, centroid, vafIndices, numSamples);
		tNode.setProfile(profile);
		for(int i = 3; i<entryParts.length; i++)
		{
			tNode.receiveSV(ssnvs.get(Integer.parseInt(entryParts[i].substring(3))-1));
		
		}
		tNode.makeStats(numSamples);
		return tNode;
                
	}
	
	private static int[] parseTreeEdge(String line) {
		String[] entryParts = line.split("->");
		if(entryParts.length != 2) {
			returnInputFileFormatError("Expecting edge line with formant 'a -> b' ", line);
		}
		Integer from = Integer.parseInt(entryParts[0].trim());
		Integer to = Integer.parseInt(entryParts[1].trim());
		return new int[] {from, to};
	}
	
	public static void returnInputFileFormatError(String desc, String entry) {
		System.err.println("[Wrong tree input file format] " + desc);
		if(entry != null) {
			System.err.println("Line: " + entry);
		}
		System.exit(1);
	}
}
