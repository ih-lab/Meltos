package test;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import lineage.LineageTree;
import lineage.SVEntry;
import lineage.SVTreeMapper;
import lineage.TreeNode;

import org.junit.Test;

public class NodeAttatchmentTests {

	@Test
	public void oneNodeOneRootTest() {
		LineageTree tree = new LineageTree(5);
		double[] centroid1 = {0,1,1,1,0};
		double[] centroid2 = {0,0,1,1,0};
		double[] centroid3 = {0,0,0,1,0};
		double[] centroid4 = {0,1,0,1,0};
		double[] centroid5 = {1,0,0,0,1};
		
		int[] indices1 = {1,2,3,0,0};
		int[] indices2 = {2,3,0,0,0};
		int[] indices3 = {3,0,0,0,0};
		int[] indices4 = {1,3,0,0,0};
		int[] indices5 = {0,4,0,0,0};
		
		TreeNode node1 = new TreeNode(1, centroid1, indices1, 3);
		TreeNode node2 = new TreeNode(2, centroid2, indices2, 2);
		TreeNode node3 = new TreeNode(3, centroid3, indices3, 1);
		TreeNode node4 = new TreeNode(4, centroid4, indices4, 2);
		TreeNode node5 = new TreeNode(5, centroid5, indices5, 2);
		
		
		tree.addNode(node1);
		tree.addNode(node2);
		tree.addNode(node3);
		tree.addNode(node4);
		tree.addNode(node5);
		
		tree.addEdge(1, 2);
		tree.addEdge(1, 5);
		tree.addEdge(1, 4);
		tree.addEdge(2, 3);
		

		
		
		double[] vaf1 = {0,.9,.9,.9,0};
		double[] vaf2 = {0,0,.9,.9,0};
		double[] vaf3 = {0,0,0,.9,0};
		double[] vaf4 = {0,.9,0,.9,0};
		double[] vaf5 = {.9,0,0,0,.9};
		
		SVEntry sv1 = new SVEntry(5, "1");
		sv1.setTemporaryVAFs(vaf1);
		SVEntry sv2 = new SVEntry(5, "2");
		sv2.setTemporaryVAFs(vaf2);
		SVEntry sv3 = new SVEntry(5, "3");
		sv3.setTemporaryVAFs(vaf3);
		SVEntry sv4 = new SVEntry(5, "4");
		sv4.setTemporaryVAFs(vaf4);
		SVEntry sv5 = new SVEntry(5, "5");
		sv5.setTemporaryVAFs(vaf5);
		
		ArrayList<SVEntry> svs = new ArrayList<SVEntry>();
		svs.add(sv1);
		svs.add(sv2);
		svs.add(sv3);
		svs.add(sv4);
		svs.add(sv5);
		
		SVTreeMapper mapper = new SVTreeMapper(svs, tree, 5);
		
		mapper.assignSVsViaClustering();
		//mapper.handleUnassigned();
		mapper.newNodeTechnique1();
		System.out.println(tree);
	}
	
	@Test
	public void twoClusterOneGroupRandomComponent() {
		
		LineageTree tree = new LineageTree(5);
		double[] centroid1 = {.5,.5,.5};
		double[] centroid2 = {.11,.11};
		double[] centroid3 = {.101};
		double[] centroid4 = {.101,.101};
		double[] centroid5 = {.101,.101};
		
		int[] indices1 = {1,2,3,0,0};
		int[] indices2 = {2,3,0,0,0};
		int[] indices3 = {3,0,0,0,0};
		int[] indices4 = {1,3,0,0,0};
		int[] indices5 = {0,4,0,0,0};
		
		TreeNode node1 = new TreeNode(1, centroid1, indices1, 3);
		TreeNode node2 = new TreeNode(2, centroid2, indices2, 2);
		TreeNode node3 = new TreeNode(3, centroid3, indices3, 1);
		TreeNode node4 = new TreeNode(4, centroid4, indices4, 2);
		TreeNode node5 = new TreeNode(5, centroid5, indices5, 2);
		
		
		tree.addNode(node1);
		tree.addNode(node2);
		tree.addNode(node3);
		tree.addNode(node4);
		tree.addNode(node5);
		
		tree.addEdge(1, 2);
		tree.addEdge(1, 5);
		tree.addEdge(1, 4);
		tree.addEdge(2, 3);
		
		
//		double[] vaf1 = {0,.398+Math.random()*.004,.398+Math.random()*.004,.398+Math.random()*.004,0};
//		double[] vaf2 = {0,.398+Math.random()*.004,.398+Math.random()*.004,.398+Math.random()*.004,0};
//		double[] vaf3 = {0,.398+Math.random()*.004,.398+Math.random()*.004,.398+Math.random()*.004,0};
//		double[] vaf4 = {0,.398+Math.random()*.004,.398+Math.random()*.004,.398+Math.random()*.004,0};
//		double[] vaf5 = {0,.398+Math.random()*.004,.398+Math.random()*.004,.398+Math.random()*.004,0};
		
//		double[] vaf6 = {0,.198+Math.random()*.004,.198+Math.random()*.004,.198+Math.random()*.004,0};
//		double[] vaf7 = {0,.198+Math.random()*.004,.198+Math.random()*.004,.198+Math.random()*.004,0};
//		double[] vaf8 = {0,.198+Math.random()*.004,.198+Math.random()*.004,.198+Math.random()*.004,0};
//		double[] vaf9 = {0,.198+Math.random()*.004,.198+Math.random()*.004,.198+Math.random()*.004,0};
//		double[] vaf10 = {0,.198+Math.random()*.004,.198+Math.random()*.004,.198+Math.random()*.004,0};
		
		double[] vaf1 = {0,0.400787626931622, 0.3985818654609271, 0.40022483770984896,0};
		double[] vaf2= {0,0.39917297597754536, 0.3990588022537328, 0.3983606782322243,0};
		double[] vaf3 = { 0, 0.3996905567350172, 0.40124533419396025, 0.40056035495865205,0 };
		double[] vaf4 = { 0,0.40028381824791837, 0.4013020582344281, 0.4018243581643057,0 };
		double[] vaf5 = { 0,0.40111578560772765, 0.398680201281717, 0.4009895245476508,0 };
		double[] vaf6 = { 0,0.20155790641963112, 0.2018683128171323, 0.19813824843536146,0 };
		double[] vaf7 = { 0,0.1994139294683819, 0.19923860085874862, 0.19927893161142116,0 };
		double[] vaf8 = { 0,0.19908162465517937, 0.2000590009268309, 0.1996350445707485,0 };
		double[] vaf9 = { 0,0.1980720771378238, 0.19871265617480252, 0.19818930735945253,0 };
		double[] vaf10 = {0,0.1995413662837255, 0.20002159333296435, 0.19836944837092466,0 };
		
		System.out.println(Arrays.toString(vaf1));
		System.out.println(Arrays.toString(vaf2));
		System.out.println(Arrays.toString(vaf3));
		System.out.println(Arrays.toString(vaf4));
		System.out.println(Arrays.toString(vaf5));
		System.out.println(Arrays.toString(vaf6));
		System.out.println(Arrays.toString(vaf7));
		System.out.println(Arrays.toString(vaf8));
		System.out.println(Arrays.toString(vaf9));
		System.out.println(Arrays.toString(vaf10));
		
		SVEntry sv1 = new SVEntry(5);
		sv1.setTemporaryVAFs(vaf1);
		SVEntry sv2 = new SVEntry(5);
		sv2.setTemporaryVAFs(vaf2);
		SVEntry sv3 = new SVEntry(5);
		sv3.setTemporaryVAFs(vaf3);
		SVEntry sv4 = new SVEntry(5);
		sv4.setTemporaryVAFs(vaf4);
		SVEntry sv5 = new SVEntry(5);
		sv5.setTemporaryVAFs(vaf5);
		
		SVEntry sv6 = new SVEntry(5);
		sv6.setTemporaryVAFs(vaf6);
		SVEntry sv7 = new SVEntry(5);
		sv7.setTemporaryVAFs(vaf7);
		SVEntry sv8 = new SVEntry(5);
		sv8.setTemporaryVAFs(vaf8);
		SVEntry sv9 = new SVEntry(5);
		sv9.setTemporaryVAFs(vaf9);
		SVEntry sv10 = new SVEntry(5);
		sv10.setTemporaryVAFs(vaf10);
		
		
		
		ArrayList<SVEntry> svs = new ArrayList<SVEntry>();
		svs.add(sv1);
		svs.add(sv2);
		svs.add(sv3);
		svs.add(sv4);
		svs.add(sv5);
		svs.add(sv6);
		svs.add(sv7);
		svs.add(sv8);
		svs.add(sv9);
		svs.add(sv10);
		
		
		SVTreeMapper mapper = new SVTreeMapper(svs, tree, 5);
		
		mapper.assignSVsViaClustering();
		
		HashMap<Integer, TreeNode> treeNodes = mapper.getNodes();
		
		assertTrue(sv1.getNodeId()==node1.getNodeId());
		assertTrue(sv2.getNodeId()==node1.getNodeId());
		assertTrue(sv3.getNodeId()==node1.getNodeId());
		assertTrue(sv4.getNodeId()==node1.getNodeId());
		assertTrue(sv5.getNodeId()==node1.getNodeId());
		
		//assertEquals(sv6.getNodeId(),6);
		//assertEquals(sv7.getNodeId(),6);
		//assertEquals(sv8.getNodeId(),6);
		//assertEquals(sv9.getNodeId(),6);
		//assertEquals(sv10.getNodeId(),6);
		
		//mapper.handleUnassigned();
		mapper.newNodeTechnique1();
		
		assertTrue(sv6.getNodeId()==6);
		assertTrue(sv7.getNodeId()==6);
		assertTrue(sv8.getNodeId()==6);
		assertTrue(sv9.getNodeId()==6);
		assertTrue(sv10.getNodeId()==6);
		
		System.out.println(tree);
		
	}

}
