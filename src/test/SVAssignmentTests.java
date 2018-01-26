package test;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import lineage.*;
import static org.junit.Assert.*;
//import static org.junit.assertTrue.*;

import org.junit.Test;

public class SVAssignmentTests {

	@Test
	public void initialTest() {
		
		LineageTree tree = new LineageTree(5);
		double[] centroid1 = {0,1,1,1,0};
		double[] centroid2 = {0,0,1,1,0};
		double[] centroid3 = {0,0,0,1,0};
		double[] centroid4 = {0,1,0,1,0};
		double[] centroid5 = {1,0,0,0,1};
		
		int[] indices = {0,1,2,3,4};
		
		TreeNode node1 = new TreeNode(1, centroid1, indices, 5);
		TreeNode node2 = new TreeNode(2, centroid2, indices, 5);
		TreeNode node3 = new TreeNode(3, centroid3, indices, 5);
		TreeNode node4 = new TreeNode(4, centroid4, indices, 5);
		TreeNode node5 = new TreeNode(5, centroid5, indices, 5);
		
		
		tree.addNode(node1);
		tree.addNode(node2);
		tree.addNode(node3);
		tree.addNode(node4);
		tree.addNode(node5);
		

		
		
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
		
		HashMap<Integer, TreeNode> treeNodes = mapper.getNodes();
		
		assertTrue(sv1.getNodeId()==node1.getNodeId());
		assertTrue(sv2.getNodeId()==node2.getNodeId());
		assertTrue(sv3.getNodeId()==node3.getNodeId());
		assertTrue(sv4.getNodeId()==node4.getNodeId());
		assertTrue(sv5.getNodeId()==node5.getNodeId());
	}
	
	@Test
	public void initialTestWithRandomComponent() {
		
		LineageTree tree = new LineageTree(5);
		double[] centroid1 = {0,.9,.9,.9,0};
		double[] centroid2 = {0,0,.9,.9,0};
		double[] centroid3 = {0,0,0,.9,0};
		double[] centroid4 = {0,.9,0,.9,0};
		double[] centroid5 = {.9,0,0,0,.9};
		
		int[] indices = {0,1,2,3,4};
		
		TreeNode node1 = new TreeNode(1, centroid1, indices, 5);
		TreeNode node2 = new TreeNode(2, centroid2, indices, 5);
		TreeNode node3 = new TreeNode(3, centroid3, indices, 5);
		TreeNode node4 = new TreeNode(4, centroid4, indices, 5);
		TreeNode node5 = new TreeNode(5, centroid5, indices, 5);
		
		
		tree.addNode(node1);
		tree.addNode(node2);
		tree.addNode(node3);
		tree.addNode(node4);
		tree.addNode(node5);
		

		
		
		double[] vaf1 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf2 = {0,0,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf3 = {0,0,0,.7+Math.random()*.2,0};
		double[] vaf4 = {0,.7+Math.random()*.2,0,.7+Math.random()*.2,0};
		double[] vaf5 = {.9,0,0,0,.7+Math.random()*.2};
		
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
		
		ArrayList<SVEntry> svs = new ArrayList<SVEntry>();
		svs.add(sv1);
		svs.add(sv2);
		svs.add(sv3);
		svs.add(sv4);
		svs.add(sv5);
		
		SVTreeMapper mapper = new SVTreeMapper(svs, tree, 5);
		
		mapper.assignSVsViaClustering();
		
		HashMap<Integer, TreeNode> treeNodes = mapper.getNodes();
		
		assertTrue(sv1.getNodeId()==node1.getNodeId());
		assertTrue(sv2.getNodeId()==node2.getNodeId());
		assertTrue(sv3.getNodeId()==node3.getNodeId());
		assertTrue(sv4.getNodeId()==node4.getNodeId());
		assertTrue(sv5.getNodeId()==node5.getNodeId());
		
	}	

	
	@Test
	public void fiveToOneRandomComponent() {
		
		LineageTree tree = new LineageTree(5);
		double[] centroid1 = {0,.9,.9,.9,0};
		double[] centroid2 = {0,0,.9,.9,0};
		double[] centroid3 = {0,0,0,.9,0};
		double[] centroid4 = {0,.9,0,.9,0};
		double[] centroid5 = {.9,0,0,0,.9};
		
		int[] indices = {0,1,2,3,4};
		
		TreeNode node1 = new TreeNode(1, centroid1, indices, 5);
		TreeNode node2 = new TreeNode(2, centroid2, indices, 5);
		TreeNode node3 = new TreeNode(3, centroid3, indices, 5);
		TreeNode node4 = new TreeNode(4, centroid4, indices, 5);
		TreeNode node5 = new TreeNode(5, centroid5, indices, 5);
		
		
		tree.addNode(node1);
		tree.addNode(node2);
		tree.addNode(node3);
		tree.addNode(node4);
		tree.addNode(node5);
		

		
		
		double[] vaf1 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf2 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf3 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf4 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf5 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		
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
		
		ArrayList<SVEntry> svs = new ArrayList<SVEntry>();
		svs.add(sv1);
		svs.add(sv2);
		svs.add(sv3);
		svs.add(sv4);
		svs.add(sv5);
		
		SVTreeMapper mapper = new SVTreeMapper(svs, tree, 5);
		
		mapper.assignSVsViaClustering();
		
		HashMap<Integer, TreeNode> treeNodes = mapper.getNodes();
		
		assertTrue(sv1.getNodeId()==node1.getNodeId());
		assertTrue(sv2.getNodeId()==node1.getNodeId());
		assertTrue(sv3.getNodeId()==node1.getNodeId());
		assertTrue(sv4.getNodeId()==node1.getNodeId());
		assertTrue(sv5.getNodeId()==node1.getNodeId());
		
	}

	@Test
	public void fiveToOneDoubledRandomComponent() {
		
		LineageTree tree = new LineageTree(5);
		double[] centroid1 = {0,.9,.9,.9,0};
		double[] centroid2 = {0,0,.9,.9,0};
		double[] centroid3 = {0,0,0,.9,0};
		double[] centroid4 = {0,.9,0,.9,0};
		double[] centroid5 = {.9,0,0,0,.9};
		
		int[] indices = {0,1,2,3,4};
		
		TreeNode node1 = new TreeNode(1, centroid1, indices, 5);
		TreeNode node2 = new TreeNode(2, centroid2, indices, 5);
		TreeNode node3 = new TreeNode(3, centroid3, indices, 5);
		TreeNode node4 = new TreeNode(4, centroid4, indices, 5);
		TreeNode node5 = new TreeNode(5, centroid5, indices, 5);
		
		
		tree.addNode(node1);
		tree.addNode(node2);
		tree.addNode(node3);
		tree.addNode(node4);
		tree.addNode(node5);
		

		
		
		double[] vaf1 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf2 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf3 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf4 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf5 = {0,.7+Math.random()*.2,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf6 = {0,0,.7+Math.random()*.2,.7+Math.random()*.2,0};
		double[] vaf7 = {0,0,0,.7+Math.random()*.2,0};
		double[] vaf8 = {0,.7+Math.random()*.2,0,.7+Math.random()*.2,0};
		double[] vaf9 = {.9,0,0,0,.7+Math.random()*.2};
		
		
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
		
		SVTreeMapper mapper = new SVTreeMapper(svs, tree, 5);
		
		mapper.assignSVsViaClustering();
		
		HashMap<Integer, TreeNode> treeNodes = mapper.getNodes();
		
		assertTrue(sv1.getNodeId()==node1.getNodeId());
		assertTrue(sv2.getNodeId()==node1.getNodeId());
		assertTrue(sv3.getNodeId()==node1.getNodeId());
		assertTrue(sv4.getNodeId()==node1.getNodeId());
		assertTrue(sv5.getNodeId()==node1.getNodeId());
		assertTrue(sv6.getNodeId()==node2.getNodeId());
		assertTrue(sv7.getNodeId()==node3.getNodeId());
		assertTrue(sv8.getNodeId()==node4.getNodeId());
		assertTrue(sv9.getNodeId()==node5.getNodeId());
		
	}
	
	@Test
	public void twoClusterOneGroupRandomComponent() {
		
		LineageTree tree = new LineageTree(5);
		double[] centroid1 = {0,.9,.9,.9,0};
		double[] centroid2 = {0,0,.9,.9,0};
		double[] centroid3 = {0,0,0,.9,0};
		double[] centroid4 = {0,.9,0,.9,0};
		double[] centroid5 = {.9,0,0,0,.9};
		
		int[] indices = {0,1,2,3,4};
		
		TreeNode node1 = new TreeNode(1, centroid1, indices, 5);
		TreeNode node2 = new TreeNode(2, centroid2, indices, 5);
		TreeNode node3 = new TreeNode(3, centroid3, indices, 5);
		TreeNode node4 = new TreeNode(4, centroid4, indices, 5);
		TreeNode node5 = new TreeNode(5, centroid5, indices, 5);
		
		
		tree.addNode(node1);
		tree.addNode(node2);
		tree.addNode(node3);
		tree.addNode(node4);
		tree.addNode(node5);
		

		
		
		double[] vaf1 = {0,.898+Math.random()*.004,.898+Math.random()*.004,.898+Math.random()*.004,0};
		double[] vaf2 = {0,.898+Math.random()*.004,.898+Math.random()*.004,.898+Math.random()*.004,0};
		double[] vaf3 = {0,.898+Math.random()*.004,.898+Math.random()*.004,.898+Math.random()*.004,0};
		double[] vaf4 = {0,.898+Math.random()*.004,.898+Math.random()*.004,.898+Math.random()*.004,0};
		double[] vaf5 = {0,.898+Math.random()*.004,.898+Math.random()*.004,.898+Math.random()*.004,0};
		
		double[] vaf6 = {0,.598+Math.random()*.004,.598+Math.random()*.004,.598+Math.random()*.004,0};
		double[] vaf7 = {0,.598+Math.random()*.004,.598+Math.random()*.004,.598+Math.random()*.004,0};
		double[] vaf8 = {0,.598+Math.random()*.004,.598+Math.random()*.004,.598+Math.random()*.004,0};
		double[] vaf9 = {0,.598+Math.random()*.004,.598+Math.random()*.004,.598+Math.random()*.004,0};
		double[] vaf10 = {0,.598+Math.random()*.004,.598+Math.random()*.004,.598+Math.random()*.004,0};
		
//		double[] vaf1 = {0.0, 0.9013557500655648, 0.9144289614032401, 0.9189269100463282, 0.0};
//		double[] vaf2= {0.0, 0.89825532532415084, 0.7935381196786324, 0.9132681561402629, 0.0};
//		double[] vaf3 = {0.0, 0.9104818537516357, 0.7925973375118796, 0.89818638154800179, 0.0};
//		double[] vaf4 = {0.0, 0.9081621155395647, 0.9118448765241177, 0.7948230069539983, 0.0};
//		double[] vaf5 = {0.0, 0.7965931309004946, 0.7965076135162366, 0.89845896051145156, 0.0};
//		double[] vaf6 = {0.0, 0.59842828465263851, 0.6118387918033983, 0.6095455497151279, 0.0};
//		double[] vaf7 = {0.0, 0.6113632231354944, 0.6179147200029306, 0.6162112886126917, 0.0};
//		double[] vaf8 = {0.0, 0.59896747993301747, 0.6048060528437476, 0.6127319503602517, 0.0};
//		double[] vaf9 = {0.0, 0.6171614591745778, 0.6109179627540178, 0.5962676246842346, 0.0};
//		double[] vaf10 = {0.0, 0.59816484366861081, 0.59824004253195237, 0.610577334601702, 0.0};
		
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
		
		assertEquals(sv1.getNodeId(),node1.getNodeId());
		assertEquals(sv2.getNodeId(),node1.getNodeId());
		assertEquals(sv3.getNodeId(),node1.getNodeId());
		assertEquals(sv4.getNodeId(),node1.getNodeId());
		assertEquals(sv5.getNodeId(),node1.getNodeId());
		
		assertEquals(sv6.getNodeId(),6);
		assertEquals(sv7.getNodeId(),6);
		assertEquals(sv8.getNodeId(),6);
		assertEquals(sv9.getNodeId(),6);
		assertEquals(sv10.getNodeId(),6);
		
	}
	
	@Test
	public void twoClusterTwoSSNVNodesOneGroupRandomComponent() {
		
		LineageTree tree = new LineageTree(5);
		double[] centroid1 = {0,.9,.9,.9,0};
		double[] centroid2 = {0,0,.9,.9,0};
		double[] centroid3 = {0,0,0,.9,0};
		double[] centroid4 = {0,.9,0,.9,0};
		double[] centroid5 = {.9,0,0,0,.9};
		
		double[] centroid6 = {0,.6,.6,.6,0};
		
		int[] indices = {0,1,2,3,4};
		
		TreeNode node1 = new TreeNode(1, centroid1, indices, 5);
		TreeNode node2 = new TreeNode(2, centroid2, indices, 5);
		TreeNode node3 = new TreeNode(3, centroid3, indices, 5);
		TreeNode node4 = new TreeNode(4, centroid4, indices, 5);
		TreeNode node5 = new TreeNode(5, centroid5, indices, 5);
		
		TreeNode node6 = new TreeNode(6, centroid6, indices, 5);
		
		
		tree.addNode(node1);
		tree.addNode(node2);
		tree.addNode(node3);
		tree.addNode(node4);
		tree.addNode(node5);
		tree.addNode(node6);
		

			

		double[] vaf1 = {0,.898+Math.random()*.004,.898+Math.random()*.004,.898+Math.random()*.004,0};
		double[] vaf2 = {0,.898+Math.random()*.004,.898+Math.random()*.004,.898+Math.random()*.004,0};
		double[] vaf3 = {0,.898+Math.random()*.004,.898+Math.random()*.004,.898+Math.random()*.004,0};
		double[] vaf4 = {0,.898+Math.random()*.004,.898+Math.random()*.004,.898+Math.random()*.004,0};
		double[] vaf5 = {0,.898+Math.random()*.004,.898+Math.random()*.004,.898+Math.random()*.004,0};
		
		double[] vaf6 = {0,.598+Math.random()*.004,.598+Math.random()*.004,.598+Math.random()*.004,0};
		double[] vaf7 = {0,.598+Math.random()*.004,.598+Math.random()*.004,.598+Math.random()*.004,0};
		double[] vaf8 = {0,.598+Math.random()*.004,.598+Math.random()*.004,.598+Math.random()*.004,0};
		double[] vaf9 = {0,.598+Math.random()*.004,.598+Math.random()*.004,.598+Math.random()*.004,0};
		double[] vaf10 = {0,.598+Math.random()*.004,.598+Math.random()*.004,.598+Math.random()*.004,0};
		
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
		
		assertTrue(sv6.getNodeId()==node6.getNodeId());
		assertTrue(sv7.getNodeId()==node6.getNodeId());
		assertTrue(sv8.getNodeId()==node6.getNodeId());
		assertTrue(sv9.getNodeId()==node6.getNodeId());
		assertTrue(sv10.getNodeId()==node6.getNodeId());
		
	}
	

	@Test
	public void clusterSplitterTest() {
		LineageTree tree = new LineageTree(5);
		double[] centroid1 = {0,.901,.901,.901,0};
		double[] centroid2 = {0,0,.9,.9,0};
		double[] centroid3 = {0,0,0,.9,0};
		double[] centroid4 = {0,.9,0,.9,0};
		double[] centroid5 = {.9,0,0,0,.9};
		
		double[] centroid6 = {0,.899,.899,.899,0};
		
		int[] indices = {0,1,2,3,4};
		
		TreeNode node1 = new TreeNode(1, centroid1, indices, 5);
		TreeNode node2 = new TreeNode(2, centroid2, indices, 5);
		TreeNode node3 = new TreeNode(3, centroid3, indices, 5);
		TreeNode node4 = new TreeNode(4, centroid4, indices, 5);
		TreeNode node5 = new TreeNode(5, centroid5, indices, 5);
		
		TreeNode node6 = new TreeNode(6, centroid6, indices, 5);
		
		
		tree.addNode(node1);
		tree.addNode(node2);
		tree.addNode(node3);
		tree.addNode(node4);
		tree.addNode(node5);
		tree.addNode(node6);
		
		double[] vaf1 = {0,.90,.91,.90,0};
		double[] vaf2 = {0,.91,.90,.90,0};
		double[] vaf3 = {0,.90,.90,.91,0};
		double[] vaf4 = {0,.90,.90,.89,0};
		double[] vaf5 = {0,.90,.89,.90,0};
		double[] vaf6 = {0,.89,.90,.90,0};
		
		
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
		SVEntry sv6 = new SVEntry(5, "6");
		sv6.setTemporaryVAFs(vaf6);
		
		ArrayList<SVEntry> svs = new ArrayList<SVEntry>();
		svs.add(sv1);
		svs.add(sv2);
		svs.add(sv3);
		svs.add(sv4);
		svs.add(sv5);
		svs.add(sv6);
		
		SVTreeMapper mapper = new SVTreeMapper(svs, tree, 5);
		
		mapper.assignSVsViaClustering();
		
		HashMap<Integer, TreeNode> treeNodes = mapper.getNodes();
		
		assertTrue(sv1.getNodeId()==node1.getNodeId());
		assertTrue(sv2.getNodeId()==node1.getNodeId());
		assertTrue(sv3.getNodeId()==node1.getNodeId());
		assertTrue(sv4.getNodeId()==node6.getNodeId());
		assertTrue(sv5.getNodeId()==node6.getNodeId());
		
		assertTrue(sv6.getNodeId()==node6.getNodeId());

	}
	
	
	
}

