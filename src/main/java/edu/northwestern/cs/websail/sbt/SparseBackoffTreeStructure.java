package edu.northwestern.cs.websail.sbt;
import gnu.trove.map.hash.TIntIntHashMap;

import java.util.*;


public class SparseBackoffTreeStructure {
	int [] _numLeaves; //the number of leaf nodes descendant from each child.  All 1s if children are leaves
	int [] _numLeavesHereAndLeft; //convenience.  _numLeavesToLeft[j] = sum of _numLeaves[i] for i<=j

	SparseBackoffTreeStructure [] _children; // if null, children are leaves (= random variable values)
	int _minGlobalIndex; //minimum global index of descendent leaf
	double _delta; //discount subtracted from every leaf and added to this node's smoother
	
	//assumes all entries of branches are > 0
	public SparseBackoffTreeStructure(ArrayList<ArrayList<Integer>> branches, int idx1, int idx2) {
		this(branches, new double[branches.size()], 0, 0, 0);
	}
	
	public SparseBackoffTreeStructure(ArrayList<ArrayList<Integer>> branches, int idx1, int idx2,int startingGlobalIndex) {
		this(branches, new double[branches.size()], idx1, idx2, startingGlobalIndex);
	}
	
	public SparseBackoffTreeStructure(ArrayList<ArrayList<Integer>> branches, double [] discounts, int idx1, int idx2,int startingGlobalIndex) {
		_minGlobalIndex = startingGlobalIndex;
		_numLeaves = new int[branches.get(idx1).get(idx2)];
		_numLeavesHereAndLeft = new int[branches.get(idx1).get(idx2)];
		_delta = discounts[idx1];
		int sum = 0;
		int begin_idx = 0;
		for (int i = 0; i < idx2 ; i++) {
			begin_idx += branches.get(idx1).get(i);
		}
		if(idx1 < branches.size() - 1) {
			_children = new SparseBackoffTreeStructure[branches.get(idx1).get(idx2)];
			for(int i=0; i<_children.length; i++) {
				_children[i] = new SparseBackoffTreeStructure(branches, discounts, idx1+1, begin_idx + i, sum + _minGlobalIndex);
				_numLeaves[i] = _children[i].sumLeaves();
				sum += _numLeaves[i];
				_numLeavesHereAndLeft[i] = sum;
			}
		}
		else {
			Arrays.fill(_numLeaves, 1);
			for(int i=0; i<_numLeavesHereAndLeft.length; i++) {
				_numLeavesHereAndLeft[i] = i+1;
			}
		}
	}

	public SparseBackoffTreeStructure SBT_grow(ArrayList<Integer> new_branch,double new_discount){
		// get the new branching list
		ArrayList<ArrayList<Integer>> branches = this.extract_list();
		int check = 0;
		for(int branching_factor : branches.get(branches.size()-1)){
			check += branching_factor;
		}
		if(check != new_branch.size()){
			System.out.println("Can't expand because of the imcompatible new_branch list");
			return null;
		}
		branches.add(new_branch);
		// get the discount array by checking a random leaf's discount trace and update it to the newest version
		int random_leaf = 0;
		int [] leaf_trace = this.getLocalIdxTrace(random_leaf);
		double []discounts = this.getDiscountTrace(leaf_trace);
		double [] whole_discounts = new double[discounts.length+1];
		for(int i = 0 ; i < discounts.length ; i ++){
			whole_discounts[i] = discounts[i];
		}
		whole_discounts[whole_discounts.length-1] = new_discount;

		SparseBackoffTreeStructure new_root = new SparseBackoffTreeStructure (branches, whole_discounts);
		return new_root;

	}
	// Extract the branching list of a given SBT
	public  ArrayList<ArrayList<Integer>> extract_list(){
		ArrayList<ArrayList<Integer>> result = new ArrayList<>();
		if(this == null) return result;
		
		ArrayList<SparseBackoffTreeStructure> list_of_node = new ArrayList<>();
		list_of_node.add(this);
		
		while(!list_of_node.isEmpty()){
			int size = list_of_node.size();
			ArrayList<Integer> sublist = new ArrayList<Integer>();
			while(size > 0){
				SparseBackoffTreeStructure node = list_of_node.remove(0);
				 if(node._children != null){
				 	sublist.add(node._children.length);
				 	for(int i = 0; i < node._children.length; i++){
				 		list_of_node.add(node._children[i]);
				 	}
				 }
				 // dealing with the second last layer
				 else{
				 	sublist.add(node._numLeaves.length);
				 }
				size--;
			}
			result.add(sublist);
		}
		return result;	
	}
	
	//assumes all entries of branches are > 0
	public SparseBackoffTreeStructure(ArrayList<ArrayList<Integer>> branches) {
		this(branches, 0, 0);
	}
	
	public SparseBackoffTreeStructure(ArrayList<ArrayList<Integer>> branches, double [] discounts) {
		this(branches, discounts, 0,0,0);
	}
	
	public int sumLeaves() {
		int sum = 0;
		for(int i=0; i<_numLeaves.length; i++)
			sum += _numLeaves[i];
		return sum;
	}
	
	//returns {childIndex, leafIndexInChild}
	public int [] getLocalIndex(int leafIndex) {
		int [] out = new int[2];
		// System.out.println(_numLeavesHereAndLeft[0]);
		if(_children == null) {
			out[0] = leafIndex;
			out[1] = 0;
		}
		else {
			for(int i=0; i<_numLeavesHereAndLeft.length; i++) {
				if(_numLeavesHereAndLeft[i] > leafIndex) {
					out[0] = i;
					if(i>0)
						out[1] = leafIndex - _numLeavesHereAndLeft[i-1];
					else 
						out[1] = leafIndex;
					break;
				}
			}
		}
		return out;
	}
	
	public int numLeaves() {
		return _numLeavesHereAndLeft[_numLeavesHereAndLeft.length - 1];
	}
	
	//only meaningful for leaf children
	public int getGlobalIndex(int childIndex) {
		return _minGlobalIndex + childIndex;
	}
	
	public int randomLeaf(Random r) {
		int ch = r.nextInt(numLeaves());
		return ch + _minGlobalIndex;
	}
	//return route in form of ancesters' index given a leaf index
	public int [] getLocalIdxTrace(int leafIndex) {
		ArrayList<Integer> localIdx = new ArrayList<Integer>();
		SparseBackoffTreeStructure struct = this;
		while(struct._children != null) {
			int [] idxs = struct.getLocalIndex(leafIndex);
			int childIdx = idxs[0];
			localIdx.add(childIdx);
			leafIndex = idxs[1];
			struct = struct._children[childIdx];
		}
		// dealing with the second last layer
		int [] idxs = struct.getLocalIndex(leafIndex);
		int childIdx = idxs[0];
		localIdx.add(childIdx);
		
		int [] out = new int[localIdx.size()];
		Iterator<Integer> it = localIdx.iterator();
		for(int i=0; i<out.length; i++) {
			out[i] = it.next();
		}
		return out;
	}
	
	// return the discount list by tracing a leaf's route
	public double [] getDiscountTrace(int [] localIdxTrace) {
		double [] out = new double[localIdxTrace.length];
		SparseBackoffTreeStructure struct = this;
		for(int i=0; i<localIdxTrace.length; i++) {
			out[i] = struct._delta;
			if(i<localIdxTrace.length-1)
				struct = struct._children[localIdxTrace[i]];
		}
		return out;
	}
	
	public static boolean testRandomLeaf() {
		ArrayList<Integer> sublist_1 = new ArrayList<Integer>(Arrays.asList(2));
		ArrayList<Integer> sublist_2 = new ArrayList<Integer>(Arrays.asList(1,2));
		ArrayList<Integer> sublist_3 = new ArrayList<Integer>(Arrays.asList(1,2,3));
		ArrayList<ArrayList<Integer>> list = new ArrayList<ArrayList<Integer>>(Arrays.asList(sublist_1,sublist_2,sublist_3));

		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(list);
		
//		Random r = new Random();
//		TIntIntHashMap cts = new TIntIntHashMap();
//		for(int i=0; i<10000; i++) {
//			int d = struct._children[1].randomLeaf(r);
//			cts.adjustOrPutValue(d,  1,  1);
//		}
		System.out.println("Tests of non_uniform branching SBT constructor");
				
		if(struct._numLeaves[0] == 1 && struct._numLeaves[1] == 5){
			System.out.println("Test 1 passed");
		}else
			System.out.println("Test 1 failed");
		if(struct._children[0]._numLeaves[0] == 1 && struct._children[1]._numLeaves[0] == 2 && struct._children[1]._numLeaves[1] == 3 )
			System.out.println("Test 2 passed");
		else
			System.out.println("Test 2 failed");
		//System.out.println(struct._numLeavesHereAndLeft[1]);
		//System.out.println("should be roughly uniform from 6 to 11:");
		//System.out.println(cts.toString());
		System.out.println("-----------------------");
		return false;
	}
	
	public static boolean testGetLocalIndex() {
		ArrayList<Integer> sublist_1 = new ArrayList<Integer>(Arrays.asList(2));
		ArrayList<Integer> sublist_2 = new ArrayList<Integer>(Arrays.asList(1,2));
		ArrayList<Integer> sublist_3 = new ArrayList<Integer>(Arrays.asList(1,2,3));
		ArrayList<ArrayList<Integer>> list = new ArrayList<ArrayList<Integer>>(Arrays.asList(sublist_1,sublist_2,sublist_3));
		double []discount_ = new double[] {1.1,1.2,1.3};
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(list,discount_);
		int [] lidx = struct.getLocalIndex(3);
		// for(int i = 0 ; i < lidx.length ; i++){
		// 	System.out.println(lidx[i]);
		// }
		System.out.println("Tests of get local index");
		boolean passed = true;
		if(lidx[0] == 1) {
			System.out.println("test 1 passed.");
		}
		else {
			System.err.println("test 1 failed: " + lidx[0] + " should be " + 1);
			passed = false;
		}
		if(lidx[1] == 2) {
			System.out.println("test 2 passed.");
		}
		else {
			System.err.println("test 2 failed: " + lidx[1] + " should be " + 2);
			passed = false;
		}
	    int leaf = 0;
	    int [] trace = struct.getLocalIdxTrace(leaf);
	    double [] discounts = struct.getDiscountTrace(trace);
	    System.out.println("-----------------------");
		return passed;
		
	}
	
	public static boolean testExtract(){
		ArrayList<Integer> sublist_1 = new ArrayList<Integer>(Arrays.asList(2));
		ArrayList<Integer> sublist_2 = new ArrayList<Integer>(Arrays.asList(1,2));
		ArrayList<Integer> sublist_3 = new ArrayList<Integer>(Arrays.asList(1,2,3));
		ArrayList<ArrayList<Integer>> list = new ArrayList<ArrayList<Integer>>(Arrays.asList(sublist_1,sublist_2,sublist_3));
		double []discount_ = new double[] {1.1,1.2,1.3};
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(list,discount_);
		ArrayList<ArrayList<Integer>> list_of_tree = new ArrayList<>();
		list_of_tree = struct.extract_list();
        
		boolean pass_ = list_of_tree.equals(list);
        System.out.println("Test of extract_()");
		if (pass_)
			System.out.println("test passed");
		else
			System.out.print(String.valueOf(list_of_tree) + "should be " + String.valueOf(list));
		// System.out.println(struct.extract_list());
		
		//test SBT_grow
		ArrayList<Integer> sublist_4 = new ArrayList<Integer>();
		int[]num = new int[]{1,2,3,4,5,6};
		for(int i = 0 ; i < num.length ; i++){
			sublist_4.add(num[i]);
		}
		list.add(sublist_4);
		System.out.println(list);
		struct = struct.SBT_grow(sublist_4, 0.5);
		try{
			pass_ = struct.extract_list().equals(list);
		}

		catch (NullPointerException nullpointer){
			System.out.println("Failed to add new branch and test failed.");
			return false;
		}
		System.out.println("Test of SBT_grow()");
		if (pass_)
			System.out.println("test passed.");
		else
			System.out.print("Test failed and" + String.valueOf(struct.extract_list()) + "should be " + String.valueOf(list));
		return false;
	}
	
	public static void main(String args[]) {
		testRandomLeaf();
		testGetLocalIndex();
		testExtract();
	}
	
}
