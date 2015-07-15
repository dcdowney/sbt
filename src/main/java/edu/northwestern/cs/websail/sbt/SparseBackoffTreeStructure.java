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
	public SparseBackoffTreeStructure(int [] branches, int idx) {
		this(branches, new double[branches.length], 0, 0);
	}
	
	public SparseBackoffTreeStructure(int [] branches, int idx, int startingGlobalIndex) {
		this(branches, new double[branches.length], idx, startingGlobalIndex);
	}
	
	public SparseBackoffTreeStructure(int [] branches, double [] discounts, int idx, int startingGlobalIndex) {
		_minGlobalIndex = startingGlobalIndex;
		_numLeaves = new int[branches[idx]];
		_numLeavesHereAndLeft = new int[branches[idx]];
		_delta = discounts[idx];
		int sum = 0;
		if(idx < branches.length - 1) {
			_children = new SparseBackoffTreeStructure[branches[idx]];
			for(int i=0; i<_children.length; i++) {
				_children[i] = new SparseBackoffTreeStructure(branches, discounts, idx+1, sum + _minGlobalIndex);
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
	
	//assumes all entries of branches are > 0
	public SparseBackoffTreeStructure(int [] branches) {
		this(branches, 0, 0);
	}
	
	public SparseBackoffTreeStructure(int [] branches, double [] discounts) {
		this(branches, discounts, 0, 0);
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
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3});
		Random r = new Random();
		TIntIntHashMap cts = new TIntIntHashMap();
		for(int i=0; i<10000; i++) {
			int d = struct._children[1].randomLeaf(r);
			cts.adjustOrPutValue(d,  1,  1);
		}
		System.out.println("should be roughly uniform from 6 to 11:");
		System.out.println(cts.toString());
		return false;
	}
	
	public static boolean testGetLocalIndex() {
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3});
		int [] lidx = struct.getLocalIndex(3);
		boolean passed = true;
		if(lidx[0]==0) {
			System.out.println("test 1 passed.");
		}
		else {
			System.err.println("test 1 failed: " + lidx[0] + " should be " + 0);
			passed = false;
		}
		if(lidx[1]==3) {
			System.out.println("test 2 passed.");
		}
		else {
			System.err.println("test 2 failed: " + lidx[1] + " should be " + 3);
			passed = false;
		}
		return passed;
		
	}
	
	public static void main(String args[]) {
		testRandomLeaf();
		testGetLocalIndex();
	}
	
}
