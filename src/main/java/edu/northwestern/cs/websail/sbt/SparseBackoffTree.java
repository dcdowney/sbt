package edu.northwestern.cs.websail.sbt;
import gnu.trove.map.hash.TIntIntHashMap;

import java.util.Arrays;
import java.util.Random;


public class SparseBackoffTree {

	SparseBackoffTree [] _children; // if null, children are leaves (= random variable values)
	double [] _childMass;
	
	double _totalMass;
	double _smooth; //prob added to every leaf
	double _delta; //discount subtracted from every leaf
	
	SparseBackoffTreeStructure _struct;
	
	
	public SparseBackoffTree(SparseBackoffTreeStructure struct) {
		_struct = struct;
		_childMass = new double[_struct._numLeaves.length];
		if(_struct._children != null) {
			_children = new SparseBackoffTree[_struct._children.length];
		}
		else {
			_children = null;
		}
	}
	

	
	//for each node at each depth level i,
	//subtracts amount d[i] for all non-zero descendants and re-distributes to all descendants
	public void applySmoother(double [] d) {
		System.err.println("not implemented.");
	}
	
	//adds the given mass after applying the given discounts
	//discounts are provided from shallowest level to deepest
	//NOTE ignores the discounts in the sbt structure!!  Use with care.
	public void smoothAndAddMass(int i, double mass, double [] ds) {
		double [] addAtLevel = new double[ds.length + 1];
		double [] totalBelow = new double[ds.length];
		for(int j=0; j<addAtLevel.length; j++) {
			if(mass > 0.0) {
				if(j<ds.length) {
					addAtLevel[j] = Math.min(mass, ds[j]);
					mass -= ds[j];
				}
				else { //leaf
					addAtLevel[j] = mass;
				}
			}
		}
		totalBelow[totalBelow.length-1] = addAtLevel[totalBelow.length];
		for(int j=totalBelow.length-2;j>=0; j--) {
			totalBelow[j] = totalBelow[j+1] + addAtLevel[j+1];
		}
		addMass(i, totalBelow, addAtLevel, 0);
	}
	
	//adds array of probability mass down to leaf
	//totalBelow is below but not including this level
	//addAtLevel includes the smoothing at this level
	public void addMass(int i, double [] totalBelow, double [] addAtLevel, int idx) {
		int [] locidx = _struct.getLocalIndex(i);
		int j = locidx[0];
		int res = locidx[1];
		_childMass[j] += totalBelow[idx];
		_totalMass += totalBelow[idx] + addAtLevel[idx];
		this._smooth += addAtLevel[idx]/this._struct.numLeaves();
		if(idx<totalBelow.length-1) {
			if(_children[j] == null) {
				_children[j] = new SparseBackoffTree(_struct._children[j]);
			}
			_children[j].addMass(res, totalBelow, addAtLevel, idx+1);
		}
	}
	
	//adds probability at rv (leaf) i
	public void addMass(int i, double d) {
		int [] locidx = _struct.getLocalIndex(i);
		int j = locidx[0];
		int res = locidx[1];
		_childMass[j] += d;
		_totalMass += d;
		if(_children[j] == null) {
			_children[j] = new SparseBackoffTree(_struct._children[j]);
		}
		_children[j].addMass(res, d);
	}
	
	//selects the state at given mass pt
	//HACK: always checks children first, so can remove top smoother by limiting range
	public int selectAt(double d, Random r) {
		int child = -1;
		for(int i=0; i<_childMass.length; i++) {
			if(d<_childMass[i]) {
				child = i;
				break;
			}
			else {
				d-=_childMass[i];
			}
		}
		if(child >= 0) {
			if(_children == null) { //leaf
				return _struct.getGlobalIndex(child);
			}
			else {
				return _children[child].selectAt(d, r);
			}
		}
		else { //smoother -- select uniform below
			int j = _struct.randomLeaf(r);
			//System.out.print(j + ".");
			return j;
			
		}
	}
	
	//selects the state at given mass pt
	public int selectAt(double d, Random r, SBTSubtractor sub, int subIdx) {
		int child = -1;
		for(int i=0; i<_childMass.length; i++) {
			double mass = _childMass[i];
			if(sub._localIdx[subIdx]==i) {
				mass -= sub._total[subIdx+1];
			}
			if(d<mass) {
				child = i;
				break;
			}
			else {
				d-=mass;
			}
		}
		if(child >= 0) {
			if(_children == null) { //leaf
				return _struct.getGlobalIndex(child);
			}
			else {
				if(child==sub._localIdx[subIdx])
					return _children[child].selectAt(d, r, sub, subIdx+1);
				else
					return _children[child].selectAt(d, r);
			}
		}
		else { //smoother -- select uniform below
			int j = _struct.randomLeaf(r);
			//System.out.print(j + ".");
			return j;
		}
	}
	
	//generates a single sample from the tree
	public int sample(Random r) {
		double ch = r.nextDouble()*_totalMass;
		return selectAt(ch, r);
	}
	
	public int sample(Random r, SBTSubtractor sub) {
		return sample(r, sub, 0);
	}
	
	public int sample(Random r, SBTSubtractor sub, int subIdx) {
		double ch = r.nextDouble()*(_totalMass - sub._total[subIdx]);
		return selectAt(ch, r, sub, subIdx);
	}
	
	public int sample(Random r, boolean noTopSmoother) {
		if(noTopSmoother) {
			double ch = r.nextDouble()*(_totalMass - this._smooth * this._struct.numLeaves());
			return selectAt(ch, r);
		}
		else {
			return sample(r);
		}
	}
	
	public int sample(Random r, SBTSubtractor sub, boolean noTopSmoother) {
		return sample(r, sub, 0, noTopSmoother);
	}
	
	public int sample(Random r, SBTSubtractor sub, int subIdx, boolean noTopSmoother) {
		if(noTopSmoother) {
			double ch = r.nextDouble()*(_totalMass - (this._smooth * this._struct.numLeaves() + sub._total[subIdx]));
			return selectAt(ch, r, sub, subIdx);			
		}
		else {
			return sample(r, sub, subIdx);
		}
	}
	
	
	public TIntIntHashMap sampleSet(int numSamples, SBTSubtractor sub) {
		TIntIntHashMap cts = new TIntIntHashMap(); 
		Random r = new Random(4);
		for(int i=0; i<numSamples; i++) {
			int s = -1;
			if(sub==null) 
				s = sample(r);
			else
				s = sample(r, sub);
			cts.adjustOrPutValue(s, 1, 1);
		}
		return cts;
	}
	
	public double [] getTotalsTrace(int [] localIdxTrace) {
		double [] out = new double[localIdxTrace.length + 1];
		SparseBackoffTree sbt = this;
		out[0] = sbt._totalMass;
		for(int i=0; i<out.length - 1; i++) {
			out[i+1] = sbt._childMass[localIdxTrace[i]];
			if(i < out.length - 2)
				sbt = sbt._children[localIdxTrace[i]];
		}
		return out;
	}
	
	public static void testSubtraction() {
		double [] ds = new double [] {0.24, 0.36, 0.3};
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3}, ds);
		SparseBackoffTree sbt = new SparseBackoffTree(struct);
		sbt.smoothAndAddMass(0, 1.0, ds);
		sbt.smoothAndAddMass(3, 2.0, ds);
		sbt.smoothAndAddMass(4, 1.0, ds);
		sbt.smoothAndAddMass(11, 1.0, ds);
		SBTSubtractor sub = new SBTSubtractor(sbt, 0);
		System.out.println("sbt mass: " + sbt._totalMass);
		TIntIntHashMap hm = sbt.sampleSet(400000, sub);
		System.out.println(hm.toString());
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3});
		SparseBackoffTree sbt = new SparseBackoffTree(struct);
		double [] ds = new double [] {0.24, 0.36, 0.3};
		sbt.smoothAndAddMass(0, 1.0, ds);
		sbt.smoothAndAddMass(3, 2.0, ds);
		sbt.smoothAndAddMass(4, 1.0, ds);
		sbt.smoothAndAddMass(11, 1.0, ds);
		System.out.println("sbt mass: " + sbt._totalMass);
		TIntIntHashMap hm = sbt.sampleSet(500000, null);
		System.out.println(hm.toString());
		testSubtraction();
	}
}
