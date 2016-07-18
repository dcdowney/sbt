package edu.northwestern.cs.websail.sbt;
import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import junit.framework.Assert;

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
	
	//returns unnormalized smoothed count at r.v. leaf value i
	public double getSmoothed(int i) {
		double [] sc = this.getSmoothAndCount(i);
		return sc[0] + sc[1];
	}
	
	//only performs smoothing if count at leaf is zero 
	public void addAndSmoothIfZero(int i, double mass, double [] ds) {
		double [] smAndCt = this.getSmoothAndCount(i);
		if(smAndCt[1] > 0.0)
			addMass(i, mass);
		else
			smoothAndAddMass(i, mass, ds);
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
		if(_children == null)
			return;
		if(_children[j] == null) {
			_children[j] = new SparseBackoffTree(_struct._children[j]);
		}
		_children[j].addMass(res, d);
	}
	
	public void addAllMass(TIntDoubleHashMap hm, double [] ds) {
		TIntDoubleIterator it = hm.iterator();
		while(it.hasNext()) {
			it.advance();
			if(ds == null)
				addMass(it.key(), it.value());
			else
				smoothAndAddMass(it.key(), it.value(), ds);
		}
	}
	
	public static SparseBackoffTree sum(SparseBackoffTree [] sbts, SparseBackoffTreeStructure struct) {
		SparseBackoffTree out = new SparseBackoffTree(struct);
		for(int i=0; i<sbts.length; i++) {
			out.add(sbts[i]);
		}
		return out;
	}
	
	//adds the given sbt's mass to this sbt
	//assumes shared struct
	public void addWeighted(SparseBackoffTree sbt, double weight) {
		if(sbt != null) {
			this._totalMass += weight*sbt._totalMass;
			this._smooth += weight*sbt._smooth;
			this._delta += weight*sbt._delta;
			for(int i=0; i<this._childMass.length; i++) {
				this._childMass[i] += weight*sbt._childMass[i];
				if(sbt._children != null) {
					if(this._children[i] == null) {
						this._children[i] = new SparseBackoffTree(this._struct._children[i]);
					}
					this._children[i].addWeighted(sbt._children[i], weight);
				}
			}
		}		
	}
	
	//adds the given sbt's mass to this sbt
	//assumes shared struct
	public void add(SparseBackoffTree sbt) {
		addWeighted(sbt, 1.0);
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
			if(sbt == null) {
				out[i+1] = 0.0;
				continue;
			}
			out[i+1] = sbt._childMass[localIdxTrace[i]];
			if(i < out.length - 2)
				sbt = sbt._children[localIdxTrace[i]];
		}
		return out;
	}
	
	public double [] getSmoothsTrace(int [] localIdxTrace) {
		double [] out = new double[localIdxTrace.length];
		SparseBackoffTree sbt = this;
		for(int i=0; i<out.length; i++) {
			if(sbt == null) //can happen if no mass here
				break;
			out[i] = sbt._smooth;
			if(i < out.length - 1)
				sbt = sbt._children[localIdxTrace[i]];
		}
		return out;
	}
	
	/**
	 * Returns a sparse hash map containing leaves with non-zero counts
	 * @return
	 */
	public TIntDoubleHashMap getLeafCounts() {
		TIntDoubleHashMap out = new TIntDoubleHashMap();
		if(this._children != null) {
			for(int i=0; i<_children.length; i++) {
				if(_children[i] != null)
					out.putAll(_children[i].getLeafCounts());
			}
		}
		else {
			for(int i=0; i<_childMass.length; i++) {
				double d = _childMass[i];
				if(d > 0.0) {
					int idx = _struct.getGlobalIndex(i);
					out.put(idx, d);
				}
			}
		}
		return out;
	}
	
	//returns a two-element array {smooth total, count total} for the given leaf
	public double [] getSmoothAndCount(int leafIdx) {
		int [] localIdxTrace = _struct.getLocalIdxTrace(leafIdx);
		double [] out = new double[2];
		SparseBackoffTree sbt = this;
		for(int i=0; i<localIdxTrace.length; i++) {
			if(sbt == null)
				break;
			out[0] += sbt._smooth;
			if(i < localIdxTrace.length - 1)
				sbt = sbt._children[localIdxTrace[i]];
			else {
				out[1] = sbt._childMass[localIdxTrace[i]];
			}
		}
		return out;
	}
	
	//divides leaf counts and re-totals
	//returns new total
	public double divideCountsBy(double [] norm) {
		for(int i=0; i<_childMass.length; i++) {
			if(_childMass[i] > 0.0) {
				double prev = _childMass[i];
				if(_children==null) { //leaf
					int gIdx = _struct.getGlobalIndex(i);
					double n = norm[gIdx];
					_childMass[i] /= n;
				}
				else {
					_childMass[i] = _children[i].divideCountsBy(norm);
				}
				_totalMass -= (prev - _childMass[i]);
			}
		}
		return _totalMass;
	}
	
	public void addInto(double [] target, int startIdx, double toAdd) {
	  if(_children==null) {
	    for(int i=0; i<_childMass.length;i++)
	      target[startIdx++] += toAdd + _childMass[i] + _smooth;
	  }
	  else {
	    toAdd += _smooth;
	    for(int i=0; i<_children.length; i++) {
	      if(_children[i] != null)
	        _children[i].addInto(target, startIdx, toAdd);
	      else
	        for(int j=startIdx;j<startIdx+_struct._numLeaves[i]; j++)
	          target[j] += toAdd;
	      startIdx += _struct._numLeaves[i];
	    }
	  }
	}
	
	public double [] toDoubleArray() {
	  double [] out = new double[_struct.numLeaves()];
	  addInto(out, 0, 0.0);
	  
	  return out;
	}
	
	public static void testDivision() {
		double [] ds = new double [] {0.24, 0.36, 0.3};
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3}, ds);
		SparseBackoffTree sbt = new SparseBackoffTree(struct);
		sbt.smoothAndAddMass(0, 1.0, ds);
		sbt.smoothAndAddMass(3, 2.0, ds);
		sbt.smoothAndAddMass(4, 1.0, ds);
		sbt.smoothAndAddMass(11, 1.0, ds);
		System.out.println("sbt mass: " + sbt._totalMass);
		sbt.divideCountsBy(new double [] {1.0, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
		System.out.println("sbt mass after: " + sbt._totalMass + " should be 4.45");
		TIntIntHashMap hm = sbt.sampleSet(445000, null);
		System.out.println(hm.toString());
		System.out.println(hm.get(1) + " should be about 36000");
		System.out.println(hm.get(3) + " should be about 101000");
		
	}
	
	public static void testSubtraction() {
		double [] ds = new double [] {0.24, 0.36, 0.3};
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3}, ds);
		SparseBackoffTree sbt = new SparseBackoffTree(struct);
		sbt.smoothAndAddMass(0, 1.0, ds);
		sbt.smoothAndAddMass(3, 2.0, ds);
		sbt.smoothAndAddMass(4, 1.0, ds);
		sbt.smoothAndAddMass(11, 1.0, ds);
		SBTSubtractor sub = new SBTSubtractor(sbt, 0, 1.0);
		System.out.println("sbt mass: " + sbt._totalMass);
		TIntIntHashMap hm = sbt.sampleSet(500000, sub);
		System.out.println(hm.toString());
	}
	
	public static void testGetSmoothAndCount() {
		double [] ds = new double [] {0.24, 0.36, 0.3};
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3}, ds);
		SparseBackoffTree sbt = new SparseBackoffTree(struct);
		sbt.smoothAndAddMass(0, 1.0, ds);
		sbt.smoothAndAddMass(3, 2.0, ds);
		sbt.smoothAndAddMass(4, 1.0, ds);
		sbt.smoothAndAddMass(11, 1.0, ds);
		double [] sac = sbt.getSmoothAndCount(3);
		System.out.println("should be {0.46, 1.1} " + Arrays.toString(sac));
	}

	public static void testAddition() {
		double [] ds = new double [] {0.24, 0.36, 0.3};
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3}, ds);
		SparseBackoffTree sbt = new SparseBackoffTree(struct);
		sbt.smoothAndAddMass(0, 1.0, ds);
		sbt.smoothAndAddMass(3, 2.0, ds);
		sbt.smoothAndAddMass(4, 1.0, ds);
		sbt.smoothAndAddMass(11, 1.0, ds);
		SparseBackoffTree sbt2 = new SparseBackoffTree(struct);
		sbt2.addWeighted(sbt, 1.0);
		double [] sac = sbt2.getSmoothAndCount(3);
		System.out.println("should be {0.46, 1.1} " + Arrays.toString(sac));
		sbt2.addWeighted(sbt2, 0.5);
		sac = sbt2.getSmoothAndCount(3);
		System.out.println("should be {0.69, 1.65} " + Arrays.toString(sac));
	}
	
	public static void testToDoubleArray() {
	  double [] ds = new double [] {0.24, 0.36, 0.3};
    SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3}, ds);
    SparseBackoffTree sbt = new SparseBackoffTree(struct);
    sbt.smoothAndAddMass(0, 1.0, ds);
    sbt.smoothAndAddMass(3, 2.0, ds);
    sbt.smoothAndAddMass(4, 1.0, ds);
    sbt.smoothAndAddMass(11, 1.0, ds);
    double [] byGet = new double[12];
    for(int i=0; i<12; i++)
      byGet[i] = sbt.getSmoothed(i);
    double [] byTDA = sbt.toDoubleArray();
    System.out.println("should be equal: ");
    System.out.println(Arrays.toString(byGet));
    System.out.println(Arrays.toString(byTDA));
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
		testDivision();
		testGetSmoothAndCount();
		testAddition();
		testToDoubleArray();
	}
}
