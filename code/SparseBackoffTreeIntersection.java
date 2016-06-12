package edu.northwestern.cs.websail.sbt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.hash.TIntIntHashMap;

import java.util.Arrays;
import java.util.Random;


public class SparseBackoffTreeIntersection {
	//input sbts:
	SparseBackoffTree [] _sbts;
	boolean _noTopSmooth = false; //only if single SBT
	
	//subtractor for skipleaf:
	SBTSubtractor [] _subs;

	SparseBackoffTreeStructure _struct;
	
	//output intersection structures:
	SparseBackoffTreeIntersection [] _childSbtis;
	SparseBackoffTreeIntersection [] _siblingSbtis;
	double [] _childMass;
	double [] _siblingMass;
	double _totalMass;
	boolean _withoutZeros;
	
	private SparseBackoffTreeIntersection() {
	  
	}
	
	public SparseBackoffTreeIntersection(SparseBackoffTree [] sbtsToIntersect, boolean withoutZeros) {
		this(sbtsToIntersect, -1, null, withoutZeros);
	}
	
	public SparseBackoffTreeIntersection(SparseBackoffTree [] sbtsToIntersect) {
		this(sbtsToIntersect, false);
	}
	
	//if skipLeaf < 1, means no skipLeaf
	public SparseBackoffTreeIntersection(SparseBackoffTree [] sbtsToIntersect, int skipLeaf, double [] amountsToSubtract) {
		this(sbtsToIntersect, skipLeaf, amountsToSubtract, false);
	}
	
	public SparseBackoffTreeIntersection(SparseBackoffTree [] sbtsToIntersect, int skipLeaf, double [] amountsToSubtract, boolean withoutZeros) {
		this(sbtsToIntersect, skipLeaf, 0, amountsToSubtract, withoutZeros);
	}
	
	private SparseBackoffTreeIntersection(SparseBackoffTree [] sbtsToIntersect, int skipLeaf, int smoothMin, double [] amountsToSubtract, boolean withoutZeros) {
		this(sbtsToIntersect, SBTSubtractor.getSubs(sbtsToIntersect, skipLeaf, amountsToSubtract), smoothMin, withoutZeros);
	}
	
	private SparseBackoffTreeIntersection(SparseBackoffTree [] sbtsToIntersect, SBTSubtractor [] subs, int smoothMin, boolean withoutZeros) {
		this(sbtsToIntersect, (subs != null), subs, 0, smoothMin, withoutZeros);
	}
	
	//returns amount to subtract from smooth min
	private int initialize(SparseBackoffTree [] sbtsToIntersect, SBTSubtractor [] subs, int subsIdx, int smoothMin, boolean withoutZeros) {
		_sbts = sbtsToIntersect;
		_subs = subs;
		_withoutZeros = withoutZeros;
		_struct = _sbts[0]._struct;
		
		boolean [] zero = new boolean[_sbts.length];
		int newSize = _sbts.length;
		for(int i=0; i<_sbts.length; i++) {
			SparseBackoffTree sbt = _sbts[i];
			if(sbt._totalMass == 0.0) {
				if(withoutZeros) {
					zero[i] = true;
					newSize--;
				}
				else {
					System.err.println("Zero total mass in SBT, indeterminate results.  Consider using withoutZeros=true.");
				}
			}
		}
		int smoothSubtract = 0;
		if(newSize != _sbts.length) {
			_sbts = new SparseBackoffTree[newSize];
			if(subs != null)
				_subs = new SBTSubtractor[newSize];
			int idx = 0;
			for(int i=0; i<sbtsToIntersect.length; i++) {
				if(sbtsToIntersect[i]._totalMass>0.0) {
					if(subs !=null)
						_subs[idx] = subs[i];
					_sbts[idx++] = sbtsToIntersect[i];
				}
				else if(i < smoothMin)
					smoothSubtract--;
			}
		}
		
		smoothMin -= smoothSubtract;
		_siblingMass = new double[_sbts.length - smoothMin];
		_siblingSbtis = new SparseBackoffTreeIntersection[_sbts.length - smoothMin];
		_childMass = new double[sbtsToIntersect[0]._struct._numLeaves.length];

		return smoothSubtract;
	}
	
	public SparseBackoffTreeIntersection cloneStartingAtSibling(int startSibling) {
	  SparseBackoffTreeIntersection out = new SparseBackoffTreeIntersection();
	  
	  
	  //ugh, is there really no way to invoke a shallow copy from here?
	  //children are shallow copy:
	  out._childMass = this._childMass;
	  out._childSbtis = this._childSbtis;
	  out._noTopSmooth = this._noTopSmooth;
	  out._sbts = this._sbts;
	  out._siblingMass = new double[this._siblingMass.length - startSibling];
	  out._siblingSbtis = new SparseBackoffTreeIntersection[this._siblingSbtis.length - startSibling];
	  out._struct = this._struct;
	  out._subs = this._subs;
	  out._withoutZeros = this._withoutZeros;
	  for(int i = startSibling; i<this._siblingMass.length; i++) {
	    out._siblingMass[i-startSibling] = this._siblingMass[i];
	    out._siblingSbtis[i-startSibling] = this._siblingSbtis[i];
	  }
	  for(int i=0; i<out._childMass.length;i++) {
	    out._totalMass += out._childMass[i];
	  }
	  for(int i=0; i<out._siblingMass.length;i++) {
	    out._totalMass += out._siblingMass[i];
	  }
	  return out;
	}
	
	private SparseBackoffTreeIntersection pruneSiblingsSbti(SparseBackoffTreeIntersection sbti, int startSibling) {
	  return sbti.cloneStartingAtSibling(startSibling);
	}
	
	 private SparseBackoffTreeIntersection(SparseBackoffTree [] sbtsToIntersect, boolean inSub, SBTSubtractor [] subs, int subsIdx, int smoothMin, boolean withoutZeros) {
	    this(sbtsToIntersect, inSub, subs, subsIdx, smoothMin, withoutZeros, new SparseBackoffTreeIntersection [sbtsToIntersect.length]);
	  }
	
	//smoothmin gives the minimum element in sbtsToIntersect that can be smoothed at this level
	//all sbts assumed to have the same number of children (deltas may differ)
	//sbtsToIntersect assumed to have at least one element
	//inSub says whether we're down a branch that should do the subtraction
	private SparseBackoffTreeIntersection(SparseBackoffTree [] sbtsToIntersect, boolean inSub, SBTSubtractor [] subs, int subsIdx, int smoothMin, boolean withoutZeros,
	    SparseBackoffTreeIntersection [] siblingSbtiCache) {
		
		smoothMin -= initialize(sbtsToIntersect, subs, subsIdx, smoothMin, withoutZeros);
		sbtsToIntersect = null;
		subs= null;
		if(_sbts.length==0) {//all have zero mass; ignore
			return;
		}

		if(_sbts.length > 1) {
			for(int i=0; i<_sbts.length - smoothMin; i++) { //try using smoother for sbt i+smoothmin
				if(siblingSbtiCache[i+smoothMin]==null) {
	        SparseBackoffTree [] sbtSubset = new SparseBackoffTree[_sbts.length - 1];
	        SBTSubtractor [] subSubset = null;
	        if(_subs != null) {
	          subSubset = new SBTSubtractor[_subs.length - 1];
	        }
	        int ct = 0;
	        for(int j=0; j<_sbts.length; j++) {
	          if((i + smoothMin)==j) continue;
	          else {
	            if(subSubset != null)
	              subSubset[ct] = _subs[j];
	            sbtSubset[ct++] = _sbts[j];
	          }
	        }
				  _siblingSbtis[i] = new SparseBackoffTreeIntersection(sbtSubset, inSub, subSubset, subsIdx, i + smoothMin, withoutZeros);
				}
				else {
//          SparseBackoffTree [] sbtSubset = new SparseBackoffTree[_sbts.length - 1];
//          SBTSubtractor [] subSubset = null;
//          if(_subs != null) {
//            subSubset = new SBTSubtractor[_subs.length - 1];
//          }
//          int ct = 0;
//          for(int j=0; j<_sbts.length; j++) {
//            if((i + smoothMin)==j) continue;
//            else {
//              if(subSubset != null)
//                subSubset[ct] = _subs[j];
//              sbtSubset[ct++] = _sbts[j];
//            }
//          }
//          SparseBackoffTreeIntersection sbtiTemp = new SparseBackoffTreeIntersection(sbtSubset, inSub, subSubset, subsIdx, i + smoothMin, withoutZeros);
				  _siblingSbtis[i] = siblingSbtiCache[i+smoothMin];
//				  System.out.println("mass of true: " + sbtiTemp._totalMass + " cache: " + _siblingSbtis[i]._totalMass);
				}
				double smoothSubtract = 0.0;
				if(inSub && _subs[i] != null) { //subs[i] may be null if we've set subtract amt on this sbt to zero
					smoothSubtract = _subs[i]._smoother[subsIdx];
				}
				_siblingMass[i] = (_sbts[i + smoothMin]._smooth - smoothSubtract) * _siblingSbtis[i]._totalMass;
				_totalMass += _siblingMass[i];
			}
			double [] restIntersectionMass = null; //optimization -- holds intersection of all SBTs except smoothmin
			double [] prevIntersectionMass = null;
			if(_sbts.length - smoothMin > 0) {
				restIntersectionMass = _siblingSbtis[0]._childMass.clone();
				prevIntersectionMass = _sbts[smoothMin]._childMass.clone();
				if(inSub) {
					prevIntersectionMass[_subs[smoothMin]._localIdx[subsIdx]] -= _subs[smoothMin]._total[subsIdx+1];
				}
			}
			else {
				restIntersectionMass = new double [_sbts[0]._childMass.length];
				Arrays.fill(restIntersectionMass, 1.0);
				prevIntersectionMass = new double[_sbts[0]._childMass.length];
				Arrays.fill(prevIntersectionMass, 1.0);
				for(int j=0; j<=Math.min(smoothMin, _sbts.length - 1); j++) {
					for(int k=0; k<prevIntersectionMass.length; k++) {
						if(inSub && _subs[j]._localIdx[subsIdx]==k) {
							prevIntersectionMass[k] *= _sbts[j]._childMass[k] - _subs[j]._total[subsIdx+1];
						}
						else {
							prevIntersectionMass[k] *= _sbts[j]._childMass[k];
						}
					}
				}
			}
//			double [] prevIntersectionMass = _sbts[smoothMin]._childMass;
//			if(smoothMin == 0) { 
//				prevIntersectionMass = _sbts[smoothMin]._childMass;
//			}
//			else {
//				prevIntersectionMass = _sbts[0]._childMass.clone();
//				for(int j=1; j<=Math.min(smoothMin, _sbts.length - 1); j++) {
//					for(int k=0; k<prevIntersectionMass.length; k++) {
//						prevIntersectionMass[k] *= _sbts[j]._childMass[k];
//					}
//				}
//			}
			
			if(_sbts[0]._children != null) {
				_childSbtis = new SparseBackoffTreeIntersection[restIntersectionMass.length];
				for(int j=0; j<restIntersectionMass.length; j++) {
					if(restIntersectionMass[j] > 0.0 && prevIntersectionMass[j] > 0.0) {
						SparseBackoffTree [] childSbts = new SparseBackoffTree[_sbts.length];
						for(int k=0; k<childSbts.length; k++) {
							childSbts[k] = _sbts[k]._children[j];
						}
						//pass the cache along if it exists, or if we just computed it (smoothMin==0)
            SparseBackoffTreeIntersection [] childSbtiCache = new SparseBackoffTreeIntersection[siblingSbtiCache.length];
//						if(siblingSbtiCache != null && siblingSbtiCache[0] != null) {
//              childSbtiCache = new SparseBackoffTreeIntersection[siblingSbtiCache.length];
//  						for(int k=0; k<childSbtiCache.length; k++) {
//  						  if(siblingSbtiCache[k] != null && siblingSbtiCache[k]._childSbtis != null)
//  						    childSbtiCache[k] = siblingSbtiCache[k]._childSbtis[j];
//  						}
//						}
//						else 
           if(smoothMin==0) {
						  childSbtiCache = new SparseBackoffTreeIntersection[_siblingSbtis.length];
              for(int k=0; k<childSbtiCache.length; k++) {
                if(_siblingSbtis[k] != null && _siblingSbtis[k]._childSbtis != null)
                  childSbtiCache[k] = pruneSiblingsSbti(_siblingSbtis[k]._childSbtis[j], k);
              }
						}
						_childSbtis[j] = new SparseBackoffTreeIntersection(childSbts, (inSub && _subs[0]._localIdx[subsIdx]==j),
							_subs, subsIdx+1, 0, withoutZeros, childSbtiCache);
						_childMass[j] = _childSbtis[j]._totalMass;
					}
					_totalMass += _childMass[j];
				}
			}
			else { //leaf
				for(int j=0; j<restIntersectionMass.length; j++) {
					if(inSub && smoothMin < _subs.length && _subs[smoothMin]._localIdx[subsIdx]==j) {
						_childMass[j] = prevIntersectionMass[j]*restIntersectionMass[j];
					}
					else {
						_childMass[j] = prevIntersectionMass[j]*restIntersectionMass[j];
					}
					_totalMass += _childMass[j];
				}
			}
		}
		else {
			_totalMass = _sbts[0]._totalMass;
			_childMass = _sbts[0]._childMass.clone();
			if(smoothMin > 0) { //no top smoother
				_totalMass -= _sbts[0]._smooth * _sbts[0]._struct.numLeaves();
				_noTopSmooth = true;
			}
			if(inSub && _subs[0] != null) { //subs[0] may be null if we've set subtract amt on this sbt to zero 
				_totalMass -= _subs[0]._total[subsIdx];
				_childMass[_subs[0]._localIdx[subsIdx]] -= _subs[0]._total[subsIdx+1];
			}
		}
	}
	
	public int sample(Random r) {
		if(_sbts.length==0 && _withoutZeros) {
			return _struct.randomLeaf(r);
		}
		double ch = r.nextDouble() * this._totalMass;
		return selectAt(ch, r, (_subs != null)?0:-1);
	}
	
	public int selectAt(double ch, Random r, int subIdx) {
		if(_sbts.length > 1) {
			for(int i=0; i<this._siblingMass.length; i++) {
				if(ch < _siblingMass[i]) {
					//return this._siblingSbtis[i].selectAt(r.nextDouble()*_siblingSbtis[i]._totalMass, r, subIdx);
					return this._siblingSbtis[i].selectAt(ch * (_siblingSbtis[i]._totalMass / _siblingMass[i]), r, subIdx);
					//sibling total mass not on same scale as this sbti's mass, so need to scale random choice for sibling 
				}
				else {
					ch -= _siblingMass[i];
				}
			}
			for(int i=0; i<_childMass.length; i++) {
				if(ch < _childMass[i]) {
					if(_childSbtis == null) { //leaf
						return _sbts[0]._struct.getGlobalIndex(i);
					}
					else {
						if(subIdx >= 0 && i==_subs[0]._localIdx[subIdx])
							return _childSbtis[i].selectAt(ch, r, subIdx+1);
						else
							return _childSbtis[i].selectAt(ch, r, -1);
					}
				}
				else {
					ch -= _childMass[i];
				}
			}
			//we should never get here except if roundoff error, so select uniformly random
			System.err.println("roundoff error");
			return _sbts[0]._struct.randomLeaf(r);
		}
		else {
			if(subIdx >= 0 && _subs[0] != null)
				return _sbts[0].sample(r, _subs[0], subIdx, _noTopSmooth);
			else
				return _sbts[0].sample(r, _noTopSmooth);
		}
	}
	
	public TIntIntHashMap sampleSet(int numSamples) {
		TIntIntHashMap cts = new TIntIntHashMap(); 
		Random r = new Random(4);
		for(int i=0; i<numSamples; i++) {
			int s = sample(r);
			cts.adjustOrPutValue(s, 1, 1);
		}
		return cts;
	}

	public static void testTwo() {
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3});
		SparseBackoffTree sbt = new SparseBackoffTree(struct);
		
		double [] ds = new double [] {0.24, 0.36, 0.3};
		sbt.smoothAndAddMass(0, 1.0, ds);
		sbt.smoothAndAddMass(3, 2.0, ds);
		sbt.smoothAndAddMass(4, 1.0, ds);
		sbt.smoothAndAddMass(11, 1.0, ds);

		SparseBackoffTree sbt2 = new SparseBackoffTree(struct);
		sbt2.smoothAndAddMass(0, 1.0, ds);
		sbt2.smoothAndAddMass(3, 2.0, ds);
		sbt2.smoothAndAddMass(4, 1.0, ds);
		sbt2.smoothAndAddMass(11, 1.0, ds);

		ds = new double [] {1.0, 1.0, 1.0};
		SparseBackoffTree sbt3 = new SparseBackoffTree(struct);
		sbt3.smoothAndAddMass(0, 3.0, ds);
		sbt3.smoothAndAddMass(1, 3.0, ds);
		sbt3.smoothAndAddMass(2, 3.0, ds);
		sbt3.smoothAndAddMass(3, 3.0, ds);
		sbt3.smoothAndAddMass(4, 3.0, ds);
		sbt3.smoothAndAddMass(5, 3.0, ds);
		sbt3.smoothAndAddMass(6, 2.0, ds);
		sbt3.smoothAndAddMass(7, 2.0, ds);
		sbt3.smoothAndAddMass(8, 2.0, ds);
		sbt3.smoothAndAddMass(9, 2.0, ds);
		sbt3.smoothAndAddMass(10, 2.0, ds);
		sbt3.smoothAndAddMass(11, 2.0, ds);
		System.out.println("sbt mass: " + sbt3._totalMass);
		TIntIntHashMap hm = sbt3.sampleSet(300000, null);
		System.out.println(hm.toString());
		
		long t = System.currentTimeMillis();
		SparseBackoffTreeIntersection sbti = new SparseBackoffTreeIntersection(new SparseBackoffTree[] {sbt, sbt2, sbt3});
		System.out.println("intersection mass: " + sbti._totalMass + " should be approx " + 10.868);
		
		double [] targets = new double [] {0.6348,0.3888,0.3888,7.3008,0.9408,0.6348,0.0392,0.0392,0.0392,0.1152,0.1152,0.2312};

		t = System.currentTimeMillis() - t;
		System.out.println("intersection time elapsed: " + t);
		int numSamples = 61650467;
		 t = System.currentTimeMillis();
		 hm = sbti.sampleSet(numSamples);
		t = System.currentTimeMillis() - t;
		System.out.println("sampling time elapsed: " + t);
		System.out.println("samples: " + numSamples);
		System.out.println(hm.toString());
		TIntIntIterator it = hm.iterator();
		int [] cts = new int[12];
		while(it.hasNext()) {
			it.advance();
			cts[it.key()] = it.value();
		}
		
		System.out.println("val\tactual\texpected");
		
		for(int i=0; i<targets.length; i++) {
			double actual = cts[i];
			double expected = (targets[i]*((double)numSamples/10.868));
			System.out.print(i + "\t" + actual + "\t" + expected);
			if((actual - expected)/expected > 0.008)
				System.out.println("\tFAIL");
			else
				System.out.println("\tSUCCESS");
				
		}
			
	}
	
	public static void testZero() {
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3});
		SparseBackoffTree sbt = new SparseBackoffTree(struct);
		
		double [] ds = new double [] {0.24, 0.36, 0.3};
		sbt.smoothAndAddMass(0, 1.0, ds);
		sbt.smoothAndAddMass(3, 2.0, ds);
		sbt.smoothAndAddMass(4, 1.0, ds);
		sbt.smoothAndAddMass(11, 1.0, ds);

		ds = new double [] {1.0, 1.0, 1.0};
		SparseBackoffTree sbt3 = new SparseBackoffTree(struct);
		sbt3.smoothAndAddMass(0, 3.0, ds);
		sbt3.smoothAndAddMass(1, 3.0, ds);
		sbt3.smoothAndAddMass(2, 3.0, ds);
		sbt3.smoothAndAddMass(3, 3.0, ds);
		sbt3.smoothAndAddMass(4, 3.0, ds);
		sbt3.smoothAndAddMass(5, 3.0, ds);
		sbt3.smoothAndAddMass(6, 2.0, ds);
		sbt3.smoothAndAddMass(7, 2.0, ds);
		sbt3.smoothAndAddMass(8, 2.0, ds);
		sbt3.smoothAndAddMass(9, 2.0, ds);
		sbt3.smoothAndAddMass(10, 2.0, ds);
		sbt3.smoothAndAddMass(11, 2.0, ds);
		System.out.println("sbt mass: " + sbt3._totalMass);
		TIntIntHashMap hm = sbt3.sampleSet(300000, null);
		System.out.println(hm.toString());
		
		long t = System.currentTimeMillis();
		SparseBackoffTreeIntersection sbti = new SparseBackoffTreeIntersection(new SparseBackoffTree[] {sbt, sbt3});
		System.out.println("intersection mass: " + sbti._totalMass + " should be approx " + 13.76);

		double [] targets = new double [] {1.38,1.08,1.08,4.68,1.68,1.38,0.28,0.28,0.28,0.48,0.48,0.68};
		
		t = System.currentTimeMillis() - t;
		System.out.println("intersection time elapsed: " + t);
		int numSamples = 1376000;
		 t = System.currentTimeMillis();
		 hm = sbti.sampleSet(numSamples);
		t = System.currentTimeMillis() - t;
		System.out.println("sampling time elapsed: " + t);
		System.out.println("samples: " + numSamples);
		System.out.println(hm.toString());
		TIntIntIterator it = hm.iterator();
		int [] cts = new int[12];
		while(it.hasNext()) {
			it.advance();
			cts[it.key()] = it.value();
		}
		
		System.out.println("val\tactual\texpected");
		
		for(int i=0; i<targets.length; i++) {
			double actual = cts[i];
			double expected = (targets[i]*((double)numSamples/13.76));
			System.out.print(i + "\t" + actual + "\t" + expected);
			if((actual - expected)/expected > 0.008)
				System.out.println("\tFAIL");
			else
				System.out.println("\tSUCCESS");
				
		}
		
	}
	
	public static void testOne() {
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
		SparseBackoffTree sbt2 = new SparseBackoffTree(struct);
		sbt2.smoothAndAddMass(0, 1.0, ds);
		sbt2.smoothAndAddMass(3, 2.0, ds);
		sbt2.smoothAndAddMass(4, 1.0, ds);
		sbt2.smoothAndAddMass(11, 1.0, ds);
		System.out.println("sbt2 mass: " + sbt2._totalMass);
		hm = sbt2.sampleSet(500000, null);
		System.out.println(hm.toString());
		SparseBackoffTree sbt3 = new SparseBackoffTree(struct);
		sbt3.smoothAndAddMass(0, 1.0, ds);
		sbt3.smoothAndAddMass(3, 2.0, ds);
		sbt3.smoothAndAddMass(4, 1.0, ds);
		sbt3.smoothAndAddMass(11, 1.0, ds);
		System.out.println("sbt3 mass: " + sbt3._totalMass);
		hm = sbt3.sampleSet(500000, null);

		SparseBackoffTree sbt4 = new SparseBackoffTree(struct);
		sbt4.smoothAndAddMass(0, 1.0, ds);
		sbt4.smoothAndAddMass(3, 2.0, ds);
		sbt4.smoothAndAddMass(4, 1.0, ds);
		sbt4.smoothAndAddMass(11, 1.0, ds);
		System.out.println("sbt4 mass: " + sbt4._totalMass);
		hm = sbt4.sampleSet(500000, null);
		
		long t = System.currentTimeMillis();
		SparseBackoffTreeIntersection sbti = new SparseBackoffTreeIntersection(new SparseBackoffTree[] {sbt, sbt2, sbt3, sbt4});
		System.out.println("intersection mass: " + sbti._totalMass + " should be approx " + 6.16504672);
		t = System.currentTimeMillis() - t;
		System.out.println("intersection time elapsed: " + t);
		int numSamples = 61650467;
		 t = System.currentTimeMillis();
		hm = sbti.sampleSet(numSamples);
		t = System.currentTimeMillis() - t;
		System.out.println("sampling time elapsed: " + t);
		System.out.println("samples: " + numSamples);
		System.out.println(hm.toString());
		TIntIntIterator it = hm.iterator();
		int [] cts = new int[12];
		while(it.hasNext()) {
			it.advance();
			cts[it.key()] = it.value();
		}
		
		double [] targets = new double [] {0.04477456,0.01679616,0.01679616,5.92240896,0.09834496,
				0.04477456,0.00038416,0.00038416,0.00038416,0.00331776,0.00331776,0.01336336};
		
		System.out.println("val\tactual\texpected");
		
		for(int i=0; i<targets.length; i++) {
			double actual = cts[i];
			double expected = (targets[i]*((double)numSamples/6.16504672));
			System.out.print(i + "\t" + actual + "\t" + expected);
			if((actual - expected)/expected > 0.008)
				System.out.println("\tFAIL");
			else
				System.out.println("\tSUCCESS");
				
		}
		
	}
	
	public static void testThree() {
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
		SparseBackoffTree sbt2 = new SparseBackoffTree(struct);
		sbt2.smoothAndAddMass(0, 1.0, ds);
		sbt2.smoothAndAddMass(3, 2.0, ds);
		sbt2.smoothAndAddMass(4, 1.0, ds);
		sbt2.smoothAndAddMass(11, 1.0, ds);
		System.out.println("sbt2 mass: " + sbt2._totalMass);
		hm = sbt2.sampleSet(500000, null);
		System.out.println(hm.toString());
		SparseBackoffTree sbt3 = new SparseBackoffTree(struct);
		sbt3.smoothAndAddMass(0, 1.0, ds);
		sbt3.smoothAndAddMass(3, 2.0, ds);
		sbt3.smoothAndAddMass(4, 1.0, ds);
		sbt3.smoothAndAddMass(11, 1.0, ds);
		System.out.println("sbt3 mass: " + sbt3._totalMass);
		hm = sbt3.sampleSet(500000, null);
		System.out.println(hm.toString());
		SparseBackoffTree sbt4 = new SparseBackoffTree(struct);
		sbt4.smoothAndAddMass(0, 1.0, ds);
		sbt4.smoothAndAddMass(3, 2.0, ds);
		sbt4.smoothAndAddMass(4, 1.0, ds);
		sbt4.smoothAndAddMass(11, 1.0, ds);
		System.out.println("sbt4 mass: " + sbt4._totalMass);
		hm = sbt4.sampleSet(500000, null);
		System.out.println(hm.toString());
		
		long t = System.currentTimeMillis();
		SparseBackoffTreeIntersection sbti = new SparseBackoffTreeIntersection(new SparseBackoffTree[] {sbt, sbt2, sbt3, sbt4},
				3, new double [] {1.0, 1.0, 1.0, 1.0});
		System.out.println("intersection mass: " + sbti._totalMass + " should be approx " + 0.3409830);
		t = System.currentTimeMillis() - t;
		System.out.println("intersection time elapsed: " + t);
		int numSamples = 34098300;
		 t = System.currentTimeMillis();
		hm = sbti.sampleSet(numSamples);
		t = System.currentTimeMillis() - t;
		System.out.println("sampling time elapsed: " + t);
		System.out.println("samples: " + numSamples);
		System.out.println(hm.toString());
		TIntIntIterator it = hm.iterator();
		int [] cts = new int[12];
		while(it.hasNext()) {
			it.advance();
			cts[it.key()] = it.value();
		}
		
		double [] targets = new double [] {0.04477456,0.01679616,0.01679616,0.09834496,0.09834496,0.04477456,0.00038416,0.00038416,0.00038416,0.00331776,0.00331776,0.01336336};
		
		System.out.println("val\tactual\texpected");
		
		for(int i=0; i<targets.length; i++) {
			double actual = cts[i];
			double expected = (targets[i]*((double)numSamples/0.3409830));
			System.out.print(i + "\t" + actual + "\t" + expected);
			if((actual - expected)/expected > 0.008)
				System.out.println("\tFAIL");
			else
				System.out.println("\tSUCCESS");
				
		}
		
	}

	public static int tinyArrayTest() {
		double [] a = new double [] {1.0 ,2.0, 3.0};
		double [] b = a.clone();
		for(int i=0; i<a.length; i++) {
			a[i] = 5.0;
		}
		System.out.println(Arrays.toString(a));
		System.out.println(Arrays.toString(b));
		return -1;
		
	}
	
	public static void sbtOfZeroTest() {
		SparseBackoffTreeStructure struct = new SparseBackoffTreeStructure(new int [] {2, 2, 3});
		SparseBackoffTree sbt = new SparseBackoffTree(struct);
		
		double [] ds = new double [] {0.24, 0.36, 0.3};
		sbt.smoothAndAddMass(0, 1.0, ds);
		sbt.smoothAndAddMass(3, 2.0, ds);
		sbt.smoothAndAddMass(4, 1.0, ds);
		sbt.smoothAndAddMass(11, 1.0, ds);

		SparseBackoffTree sbt2 = new SparseBackoffTree(struct);
		System.out.println("sbt mass: " + sbt2._totalMass);
		
		ds = new double [] {1.0, 1.0, 1.0};
		SparseBackoffTree sbt3 = new SparseBackoffTree(struct);
		sbt3.smoothAndAddMass(0, 3.0, ds);
		sbt3.smoothAndAddMass(1, 3.0, ds);
		sbt3.smoothAndAddMass(2, 3.0, ds);
		sbt3.smoothAndAddMass(3, 3.0, ds);
		sbt3.smoothAndAddMass(4, 3.0, ds);
		sbt3.smoothAndAddMass(5, 3.0, ds);
		sbt3.smoothAndAddMass(6, 2.0, ds);
		sbt3.smoothAndAddMass(7, 2.0, ds);
		sbt3.smoothAndAddMass(8, 2.0, ds);
		sbt3.smoothAndAddMass(9, 2.0, ds);
		sbt3.smoothAndAddMass(10, 2.0, ds);
		sbt3.smoothAndAddMass(11, 2.0, ds);
		System.out.println("sbt mass: " + sbt3._totalMass);
		TIntIntHashMap hm = sbt3.sampleSet(300000, null);
		System.out.println(hm.toString());
		
		
		long t = System.currentTimeMillis();
		SparseBackoffTreeIntersection sbti = new SparseBackoffTreeIntersection(new SparseBackoffTree[] {sbt, sbt2, sbt3}, true);
		System.out.println("intersection mass: " + sbti._totalMass + " should be approx " + 13.76);

		double [] targets = new double [] {1.38,1.08,1.08,4.68,1.68,1.38,0.28,0.28,0.28,0.48,0.48,0.68};
		
		t = System.currentTimeMillis() - t;
		System.out.println("intersection time elapsed: " + t);
		int numSamples = 1376000;
		 t = System.currentTimeMillis();
		 hm = sbti.sampleSet(numSamples);
		t = System.currentTimeMillis() - t;
		System.out.println("sampling time elapsed: " + t);
		System.out.println("samples: " + numSamples);
		System.out.println(hm.toString());
		TIntIntIterator it = hm.iterator();
		int [] cts = new int[12];
		while(it.hasNext()) {
			it.advance();
			cts[it.key()] = it.value();
		}
		
		System.out.println("val\tactual\texpected");
		
		for(int i=0; i<targets.length; i++) {
			double actual = cts[i];
			double expected = (targets[i]*((double)numSamples/13.76));
			System.out.print(i + "\t" + actual + "\t" + expected);
			if((actual - expected)/expected > 0.008)
				System.out.println("\tFAIL");
			else
				System.out.println("\tSUCCESS");
				
		}		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		testZero();
		testOne();
		testTwo();
		tinyArrayTest();
		testThree();
		sbtOfZeroTest();
	}

}
