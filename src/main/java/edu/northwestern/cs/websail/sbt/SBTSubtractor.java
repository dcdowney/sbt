package edu.northwestern.cs.websail.sbt;
import java.util.Arrays;


public class SBTSubtractor {
	public double [] _smoother; //amt to subtract from the smoother at node i
	public double [] _total; //amt to subtract from totals
	public int _leafIdx;
	public int [] _localIdx;
	
	//NOTE!!  Performs an approximation to the true subtraction, because takes from leaf-most smoothers even if not all contributed by the skipLeaf 
	public SBTSubtractor(SparseBackoffTree sbt, int skipLeaf) {
		_leafIdx = skipLeaf;

		SparseBackoffTreeStructure struct = sbt._struct;
		_localIdx = struct.getLocalIdxTrace(skipLeaf);
		double [] curSmoothers = struct.getDiscountTrace(_localIdx);
		double [] curTotals = sbt.getTotalsTrace(_localIdx);
		
		_smoother = new double[curSmoothers.length];
		_total = new double[curTotals.length];
		if(curTotals[curTotals.length - 1] < 1.0) {
			double amtSubtracted = curTotals[curTotals.length - 1];
			_total[curTotals.length - 1] = amtSubtracted;
			for(int i=_total.length - 2; i >=0; i--) {
				if(amtSubtracted < 1.0) {
					_smoother[i] = Math.min(Math.min(1.0 - amtSubtracted, curTotals[i]), curSmoothers[i]);
					amtSubtracted += _smoother[i];
				}
				_total[i] = amtSubtracted;
			}
		}
		else {
			Arrays.fill(_total, 1.0); //subtract 1.0 from leaf
		}
		//normalize smoothers to be per-leaf
		for(int i=0; i<_smoother.length - 1;i++) {
			_smoother[i] /= struct.numLeaves();
			struct = struct._children[_localIdx[i]];
		}
		_smoother[_smoother.length - 1] /= struct.numLeaves();
	}
	
	public static SBTSubtractor [] getSubs(SparseBackoffTree [] sbts, int skipLeaf) {
		if(skipLeaf < 0)
			return null;
		SBTSubtractor [] out = new SBTSubtractor[sbts.length];
		for(int i=0; i<sbts.length; i++) {
			out[i] = new SBTSubtractor(sbts[i], skipLeaf);
		}
		return out;
	}
}
