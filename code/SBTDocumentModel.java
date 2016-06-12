package edu.northwestern.cs.websail.sbt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.lang.reflect.Method;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.Iterator;
import java.util.Random;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import java.util.Map;
import java.util.Map.Entry;
import java.util.HashMap;

import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntDoubleHashMap;

public class SBTDocumentModel implements Serializable {

	private static final long serialVersionUID = 2L;

	private class GibbsDoer implements Runnable {
    	int _start;
    	int _end;
    	long _changes = 0;
    	
    	public void run() {
    		_changes = sampleRange(_start, _end);
    	}
    }
	
    private class TestDoer implements Runnable {
    	int _doc;
    	double [] _ds;
    	double [] _wordTopicMarginal;
    	double [] _res;
    	
    	public TestDoer(int i, double [] ds, double [] wordTopicMarginal) {
    		_doc = i;
    		_ds = ds;
    		_wordTopicMarginal = wordTopicMarginal;
    	}
    	
    	public void run() {
    		_res = testOnDocFullPpl(_doc, _ds, _wordTopicMarginal);
    		System.out.println(_doc + "\t" + _res[0] + "\t" + _res[1]);
    	}
    }
	
	Random _r = new Random();
	
	private Corpus _c;
	
	//Model information:
	int [] _branchingFactors;
	SparseBackoffTreeStructure _struct = null;
	SparseBackoffTree [] _topicGivenWord = null; 
	double [] _topicMarginal;
	double [] _docsDelta = null;
	double [] _wordsDelta = null;
	

	//threading:
	private int _NUMTHREADS = -1;
	private int[] _THREADBREAKS = null; //inclusive doc indices where threads *end* (initial thread starts
							//at 0 implicitly)
	
	//training config:
	private int _NUMITERATIONS = -1;
	private int _OPTIMIZEINTERVAL = -1;
	private int [] _EXPANSIONSCHEDULE = null; //specifies indexes of branching factor to introduce incrementally
	private int _USEEXPANSION = 1; //1 for yes, 0 for no
	private int [] _expansionBranchFactors = null; //when expansion is used, holds the full branching factors
	
	//testing config:
	private int _NUMPARTICLES = -1;
	private int _MAXTESTDOCSIZE = -1;
	
	public SBTDocumentModel(String configFile) throws Exception {
		_c = new Corpus();
		readConfigFile(configFile);
		if(_USEEXPANSION == 1) {
			if(_EXPANSIONSCHEDULE == null) {
				_EXPANSIONSCHEDULE = defaultExpansion(_branchingFactors.length);
			}
			_expansionBranchFactors = Arrays.copyOf(_branchingFactors, _branchingFactors.length);
			_branchingFactors = Arrays.copyOf(_branchingFactors, _EXPANSIONSCHEDULE[0]+1);
		}
		_struct = new SparseBackoffTreeStructure(_branchingFactors);
		initDeltas();
	}
	
	public void initDeltas() {
		_wordsDelta = new double[_branchingFactors.length];
		_docsDelta = new double[_branchingFactors.length];
		double initDelta = 0.9 / (double)_branchingFactors.length;
		Arrays.fill(_wordsDelta, initDelta);
		Arrays.fill(_docsDelta, initDelta);
	}
	
	public int get_USEEXPANSION() {
		return _USEEXPANSION;
	}

	public void set_USEEXPANSION(int _USEEXPANSION) {
		this._USEEXPANSION = _USEEXPANSION;
	}

	public int[] get_branchingFactors() {
		return _branchingFactors;
	}


	public void set_branchingFactors(int[] _branchingFactors) {
		this._branchingFactors = _branchingFactors;
	}


	public int[] get_EXPANSIONSCHEDULE() {
		return _EXPANSIONSCHEDULE;
	}

	public void set_EXPANSIONSCHEDULE(int[] _EXPANSIONSCHEDULE) {
		this._EXPANSIONSCHEDULE = _EXPANSIONSCHEDULE;
	}

	public int get_NUMPARTICLES() {
		return _NUMPARTICLES;
	}


	public void set_NUMPARTICLES(int _NUMPARTICLES) {
		this._NUMPARTICLES = _NUMPARTICLES;
	}


	public int get_MAXTESTDOCSIZE() {
		return _MAXTESTDOCSIZE;
	}


	public void set_MAXTESTDOCSIZE(int _MAXTESTDOCSIZE) {
		this._MAXTESTDOCSIZE = _MAXTESTDOCSIZE;
	}
	
	public int get_VOCABSIZE() {
		return _c._VOCABSIZE;
	}

	public void set_VOCABSIZE(int _VOCABSIZE) {
		this._c._VOCABSIZE = _VOCABSIZE;
	}

	public int get_NUMDOCS() {
		return _c._NUMDOCS;
	}

	public void set_NUMDOCS(int _NUMDOCS) {
		this._c._NUMDOCS = _NUMDOCS;
	}

	public long get_NUMTOKENS() {
		return _c._NUMTOKENS;
	}

	public void set_NUMTOKENS(long _NUMTOKENS) {
		this._c._NUMTOKENS = _NUMTOKENS;
	}

	public int get_NUMTHREADS() {
		return _NUMTHREADS;
	}

	public void set_NUMTHREADS(int _NUMTHREADS) {
		this._NUMTHREADS = _NUMTHREADS;
	}

	public int[] get_THREADBREAKS() {
		return _THREADBREAKS;
	}

	public void set_THREADBREAKS(int[] _THREADBREAKS) {
		this._THREADBREAKS = _THREADBREAKS;
	}

	public int get_NUMITERATIONS() {
		return _NUMITERATIONS;
	}

	public void set_NUMITERATIONS(int _NUMITERATIONS) {
		this._NUMITERATIONS = _NUMITERATIONS;
	}

	public int get_OPTIMIZEINTERVAL() {
		return _OPTIMIZEINTERVAL;
	}

	public void set_OPTIMIZEINTERVAL(int _OPTIMIZEINTERVAL) {
		this._OPTIMIZEINTERVAL = _OPTIMIZEINTERVAL;
	}

	private void setThreadBreaks() {
		_THREADBREAKS = new int[_NUMTHREADS];
		long approxToks = _c._NUMTOKENS / (long)_NUMTHREADS;
		long ct = 0;
		int thread = 0;
		for(int i=0; i<_c._NUMDOCS; i++ ) {
			ct += _c._docs[i].size();
			if(ct > approxToks) {
				_THREADBREAKS[thread++] = i;
				ct = 0;
			}
		}
		//extra goes in last thread:
		_THREADBREAKS[_NUMTHREADS - 1] = _c._NUMDOCS - 1;
	}

	public int [] defaultExpansion(int numLevels) {
		int [] out = new int[numLevels];
		for(int i=0; i<out.length; i++) {
			out[i] = i;
		}
		return out;
	}
	
	//config file has lines of form <param-name without underscore>\t<integer value> [<integer value> ...]
	public void readConfigFile(String inFile) throws Exception {
		System.out.println("reading config file.");
		BufferedReader brIn = new BufferedReader(
				new InputStreamReader(new FileInputStream(inFile), "UTF8" ));
		String sLine;
		Method [] ms = this.getClass().getMethods();
		while((sLine = brIn.readLine())!=null) {
			if(sLine.startsWith("#"))
				continue;
			System.out.println(sLine);
			String [] fields = sLine.split("\t");
			//HACK; this only works given the get/set naming convention
			Method m = null;
			for(int i=0; i<ms.length; i++) {
				if(ms[i].getName().equals("set_" + fields[0])) {
					m = ms[i];
					break;
				}
			}
			Class<?> c = m.getParameterTypes()[0];
			if(c.getName().equals("int")) {
				int arg = Integer.parseInt(fields[1]);
				System.out.println(m.getName() + "\t" + arg);
				m.invoke(this, arg);
			}
			else if(c.isArray() && c.getComponentType().getName().equals("int")) {
				String [] vals = fields[1].split(" ");
				int [] arg = new int[vals.length];
				for(int i=0; i<arg.length; i++) {
					arg[i] = Integer.parseInt(vals[i]);
				}
				System.out.println(m.getName() + "\t" + Arrays.toString(arg));
				m.invoke(this, arg);
			}
		}
		brIn.close();
	}
	
	private TIntDoubleHashMap aggregateCounts(TIntArrayList occurs) {
		TIntDoubleHashMap out = new TIntDoubleHashMap();
		TIntIterator it = occurs.iterator();
		while(it.hasNext()) {
			int i = it.next();
			out.adjustOrPutValue(i, 1.0, 1.0);
		}
		return out;
	}
	
	private TIntDoubleHashMap [] aggregateDocCounts() {
		TIntDoubleHashMap [] out = new TIntDoubleHashMap[_c._docs.length];
		for(int i=0; i<_c._docs.length; i++) {
			out[i] = aggregateCounts(_c._z[i]);
		}
		return out;
	}
	
	//scans doc, resampling each variable.  Not used at test time.
	//returns number of changes
	private int sampleDoc(int docId) {
		int changes = 0;
		TIntArrayList zs = _c._z[docId];
		TIntArrayList doc = _c._docs[docId];
		TIntDoubleHashMap docCts = aggregateCounts(zs);
		SparseBackoffTree sbtDoc = new SparseBackoffTree(_struct);
		//System.out.println("adding " + docCts.toString());
		sbtDoc.addAllMass(docCts, this._docsDelta);
		for(int i=0; i<doc.size(); i++) {
			int newZ = sampleZ(docId, i, false, sbtDoc);
			if(newZ != zs.get(i))
				changes++;
			zs.set(i, newZ);
		}
		return changes;
	}
	
	public int sampleZ(int doc, int pos, boolean testing, SparseBackoffTree topicGivenDoc) {
		int w = _c._docs[doc].get(pos);
		int curZ = _c._z[doc].get(pos);
		double bMarginalInv = 1.0f;
		if(curZ >= 0)
			bMarginalInv = 1.0 / _topicMarginal[curZ];
		SparseBackoffTree [] sbtPair = new SparseBackoffTree [] {topicGivenDoc, _topicGivenWord[w]};
		double [] subPair = new double [] {1.0, bMarginalInv};
		if(testing) //don't subtract from word distribution at test time:
			subPair[1] = 0.0;
		SparseBackoffTreeIntersection sbti = new SparseBackoffTreeIntersection(sbtPair, curZ, subPair);
		return sbti.sample(_r);
	}
	
    private long sampleRange(int start, int end) {
    	
    	long chg = 0;
    	
		for(int i=start; i<=end; i++) {
			chg += (long)sampleDoc(i);
		}

		return chg;
    }	
	
	
	public void initZRandom(TIntArrayList [] ds, int maxVal) {
		_c._z = new TIntArrayList[ds.length];
		for(int i=0; i<ds.length; i++) {
			_c._z[i] = new TIntArrayList(ds[i].size());
			for(int j=0; j<ds[i].size(); j++)
				_c._z[i].add(_r.nextInt(maxVal));
		}
	}
	
	public TIntDoubleHashMap [] aggregateWordCounts() {
		TIntDoubleHashMap [] hm = new TIntDoubleHashMap[_c._VOCABSIZE];
		for(int i=0; i<_c._z.length; i++) {
			for(int j=0; j<_c._z[i].size(); j++) {
				int wordID = _c._docs[i].get(j);
				int z = _c._z[i].get(j);
				if(hm[wordID] == null)
					hm[wordID] = new TIntDoubleHashMap();
				hm[wordID].adjustOrPutValue(z, 1.0, 1.0);
			}
		}
		return hm;
	}
	
	
	
	public SparseBackoffTree [] getWordParamsFromZs(double [] dsWords) {
		SparseBackoffTree [] out = new SparseBackoffTree[_c._VOCABSIZE];
		for(int i=0; i<_c._VOCABSIZE; i++) {
			out[i] = new SparseBackoffTree(_struct);
		}
		//aggregate:
		TIntDoubleHashMap [] hm = aggregateWordCounts();
		//add:
		for(int i=0; i<_c._VOCABSIZE; i++) {
			if(hm[i] == null)
				continue;
			TIntDoubleIterator it = hm[i].iterator();
			while(it.hasNext()) {
				it.advance();
				int z = it.key();
				double val = it.value();
				out[i].smoothAndAddMass(z, val, dsWords);
			}
		}
		return out;
	}
	
	public static double [] interpolateEndPoints(double [] endPts, int length) {
		double [] out = new double[length];
		double delt = (endPts[1] - endPts[0])/(double)(length*length - 1);
		for(int i=0; i<length; i++) {
			out[i] = delt*i*i + endPts[0];
		}
		return out;
	}
	
	//reads the corpus and initializes zs and model
	public int initializeForCorpus(String inFile, int maxVal) throws Exception {
		int toks = _c.readCorpusDat(inFile, true);
		setThreadBreaks();
		initZRandom(_c._docs, maxVal);
		updateModel();
		//updateWordParamsFromZs(interpolateEndPoints(this._wordsDeltaEndPts, _branchingFactors.length));
		return toks;
	}
	
    private int gibbsPass() { 
    	int changes= 0 ;
	    GibbsDoer [] gds = new GibbsDoer[_NUMTHREADS];
	    long stTime = System.currentTimeMillis();
    		ExecutorService e = Executors.newFixedThreadPool(_NUMTHREADS);
    		for(int i=0; i<_NUMTHREADS;i++) {
    			gds[i] = new GibbsDoer();
    			gds[i]._start = 0;
    			if(i > 0)
    				gds[i]._start = _THREADBREAKS[i-1] + 1;
    			gds[i]._end = _THREADBREAKS[i];
    			e.execute(gds[i]);
    		}
    		e.shutdown();
    		boolean terminated = false;
    		while(!terminated) {
	    		try {
	    			terminated = e.awaitTermination(60,  TimeUnit.SECONDS);
	    		}
	    		catch (InterruptedException ie) {
	    			
	    		}
    		}
    		for(int i=0; i<_NUMTHREADS; i++) {
    			changes += gds[i]._changes;
    		}
    	stTime = System.currentTimeMillis() - stTime;
    	System.out.println("\ttime: " + stTime + "\tchanges: " + changes);
    	//System.out.println(Arrays.toString(_topicMarginal));
    	return changes;
    }
	
    
    
	//returns array of amt_i
	//such that if we divide leaf counts by amt_i, we get smoothing with marginal P(z)
    //(see paper)
	public static double [] getNormalizers(SparseBackoffTree [] shds, SparseBackoffTreeStructure struct) {
		int numStates = struct.numLeaves();
		double [] count = new double[numStates];
		double [] smoothing = new double[numStates];
		SparseBackoffTree shdAgg = SparseBackoffTree.sum(shds, struct);
		double maxSmooth = 0.0f;
		double sumCount = 0.0f;
		double sumSmoother = 0.0f;
		for(int i=0; i<numStates; i++) {
			double [] smoothAndCount = shdAgg.getSmoothAndCount(i);
			smoothing[i] = smoothAndCount[0];
			if(smoothing[i] > maxSmooth) {
				maxSmooth = smoothing[i];
			}
			count[i] = smoothAndCount[1];
			sumCount += count[i];
			sumSmoother += smoothing[i];
		}
		double [] out = new double[numStates];
		double target = Math.max(maxSmooth + 1.0f, (sumSmoother + sumCount)/numStates);
		for(int i=0; i<out.length; i++) {
			out[i] = count[i]/(target + 0.001f - smoothing[i]);
			if(out[i] < 0.0f)
				System.out.println("zero or negative normalizer!");
		}
		return out;
	}
	
	
	
	public static double [] getCheckMarginal(SparseBackoffTree [] shds, SparseBackoffTreeStructure struct) {
		int numStates = struct.numLeaves();
		double [] out = new double[numStates];
		SparseBackoffTree shdAgg = SparseBackoffTree.sum(shds, struct);
		for(int i=0; i<numStates; i++) {
			double [] smoothAndCount = shdAgg.getSmoothAndCount(i);
			out[i] += smoothAndCount[0] + smoothAndCount[1];
		}
		return out;
	}
	
	/**
	 * Updates the model given the topic assignments (_z) and divides by marginal for next sampling pass
	 */
    public void updateModel() {
		System.out.println("first doc samples: " + _c._z[0].toString());
    	_topicGivenWord = getWordParamsFromZs(_wordsDelta);
		_topicMarginal = getNormalizers(_topicGivenWord, _struct);
		for(int i=0; i<_topicGivenWord.length; i++) 
			_topicGivenWord[i].divideCountsBy(_topicMarginal);
		System.out.println("\tdiscounts doc: " + Arrays.toString(_docsDelta) + "\tword: " + Arrays.toString(_wordsDelta));
    }
    
    /**
     * Adds the gradient of the log likelihood for the given observations
     * Assumes sum of deltas < 1
     * @param grad	the gradient is added here
     * @param hm	the observations
     * @param curDeltas	the current hyperparameters
     * @param incremental	if true, measure log likelihood when introducing observations one-by-one in random order;
     * otherwise, use leave-one-out cross validation
     * @return	log likelihood using current parameters
     */
    public double updateGradient(double [] grad, TIntDoubleHashMap hm, double [] curDeltas, boolean incremental) {
    	SparseBackoffTree sbt = new SparseBackoffTree(_struct);
    	double [] leafCount = new double[curDeltas.length];
    	leafCount[leafCount.length - 1] = _branchingFactors[leafCount.length - 1];
    	for(int i=leafCount.length - 2; i>=0;i--) {
    		leafCount[i] = leafCount[i+1]*_branchingFactors[i];
    	}
    	if(hm.size()==1 && hm.containsValue(1.0)) //can't do LOOCV with one observation
    		return 0.0;

    	ArrayList<Integer> zOrder = new ArrayList<Integer>();
    	if(!incremental) {
    		sbt.addAllMass(hm, curDeltas);
    	}
    	else {
    		TIntDoubleIterator it = hm.iterator();
    		while(it.hasNext()) {
        		it.advance();
        		int z = it.key();
        		double v = it.value();
        		//TODO: why is v a double?
    			for(int i=0; i<v; i++) {
    				zOrder.add(z);
    			}
    		}
    		Collections.shuffle(zOrder);
    	}
    	double out = 0.0;
    	double sumDelta = 0.0;
    	double singletonSmoothAdd = 0.0;
    	if(incremental)
    		sbt.smoothAndAddMass(zOrder.get(0), 1.0, curDeltas);
    	
    	for(int i=0; i<curDeltas.length; i++) {
    		sumDelta += curDeltas[i];
    		singletonSmoothAdd += curDeltas[i] / leafCount[i];
    	}

		if(!incremental) {
			TIntDoubleIterator it = hm.iterator();
	    	while(it.hasNext()) {
	    		it.advance();
	    		int z = it.key();
	    		double v = it.value();    			

	    		//TODO: optimize around these expensive steps:
	    		int [] localIdxTrace = _struct.getLocalIdxTrace(z);
	    		double [] smooths = sbt.getSmoothsTrace(localIdxTrace);
	    		
	    		double smoothed = sbt.getSmoothed(z);
	    		//gradient of log likelihood w.r.t. delta_i is:
	    		//if num_obs > 1:
	    		//  num_obs * (smooth_i/delta_i - 1) / [(get_smoothed_i - 1)]
	    		//else
	    		//  (num_leaves * smooth_i/delta - 1) / (num_leaves * (get_smoothed_i + sum_deltas - ))
	    		for(int i=0; i<grad.length; i++) {
		    		if(v > 1.0) {
		    			grad[i] += v * (smooths[i]/curDeltas[i] - 1.0) / (smoothed - 1.0);
		    			out += v * Math.log((smoothed - 1.0) / (sbt._totalMass - 1.0));
		    		}
		    		else {
		    			grad[i] += (leafCount[i]*smooths[i]/curDeltas[i] - 1.0) / (leafCount[i] * (smoothed + sumDelta - 1.0 - singletonSmoothAdd));
		    			out += Math.log((smoothed + sumDelta - 1.0 - singletonSmoothAdd) / (sbt._totalMass - 1.0));
		    		}
	    		}
	    	}
		}
		else {
			Iterator<Integer> it = zOrder.iterator();
			//skip the first one:
			it.next();
			while(it.hasNext()) {
	    		int z = it.next();    			

	    		//TODO: optimize around these expensive steps:
	    		int [] localIdxTrace = _struct.getLocalIdxTrace(z);
	    		double [] smooths = sbt.getSmoothsTrace(localIdxTrace);
	    		
	    		double [] countSmoothed = sbt.getSmoothAndCount(z);
	    		double smoothed = countSmoothed[0] + countSmoothed[1];
	    		//gradient of log likelihood w.r.t. delta_i is:
	    		//if num_obs > 1:
	    		//  num_obs * (smooth_i/delta_i - 1) / [(get_smoothed_i - 1)]
	    		//else
	    		//  (num_leaves * smooth_i/delta - 1) / (num_leaves * (get_smoothed_i + sum_deltas - ))
	    		for(int i=0; i<grad.length; i++) {
		    		if(countSmoothed[1] > 0.0) { //delta has negative as well as positive influence
		    			grad[i] += (smooths[i]/curDeltas[i] - 1.0) / smoothed;
		    		}
		    		else {
		    			grad[i] += (smooths[i]/curDeltas[i]) / smoothed;
		    		}
		    		out += Math.log(smoothed / sbt._totalMass);
	    		}
	    		if(countSmoothed[1] > 0.0)
	    			sbt.addMass(z, 1.0);
	    		else
	    			sbt.smoothAndAddMass(z, 1.0, curDeltas);
	    	}
		}
    	return out;
    }
    
    /**
     * returns gradient normalized to sum to stepSize 
     * (or less, if needed to ensure gradient + curDeltas sums to < 1.0 and has no negative elements)
     * @param gradient	gradient to normalize
     * @param curDeltas	current hyperparameters to which gradient will be applied
     * @param stepSize	how far to move (in L1)
     * @return
     */
    public double [] normalizeAndCheckBounds(double [] gradient, double [] curDeltas, double stepSize) {

    	double normalizer = 1.0 / stepSize;
    	double [] out = new double[gradient.length];

    	while(true) {
	    	double sum = 0.0;
	    	for(int i=0; i<gradient.length; i++) {
	    		sum += Math.abs(gradient[i]);
	    	}
	    	sum *= normalizer;
	    	//sum = normalizer;
	    	
	    	boolean outOfBounds = false;
	    	double newSum = 0.0;
	    	for(int i=0; i<gradient.length; i++) {
	    		out[i] = gradient[i]/sum;
	    		if(out[i] + curDeltas[i] < 1E-5) {
	    			out[i] = 1E-5 - curDeltas[i];
	    		}
	    		newSum += out[i] + curDeltas[i];
	    	}
	    	if(newSum >= 1.0)
	    		normalizer *= 2.0;
	    	else {
	    		if(newSum < stepSize) {
	    			for(int i=0; i<gradient.length; i++) {
	    				out[i] *= (stepSize / newSum);
	    			}
	    		}
	    		break;
	    	}
    	}
    	
    	return out;
    }
    
    public void addAtoB(double [] a, double [] b) {
    	for(int i=0; i<a.length; i++) 
    		b[i] += a[i];
    }
    
    /**
     * Applies and returns gradient
     * @param deltas	hyperparameters to start from
     * @param hm	observation counts
     * @param stepSize	how far to move in the step, in L1 norm
     * @param incremental	if true, measure log likelihood when introducing observations one-by-one in random order;
     * otherwise, use leave-one-out cross validation
     * @return
     */
    public double [] gradientStep(double [] deltas, TIntDoubleHashMap [] hm, double stepSize, boolean incremental) {
    	double [] gradient = new double[deltas.length];
    	double LL = 0.0;
    	for(int i=0; i<hm.length; i++) {
    		if(hm[i] != null)
    			LL += updateGradient(gradient, hm[i], deltas, incremental);
    	}
    	System.out.println("LL: " + LL);
    	System.out.println("got grad " + Arrays.toString(gradient));
    	double [] newGradient = this.normalizeAndCheckBounds(gradient, deltas, stepSize);
    	System.out.println("applying " + Arrays.toString(newGradient));
    	addAtoB(newGradient, deltas);
    	return gradient;
    }
    
    /**
     * Optimizes parameters using gradient ascent in log likelihood
     */
    public void optimizeParameters() {
//    	for(int i=0; i<_wordsDelta.length; i++) {
//    		_wordsDelta[i] *= 0.9;
//    		_docsDelta[i] *= 0.9;
//    	}
    	//TODO: multi-thread
//    	double STEPSIZE = 0.01; //start stepping this far in L1 norm
//    	double STEPDEC = 0.95; //decrease step size this much each step
//    	int NUMSTEPS = 20; //number of steps to take
    	double STEPSIZE = 0.02; //start stepping this far in L1 norm
    	double STEPDEC = 0.8; //decrease step size this much each step
    	int NUMSTEPS = 10; //number of steps to take

    	TIntDoubleHashMap [] hm = aggregateWordCounts();
    	System.out.println("words:");
    	double step = STEPSIZE;
    	for(int i=0; i<NUMSTEPS; i++) {
        	gradientStep(_wordsDelta, hm, step, false);
        	step *= STEPDEC;
    	}
    	hm = aggregateDocCounts();
    	System.out.println("docs:");
    	step = STEPSIZE;
    	for(int i=0; i<NUMSTEPS; i++) {
        	gradientStep(_docsDelta, hm, step, true);
        	step *= STEPDEC;
    	}
    	System.out.println("full word params: " + Arrays.toString(_wordsDelta));
    	System.out.println("full doc params: " + Arrays.toString(_docsDelta));
    }
    
    /**
     * Expands the model
     * @param expansionIndex	The desired index in _EXPANSIONSCHEDULE to expand to.
     * @return	whether the model was expanded
     */
    public boolean expandModel(int expansionIndex) {
    	if(expansionIndex < _EXPANSIONSCHEDULE.length) {
    		int curLeaves = _struct.numLeaves(); 
    		//get new structure
    		_branchingFactors = Arrays.copyOf(this._expansionBranchFactors, _EXPANSIONSCHEDULE[expansionIndex] + 1);
    		_struct = new SparseBackoffTreeStructure(_branchingFactors);
    		initDeltas();
    		int multiplier = _struct.numLeaves()/curLeaves;
    		System.out.println("Expanding to " + Arrays.toString(_branchingFactors));
    		for(int i=0; i<_c._z.length;i++) {
    			for(int j=0; j<_c._z[i].size(); j++) {
    				int z = multiplier * _c._z[i].get(j) + _r.nextInt(multiplier);
    				_c._z[i].set(j,  z);
    			}
    		}
    		updateModel();
    		return true;
    	}
    	else
    		return false;
    }
    
    /**
     * Trains the model for a specified number of iterations.
     *
     * @param  iterations  the number of Gibbs passes to perform
     * @param  updateInterval	number of iterations between hyperparameter updates
     */
	public void trainModel(int iterations, int updateInterval) {
		boolean done = false;
		int j=0;
		while(!done) {
			for(int i=1; i<=iterations; i++) {
				System.out.print(i + ": ");
				gibbsPass();
				if((i % updateInterval)==0) { 
					optimizeParameters();
				}
				updateModel();
			}
			if(this._USEEXPANSION == 1) {
				j++;
				done = !expandModel(j);
			}
			else 
				done = true;
		}		
	}
	
	public void writeModel(String outFile) throws Exception {
		SparseBackoffTree [] sbt = this._topicGivenWord;
		SparseBackoffTreeStructure struct = this._struct;
		_topicGivenWord = null;
		_struct = null;
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outFile));
		oos.writeObject(this);
		oos.close();
		_topicGivenWord = sbt;
		_struct = struct;
	}
	
	public static SBTDocumentModel readModel(String inFile) throws Exception {
		ObjectInputStream ois = new ObjectInputStream(new FileInputStream(inFile));
		SBTDocumentModel out = (SBTDocumentModel) ois.readObject();
		ois.close();
		out._struct = new SparseBackoffTreeStructure(out._branchingFactors);
		out.updateModel();
		//out.updateWordParamsFromZs(interpolateEndPoints(out._wordsDeltaEndPts, out._branchingFactors.length));
		return out;
	}
	
	
	//tests on a single document using left-to-right method
	//word topic marginal (i.e. P(topic) computed from P(topic, word)) is supplied to ensure evaluated distribution sums to one
	public double [] testOnDocFullPpl(int doc, double [] ds, double [] wordTopicMarginal) {
		boolean CHECKSUMTOONE = false; //slow.  Can use to confirm that test code uses a valid probability distribution over words
		double LL = 0.0;
		double numWords = 0.0;
		ArrayList<Integer> wordOrder = new ArrayList<Integer>();
		for(int i=0; i<_c._docs[doc].size(); i++) {
			wordOrder.add(i);
		}
		Collections.shuffle(wordOrder);
		double [] ps = new double[wordOrder.size()];
		for(int j=0; j<_NUMPARTICLES; j++) {
			SparseBackoffTree topicGivenDoc = new SparseBackoffTree(_struct);
			numWords = 0.0;
			for(int i=0; i<_c._docs[doc].size(); i++) {
				_c._z[doc].set(i, -1);
			}
			for(int i=0; i<_c._docs[doc].size(); i++) {
				int widx = wordOrder.get(i);
				int w = _c._docs[doc].get(widx);
				double p = 0.0;
				//only compute for seen words:
				if(_c._pWord[w] > 0) {
					//init topicGivenDoc with all samples:
					topicGivenDoc = new SparseBackoffTree(_struct);
					for(int k =0; k<_c._z[doc].size(); k++) {
						int r = _c._z[doc].get(k);
						if(r != -1)
							topicGivenDoc.addAndSmoothIfZero(r, 1.0, ds);
					}
						//resample everything before:
					for(int k=0; k<i; k++) {
						int wkidx = wordOrder.get(k);
						if(_c._z[doc].get(wkidx) != -1) { //seen word
							int r = sampleZ(doc, wkidx, true, topicGivenDoc);
							_c._z[doc].set(wkidx, r);
						}
					}

					//test on this word:
					double [] pWordGivenState = new double[wordTopicMarginal.length];
					double [] pStateGivenDoc = new double[wordTopicMarginal.length];
					float checkSum = 0.0f;
					for(int t=0;t<wordTopicMarginal.length; t++) {
						
//						if(wordTopicMarginal[t]==0.0) { //added and commented 9/11/2014
//							continue;
//						}
						//dividing by wordTopicMarginal ensures we use a distribution P(w | t) that sums to one
						pWordGivenState[t] = _topicGivenWord[w].getSmoothed(t) / wordTopicMarginal[t];
						double numt = topicGivenDoc.getSmoothed(t);
//						if(numt==0.0) {
//							continue;
//						}
						if(numWords==0) {
							pStateGivenDoc[t] = 1.0f / wordTopicMarginal.length;
						}
						else {
							pStateGivenDoc[t] = numt / topicGivenDoc._totalMass;
						}
						checkSum += numt;
						p += (pWordGivenState[t]*pStateGivenDoc[t]);
						if(Double.isNaN(p))
							System.err.println("nan.");
						if(p < 0.0)
							System.err.println("negative.");
					}
					if((i==0 || (i== _c._docs[doc].size()-1)) && (CHECKSUMTOONE)) {
						float sDoc = 0.0f;
						float sWord = 0.0f;
						float totalTopicMarginal = 0.0f;
						for(int t=0; t<pWordGivenState.length; t++) {
							sDoc += pStateGivenDoc[t];
							sWord += pWordGivenState[t]*wordTopicMarginal[t];
							totalTopicMarginal += wordTopicMarginal[t];
						}
						sWord /= totalTopicMarginal;
						System.out.println("Check SUM: " + i + "\t" + p + "\t" + sDoc + "\t" + sWord + "\t" + _c._pWord[w]);
					}
					ps[widx] += p;
					//sample this word:
					int r = -1;
					if(numWords>0)
						r = sampleZ(doc, widx, true, topicGivenDoc);
					else
						r = _topicGivenWord[w].sample(_r);
					_c._z[doc].set(widx, r);
					//TODO: make this less wildly inefficient
//					topicGivenDoc = new SparseBackoffTree(_struct);
//					for(int k =0; k<_z[doc].size(); k++) {
//						r = _z[doc].get(k);
//						if(r != -1)
//							topicGivenDoc.smoothAndAddMass(r, 1.0, ds);
//					}
					numWords++;
				}
			}
		}
		numWords = 0.0;
		for(int i=0; i<_c._docs[doc].size(); i++) {
			int widx = wordOrder.get(i);
			int w = _c._docs[doc].get(widx);
			if(_c._pWord[w] > 0) {
				LL += Math.log(ps[widx]/(double)_NUMPARTICLES);
				if(Double.isNaN(LL))
					System.out.println("got NaN with " + ps[widx]);
				numWords++;
			}
		}
		return new double [] {LL, numWords};
	}

	
	public static void train(String inputFile, String outputFile, String configFile) throws Exception {
		SBTDocumentModel sbtdm = new SBTDocumentModel(configFile);
		sbtdm.initializeForCorpus(inputFile, sbtdm._struct.numLeaves());
		sbtdm.trainModel(sbtdm._NUMITERATIONS, sbtdm._OPTIMIZEINTERVAL);
		sbtdm.writeModel(outputFile);
	}
	
	//computes log likelihood of model using left-to-right method
	//returns {ppl, number of tested words}
	//numDocs is number in test file
	//maxDocs is number actually tested (starting from the beginning of the file)
	public static double [] testModel(String modelFile, String testFile, int numDocs, int maxDocs, String configFile) throws Exception {
		boolean CHECKWORDDISTRIBUTION = false;
		
		SBTDocumentModel sbtdm = readModel(modelFile);
		sbtdm.readConfigFile(configFile);
		double [] wordTopicNormalizer = getCheckMarginal(sbtdm._topicGivenWord, sbtdm._struct); 
		int numStates = wordTopicNormalizer.length;

		float [] pWordGivenState = new float[wordTopicNormalizer.length] ; 
		if(CHECKWORDDISTRIBUTION) {
			for(int w=0; w<sbtdm._topicGivenWord.length; w++) {
				for(int t=0;t<wordTopicNormalizer.length; t++) {
					pWordGivenState[t] += sbtdm._topicGivenWord[w].getSmoothed(t) / wordTopicNormalizer[t];
				}
			}
			float sum = 0.0f;
			for(int i=0; i<pWordGivenState.length; i++) {
				sum += pWordGivenState[i];
			}
			
			System.out.println("Word sum: " + sum + " should be about " + numStates);
		}
		
		sbtdm._c.reInitForCorpus(testFile, numDocs);
		sbtdm.setThreadBreaks();
		double LL = 0.0;
		double numWords = 0.0;
		TestDoer [] tds = new TestDoer[maxDocs];
		ExecutorService e = Executors.newFixedThreadPool(sbtdm._NUMTHREADS);
		//ExecutorService e = Executors.newFixedThreadPool(1);
		for(int i=0; i<tds.length;i++) {
			if(sbtdm._MAXTESTDOCSIZE < 0 || sbtdm._c._docs[i].size() > sbtdm._MAXTESTDOCSIZE) {
				System.out.println("skipping " + i + " size " + sbtdm._c._docs[i].size());
				continue;
			}
			tds[i] = sbtdm.new TestDoer(i, sbtdm._docsDelta, wordTopicNormalizer);
			e.execute(tds[i]);
		}
		e.shutdown();
		boolean terminated = false;
		while(!terminated) {
    		try {
    			terminated = e.awaitTermination(60,  TimeUnit.SECONDS);
    		}
    		catch (InterruptedException ie) {
    			
    		}
		}
		for(int i=0; i<tds.length; i++) {
			if(tds[i]==null)
				continue;
			LL += tds[i]._res[0];
			numWords += tds[i]._res[1];
		}
		return new double [] {LL, numWords};
	}
	
	public static void test(String modelFile, String inputFile, int numDocsInFile, int numDocsToTest, String configFile) throws Exception {
		double [] ll = testModel(modelFile, inputFile, numDocsInFile, numDocsToTest, configFile);
		System.out.println("Test ppl: " + Arrays.toString(ll));
	}
	
	/**
	 * Outputs the sparse word-to-topic vector proportional to P(t | w)/P(t) for words and topics with positive counts
	 * @param modelFile	Contains the model
	 * @param outputFile	Where to write the output
	 * @throws Exception
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public static void outputWordToTopicFile(String modelFile, String outputFile) throws Exception {
		BufferedWriter bwOut = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), "UTF8"));
		SBTDocumentModel sbtdm = readModel(modelFile);
//		for(int i=0; i<sbtdm._c._pWord.length; i++) {
//			if(sbtdm._c._pWord[i] > 0.0) {
//				bwOut.write(i + "\t");
//				TIntDoubleHashMap hm = sbtdm._topicGivenWord[i].getLeafCounts();
//				TIntDoubleIterator it = hm.iterator();
//				while(it.hasNext()) {
//					it.advance();
//					int wId = it.key();
//					double val = it.value();
//					bwOut.write(wId + ":" + val + " ");
//				}
//				bwOut.write("\r\n");
//			}
//		}
		SparseBackoffTreeStructure root = sbtdm._struct;
		int topSum = 0;
		for(int leaf : root._numLeaves){
			topSum += leaf;
		}

		HashMap<Integer,TreeMap<Integer,Double>> topic =  word_distribute(modelFile);
		
		
		
		
		helper(bwOut, sbtdm, root, 0, topic);

		bwOut.close();
	}
	// this function writes json
	public static void helper(BufferedWriter bwOut, SBTDocumentModel sbtdm, SparseBackoffTreeStructure curr, int indent, HashMap<Integer,TreeMap<Integer,Double>> topic) throws Exception {
		String ind = "";
		for (int i = 0; i < indent; i++){
			ind += "\t";
		}

		bwOut.write(ind + "{\"name\": " + "\""+ curr._delta + "\"," + "\n");
		if(curr._children != null){
			bwOut.write(ind + "\"children\": " + "[\n");
			int lenOfChildren = curr._children.length;
			for(int i = 0; i < lenOfChildren - 1; i++){
				helper(bwOut, sbtdm, curr._children[i], indent + 1, topic);
				bwOut.write(",\n");
			}
			helper(bwOut, sbtdm, curr._children[lenOfChildren - 1], indent + 1, topic);
			bwOut.write("\n" + ind + "\t]");
		}
		// if we reach the second last layer
		else {
			bwOut.write(ind + "\"children\": " + "[\n");
			int len = curr._numLeaves.length;
			// output the topics
			int per = 10;
			for(int j = 0; j < len - 1 ; j++){
				bwOut.write(ind + "\t{\"name\": " + "\""+ curr._delta + "\"");
				bwOut.write(",\n" + ind + "\t\"children\": " + "[\n");
				
				// output words
//				 String key : team1.keySet()
				int topNum = j + curr._minGlobalIndex;
				int idx = 0;
				for(int k : topic.get(topNum).keySet()){
					if(idx < per - 1){
						bwOut.write(ind + "\t\t{\"name\": " + k + "},\n");
					}
					else if(idx == per - 1){
						// deal with last one
						bwOut.write(ind + "\t\t{\"name\": " + k + "}\n");
					}
					else{
						break;
					}
					idx ++;
				}
				
				bwOut.write("\n" + ind + "\t]");
				
				bwOut.write("\n" + ind + "\t},\n");
			}
			
			// deal with last topic
			bwOut.write(ind + "\t{\"name\": " + "\""+ curr._delta + "\"");
			bwOut.write(",\n" + ind + "\t\"children\": " + "[\n");
			int topNum = len - 1 + curr._minGlobalIndex;
			int idx = 0;
			for(int k : topic.get(topNum).keySet()){

				if(idx < per - 1){
					bwOut.write(ind + "\t\t{\"name\": " + k + "},\n");
				}
				else if(idx == per - 1){
					// deal with last one
					bwOut.write(ind + "\t\t{\"name\": " + k + "}\n");
				}
				else{
					break;
				}
				idx ++;
			}
			
			bwOut.write("\n" + ind + "\t]");
			
			bwOut.write("\n" + ind + "\t}\n");
			
			bwOut.write("\n" + ind + "\t]");
		}
		bwOut.write("\n");
		bwOut.write(ind + "}");
	}
	
	public static  HashMap<Integer,TreeMap<Integer,Double>> word_distribute(String modelFile) throws Exception{
	// Set up a hash table of the word distribution given topic
		HashMap<Integer,HashMap<Integer,Double>> word_givenTopic = new HashMap();
		SBTDocumentModel sbtdm = readModel(modelFile);
		for(int i = 0; i < sbtdm._c._pWord.length; i++) {
			if(sbtdm._c._pWord[i] > 0.0) {
				TIntDoubleHashMap hm = sbtdm._topicGivenWord[i].getLeafCounts();
				TIntDoubleIterator it = hm.iterator();
				while(it.hasNext()) {
					it.advance();
					int wId = it.key();
					double val = it.value();
					if(word_givenTopic.containsKey(wId)){
                        HashMap<Integer,Double> sublist = word_givenTopic.get(wId);
                        sublist.put(i,val);
					}else{
						HashMap<Integer,Double> sublist = new HashMap();
						sublist.put(i,val);
						word_givenTopic.put(wId,sublist);
					}
				}
			}
		}

	//sort the hash map by value.
	@SuppressWarnings("rawtypes")
	HashMap<Integer,TreeMap<Integer,Double>> final_result = new HashMap();
	for(Entry<Integer, HashMap<Integer, Double>> entry : word_givenTopic.entrySet()) {
			int key = entry.getKey();
    		HashMap<Integer,Double> word_distribution = entry.getValue();
    		ByValueComparator bvc = new ByValueComparator(word_distribution);
    		TreeMap<Integer, Double> sorted_word_distribution = new TreeMap<Integer, Double>(bvc);
            sorted_word_distribution.putAll(word_distribution);
            final_result.put(key,sorted_word_distribution);
    }    
    return final_result;
}
static class ByValueComparator implements Comparator<Integer> {
    HashMap<Integer, Double> base_map;
    
    public ByValueComparator(HashMap<Integer, Double>  base_map) {
        this.base_map = base_map;
    }

    public int compare(Integer arg0, Integer arg1) {
        if (!base_map.containsKey(arg0) || !base_map.containsKey(arg1)) {
            return 0;
        }

        if (base_map.get(arg0) < base_map.get(arg1)) {
            return 1;
        } else if (base_map.get(arg0) == base_map.get(arg1)) {
            return 0;
        } else {
            return -1;
        }
    }
}
	public static void main(String[] args) throws Exception {
		if(args.length > 0 && args[0].equalsIgnoreCase("train")) {
			if(args.length != 4) {
				System.err.println("Usage: train <input_file> <model_output_file> <configuration_file>");
				return;
			}
			train(args[1], args[2], args[3]);
		} 
		else if(args.length > 0 && args[0].equalsIgnoreCase("test")) {
			if(args.length != 6) {
				System.err.println("Usage: test <model_file> <test_file> <configuration_file> <num_docs_in_test_file> <num_docs_to_test>");
			}
			test(args[1], args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3]);
		}
		else if(args.length > 0 && args[0].equalsIgnoreCase("wordreps")) {
			if(args.length != 3) {
				System.err.println("Usage: wordreps <model_file> <output_file>");
			}
			outputWordToTopicFile(args[1], args[2]);
		}
		else {
			System.err.println("Usage: <train|test|wordreps> <args...>");
		}
	}

}
