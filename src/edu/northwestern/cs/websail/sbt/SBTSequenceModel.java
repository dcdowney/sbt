package edu.northwestern.cs.websail.sbt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntDoubleHashMap;

//TODO: factor out all the gross duplication between this class and SBTDocumentModel

public class SBTSequenceModel implements Serializable {

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
    	double [] _wordTopicMarginal;
    	double [] _res;
    	
    	public TestDoer(int i, double [] wordTopicMarginal) {
    		_doc = i;
    		_wordTopicMarginal = wordTopicMarginal;
    	}
    	
    	public void run() {
    		_res = testOnDocFullPpl(_doc, _wordTopicMarginal);
    		//_res = testOnDocFullPplExact(_doc, _wordTopicMarginal);
    		System.out.println(_doc + "\t" + _res[0] + "\t" + _res[1]);
    	}
    }

    public static class TestGetProbsDoer extends TestEnsembleExactDoer {
    	double [][] _probs;
    	
    	public TestGetProbsDoer(int i, double [][] wordTopicMarginal, SBTSequenceModel [] models, double [][][] tm, double [] ws) {
    		super(i, wordTopicMarginal, models, tm, ws);
    	}
    	
    	@Override
    	public void run() {
    		_probs = getEnsembleProbsForDoc(_doc, _wordTopicMarginal, _sbtsm, _tm);
    	}
    }
    
    private class TestExactDoer implements Runnable {
     	int _doc;
    	double [] _wordTopicMarginal;
    	double [] _res;
    	double [][] _tm;
    	
    	public TestExactDoer(int i, double [] wordTopicMarginal, double [][] transModel) {
    		_doc = i;
    		_wordTopicMarginal = wordTopicMarginal;
    		_tm = transModel;
    	}
    	
    	public void run() {
    		//_res = testOnDocFullPpl(_doc, _wordTopicMarginal);
    		_res = testOnDocFullPplExact(_doc, _wordTopicMarginal, _tm);
    		System.out.println(_doc + "\t" + _res[0] + "\t" + _res[1]);
    	} 	
    }
    
    public static class TestEnsembleDoer implements Runnable {
    	int _doc;
    	double [][] _wordTopicMarginal;
    	SBTSequenceModel [] _sbtsm;
    	double [] _res;
    	
    	public TestEnsembleDoer(int i, double [][] wordTopicMarginal, SBTSequenceModel [] models) {
    		_doc = i;
    		_wordTopicMarginal = wordTopicMarginal;
    		_sbtsm = models;
    	}
    	
    	public void run() {
    		_res = testEnsembleOnDocFullPpl(_doc, _wordTopicMarginal, _sbtsm);
    		System.out.println(_doc + "\t" + _res[0] + "\t" + _res[1]);
    	}
    }
	
    public static class TestEnsembleExactDoer implements Runnable {
    	int _doc;
    	double [][] _wordTopicMarginal;
    	double [][][] _tm;
    	SBTSequenceModel [] _sbtsm;
    	double [] _res;
    	double [] _ws;
    	double [] _ps;
    	
    	public TestEnsembleExactDoer(int i, double [][] wordTopicMarginal, SBTSequenceModel [] models, double [][][] tm,
    			double [] ws) {
    		_doc = i;
    		_wordTopicMarginal = wordTopicMarginal;
    		_sbtsm = models;
    		_tm = tm;
    		_ws = ws;
    	}
    	
    	public void run() {
    		_ps = new double[_sbtsm[0]._c._docs[_doc].size()];
    		_res = testEnsembleOnDocExact(_doc, _wordTopicMarginal, _sbtsm, _tm, _ws, _ps); //side effect, fills _ps with probs per tok
    		System.out.println(_doc + "\t" + _res[0] + "\t" + _res[1]);
    	}
    }
    
	//Model information:
    //state index zero is special start-of-sentence state
	int [] _branchingFactors;
	SparseBackoffTree [] _forward;
	SparseBackoffTree _startState;
	SparseBackoffTree _endState;
	SparseBackoffTree [] _backward;
	SparseBackoffTree [] _wordToState;
	double [] _wordStateMarginal;
	double [] _backwardStateMarginal;
	
	SparseBackoffTreeStructure _struct;
	double [] _forwardDelta = null;
	double [] _backwardDelta = null;
	double [] _wordDelta = null;
	
	Corpus _c;
	
	double [][] baseProbs = null; //holds existing model probabilities when training a layer of an ensemble
	
	Random _r = new Random(4);
	
	//TODO: think about wrapping all of threading/train config up into a trainer class
	
	//threading:
	private int _NUMTHREADS = -1;
	private int[] _THREADBREAKS = null; //inclusive doc indices where threads *end* (initial thread starts
							//at 0 implicitly)
	
	//training config:
	protected int _NUMITERATIONS = -1;
	protected int _OPTIMIZEINTERVAL = -1;
	private int [] _EXPANSIONSCHEDULE = null; //specifies indexes of branching factor to introduce incrementally
	private int _USEEXPANSION = 1; //1 for yes, 0 for no
	private int [] _expansionBranchFactors = null; //when expansion is used, holds the full branching factors
	private final boolean _SAMPLEONSAMPLES = false;
	
	//testing config:
	protected int _NUMPARTICLES = -1;
	private int _MAXTESTDOCSIZE = -1;
	
	private int _numSkipWords = 0;
	
	public SBTSequenceModel(String configFile) throws Exception {
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
		_forwardDelta = new double[_branchingFactors.length];
		_backwardDelta = new double[_branchingFactors.length];
		_wordDelta =new double[_branchingFactors.length];
		double initDelta = 0.9 / (double)_branchingFactors.length;
		Arrays.fill(_wordDelta, initDelta);
		Arrays.fill(_forwardDelta, initDelta);
		Arrays.fill(_backwardDelta, initDelta);
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
		
		long approxToks = _c._NUMTOKENS / _NUMTHREADS;
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

	private double [][] readBaseProbs(String inFile) throws Exception {
		double [][] out = new double[_c._docs.length][];
		BufferedReader brIn = new BufferedReader(
				new InputStreamReader(new FileInputStream(inFile), "UTF8" ));
		String sLine;

		for(int i=0; i<out.length; i++) {
			sLine = brIn.readLine();
			String [] sps = sLine.split(" ");
			double [] ps = new double[_c._docs[i].size()];
			for(int j=0; j<ps.length; j++) {
				ps[j] = Double.parseDouble(sps[j]);
			}
			out[i] = ps;
		}
		brIn.close();
		return out;
	}
	
	private static void writeBaseProbs(String outFile, double [][] ps) throws Exception {
		BufferedWriter bwOut = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "UTF8"));
		for(int i=0; i<ps.length; i++) {
			for(int j=0; j<ps[i].length;j++) {
				if(j>0)
					bwOut.write(" ");
				bwOut.write(ps[i][j] + "");
			}
			bwOut.write("\r\n");
		}
		bwOut.close();
	}
	
	private double [][] transModel() {
		int numStates = this._wordStateMarginal.length;
		double [][] out = new double[numStates][numStates];
		double diff = 0.0;
		double total = 0.0;
		double [] nextStateMarginal = new double[numStates];
		double [] invBackwardMarginal = new double[numStates];
		for(int i=0; i<invBackwardMarginal.length; i++) {
			invBackwardMarginal[i] = 1.0 / _backwardStateMarginal[i];
		}
		for(int i=0; i<_backward.length; i++) {
			_backward[i].divideCountsBy(invBackwardMarginal);
		}
		//first compute estimates based on forward only:
		for(int i=0; i<numStates; i++) {
			for(int j=0; j<numStates;j++) {
				double forwardEst = 0.0;
				double countMass = _forward[i].getSmoothed(j); 
				if(_forward[i]._totalMass > 0.0)
					forwardEst = countMass /_forward[i]._totalMass;
				else
					forwardEst = 1.0/(double)numStates;
				out[i][j] = forwardEst;
				nextStateMarginal[j] += countMass;
				total += countMass;
			}
		}
		total /= (double)numStates;
		for(int i=0; i<numStates; i++) {
			nextStateMarginal[i] /= total;
		}
		
		double [][] be = new double [numStates][numStates];
		for(int i=0; i<numStates; i++) {
			double jSum = 0.0;
			for(int j=0; j<numStates;j++) {
				double backwardEst = 0.0;
				if(_backward[j]._totalMass > 0.0)
					//backwardEst = nextStateMarginal[j] * _backward[j].getSmoothed(i) / _backward[j]._totalMass;
					backwardEst = _backward[j].getSmoothed(i);
				else
					backwardEst = _backwardStateMarginal[i] / (double)numStates; //not sure about this
				be[i][j] = backwardEst;
				jSum += backwardEst;
			}
			for(int j=0; j<numStates;j++) {
				be[i][j] /= jSum;
			}
		}
		total = 0.0;
		for(int i=0; i<numStates; i++) {
			for(int j=0; j<numStates;j++) {
				double forwardEst = out[i][j];
				double backwardEst = be[i][j];
				diff += Math.abs(forwardEst - backwardEst);
				if(Math.abs(forwardEst - backwardEst) > 0.001) {
					System.out.println("diff on " +i + "," + j + " is " + Math.abs(forwardEst - backwardEst));
				}
				backwardEst = forwardEst;
				out[i][j] = (forwardEst + backwardEst)/2.0;
				total += out[i][j];
			}
		}
		System.out.println("trans model diff: " + diff);
		System.out.println("trans model total: " + total);
		return out;
	}
	

	protected TIntDoubleHashMap aggregateCounts(TIntArrayList occurs) {
		TIntDoubleHashMap out = new TIntDoubleHashMap();
		TIntIterator it = occurs.iterator();
		while(it.hasNext()) {
			int i = it.next();
			out.adjustOrPutValue(i, 1.0, 1.0);
		}
		return out;
	}

	//scans doc, resampling each variable.  Not used at test time.
	//returns number of changes
	protected int sampleDoc(int docId) {
		int changes = 0;
		TIntArrayList zs = _c._z[docId];
		TIntArrayList scratchZs = _c._scratchZ[docId];
		TIntArrayList doc = _c._docs[docId];
		//System.out.println("adding " + docCts.toString());
		for(int i=0; i<doc.size(); i++) {
			int newZ = sampleZ(docId, i, false);
			if(newZ != zs.get(i))
				changes++;
			scratchZs.set(i, newZ);
		}
		return changes;
	}
	
	public int sampleZ(int doc, int pos, boolean testing) {
		int w = _c._docs[doc].get(pos);
		int curZ = _c._z[doc].get(pos);
		SparseBackoffTree sbtForward = null;
		if(pos==0) {
			sbtForward = this._startState;
		}
		else {
			if(!testing && _SAMPLEONSAMPLES) { //HACK:we must sample left-to-right for this to work
				sbtForward = this._forward[_c._scratchZ[doc].get(pos-1)];
			}
			else 
				sbtForward = this._forward[_c._z[doc].get(pos-1)];
		}
		SparseBackoffTree sbtBackward = null;
		boolean forward = (testing && (pos==_c._docs[doc].size() - 1 || _c._z[doc].get(pos+1)==-1)); 
		if(!forward) { //don't use backward distribution if doing forward prediction
			if(pos==_c._docs[doc].size() - 1) {  
				sbtBackward = this._endState;
			}
			else {
				int nextZ = _c._z[doc].get(pos+1);
				if(nextZ >= 0) {
					sbtBackward = this._backward[nextZ];
				}
			}
		}
		SparseBackoffTree sbtWord = this._wordToState[w]; 
		boolean unseen = false;
		//if unseen, using simple smoothing distribution for the purpose of sampling:
		if(sbtWord._totalMass == 0.0) {
			unseen = true;
		}
		double bMarginalInv = 1.0f;
		double wMarginalInv = 1.0f;
		if(curZ >= 0) {
			bMarginalInv = 1.0 / this._backwardStateMarginal[curZ];
			wMarginalInv = 1.0 / this._wordStateMarginal[curZ];
		}
		
		SparseBackoffTree [] sbts;
		double [] subs;
		if(!forward) {
			if(unseen) {
				sbts = new SparseBackoffTree [] {sbtForward, sbtBackward};
				subs = new double [] {1.0, wMarginalInv};
			}
			else {
				sbts = new SparseBackoffTree [] {sbtWord, sbtForward, sbtBackward};
				subs = new double [] {wMarginalInv, 1.0, bMarginalInv};
			}
		}
		else {
			if(unseen) {
				sbts = new SparseBackoffTree [] {sbtForward};
				subs = new double [] {1.0};
			}
			else {
				sbts = new SparseBackoffTree [] {sbtForward, sbtWord};
				subs = new double [] {1.0, wMarginalInv};
			}
		}
		if(testing) //don't subtract at test time:
			subs = null;
		if(!testing && _SAMPLEONSAMPLES && pos > 0) {
			if(_c._z[doc].get(pos-1)==_c._scratchZ[doc].get(pos-1))
				subs[0] = 0.0;
		}
		SparseBackoffTreeIntersection sbti;
		if(subs!=null) 
			sbti = new SparseBackoffTreeIntersection(sbts, curZ, subs, true);
		else
			sbti = new SparseBackoffTreeIntersection(sbts, true);
		int sample = 0;
		if(baseProbs !=null && !unseen && !forward) {
			double wordOutProb = sbti._totalMass / sbti._siblingSbtis[0]._totalMass;
			wordOutProb *= (double)this._struct.numLeaves();
			wordOutProb /= this.get_NUMTOKENS();
			double rng = wordOutProb + baseProbs[doc][pos];
			double choice = _r.nextDouble() * rng;
			if(choice < wordOutProb)
				sample = sbti.sample(_r);
			else {
				sample = sbti._siblingSbtis[0].sample(_r); //note -- dangerous, depends on this siblingSbti data structure
					//containing the intersection of forward and backward distributions
				_numSkipWords++;
			}
		}
		else
			sample = sbti.sample(_r);
		return sample;
	}
	
    private long sampleRange(int start, int end) {
    	
    	long chg = 0;
    	
		for(int i=start; i<=end; i++) {
			chg += (long)sampleDoc(i);
		}

		return chg;
    }	
	

	public void initZ(TIntArrayList [] ds, int maxVal, boolean random, boolean scratch) {
		TIntArrayList [] zs;
		if(scratch) {
			_c._scratchZ = new TIntArrayList[ds.length];
			zs = _c._scratchZ;
		}
		else {
			_c._z = new TIntArrayList[ds.length];
			zs = _c._z;
		}
		for(int i=0; i<ds.length; i++) {
			zs[i] = new TIntArrayList(ds[i].size());
			for(int j=0; j<ds[i].size(); j++)
				zs[i].add(random ? _r.nextInt(maxVal) : 0); //don't use zero, reserved for start/end
		}
	}

	/**
	 * Returns a an array of sparse hash maps, one for each X, saying how many times the X maps to each state.
	 * if from = -2 then X = start of sentence
	 * if from = -1 then X = previous state (get state-to-next-state counts)
	 * if from = 0 then X = words (get word-to-state counts)
	 * if from = 1 then X = next state (get state-to-previous-state counts)
	 * if from = 2 then X = end of sentence (after </s> marker if any)
	 * @param from  See above
	 * @return
	 */
	public TIntDoubleHashMap [] aggregateCounts(int from, TIntArrayList [] zs, TIntArrayList [] ws) {
		TIntDoubleHashMap [] hm;
		
		//handle start/end as special case
		if(from==2 || from==-2) { 
			hm = new TIntDoubleHashMap[1];
			hm[0] = new TIntDoubleHashMap();
			for(int i=0; i<zs.length; i++) {
				if(from==-2) {
					hm[0].adjustOrPutValue(zs[i].get(0), 1.0, 1.0);
				}
				else {
					hm[0].adjustOrPutValue(zs[i].get(zs[i].size()-1), 1.0, 1.0);
				}
			}
			return hm;
		}
		
		if(from==0)
			hm = new TIntDoubleHashMap[_c._VOCABSIZE];
		else
			hm = new TIntDoubleHashMap[_struct.numLeaves()];
		
		//forward steps from word j to j+1
		//backward steps from j+1 to j
		//word steps to j
		for(int i=0; i<zs.length; i++) {
			for(int j=0; j<zs[i].size(); j++) {
				int xID = -1;
				int z = -1;
				if(j==zs[i].size()-1 && from !=0) //no transition here
					continue;
				if(from == -1) {
					xID = zs[i].get(j);
					z = zs[i].get(j+1);
				}
				else if(from == 1) {
					xID = zs[i].get(j+1);
					z = zs[i].get(j);
				}
				else {
					xID = ws[i].get(j);
					z = zs[i].get(j);
				}
				if(hm[xID] == null)
					hm[xID] = new TIntDoubleHashMap();
				hm[xID].adjustOrPutValue(z, 1.0, 1.0);
			}
		}
		return hm;
	}
	
	/**
	 * Returns sbts with given discounts of given type based on current sampler state
	 * if from = -2 then returns start-state parameteres
	 * if from = -1 then returns forward parameters
	 * if from = 0 then returns word-to-state parameters
	 * if from = 1 then backward parameters
	 * if from = 2 then returns end-state parameters
	 * @param dsWords	The discounts to use in the returned parameters
	 * @param from	see above
	 * @return
	 */
	public SparseBackoffTree [] getParamsFromZs(double [] dsWords, int from, TIntArrayList [] zs,
			TIntArrayList [] ws) {
		//aggregate:
		TIntDoubleHashMap [] hm = aggregateCounts(from, zs, ws);
		SparseBackoffTree [] out = new SparseBackoffTree[hm.length];
		
		for(int i=0; i<hm.length; i++) {
			out[i] = new SparseBackoffTree(_struct);
		}
		
		//add:
		for(int i=0; i<hm.length; i++) {
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
	
	/**
	 * reads the corpus and initializes zs and model
	 * @param inFile
	 * @param maxVal
	 * @return
	 * @throws Exception
	 */
	public int initializeForCorpus(String inFile, int maxVal) throws Exception {
		int toks = _c.readCorpusDat(inFile, true);
		setThreadBreaks();
		initZ(_c._docs, maxVal, true, false);
		updateModel(_c._z);
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
    	System.out.println("\ttime: " + stTime + "\tchanges: " + changes + "\tskipWords: " + _numSkipWords);
    	_numSkipWords = 0;
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
    public void updateModel(TIntArrayList [] zs) {
    	System.out.println("arch: " + Arrays.toString(this._branchingFactors));
		System.out.println("second doc samples: " + zs[1].toString());
		this._startState = getParamsFromZs(_forwardDelta, -2, zs, _c._docs)[0];
		this._forward = getParamsFromZs(_forwardDelta, -1, zs, _c._docs);
		this._backward = getParamsFromZs(_backwardDelta, 1, zs, _c._docs);
		this._endState = getParamsFromZs(_backwardDelta, 2, zs, _c._docs)[0];
		this._wordToState = getParamsFromZs(_wordDelta, 0, zs, _c._docs);
		SparseBackoffTree [] comb = Arrays.copyOf(_backward, _backward.length + 1);
		comb[_backward.length] = _endState;
		this._backwardStateMarginal = getNormalizers(comb, _struct);
		this._wordStateMarginal = getNormalizers(_wordToState, _struct);
		for(int i=0; i<_backward.length; i++) {
			_backward[i].divideCountsBy(_backwardStateMarginal);
		}
		this._endState.divideCountsBy(_backwardStateMarginal);
		for(int i=0; i<_wordToState.length; i++) {
			_wordToState[i].divideCountsBy(_wordStateMarginal);	
		}
		System.out.println("\tdiscounts forward: " + Arrays.toString(_forwardDelta) + 
				"\tbackward " + Arrays.toString(_backwardDelta) +
						"\tword " + Arrays.toString(_wordDelta));
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
	    	if(sum==0)
	    		return gradient;
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
	    				if(out[i] + curDeltas[i] - 1E-5 > 0.001) //don't do when hitting floor
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
    public void optimizeParameters(TIntArrayList [] zs) {
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

    	TIntDoubleHashMap [] hm = aggregateCounts(0, zs, _c._docs);
    	System.out.println("words:");
    	double step = STEPSIZE;
    	for(int i=0; i<NUMSTEPS; i++) {
        	gradientStep(_wordDelta, hm, step, false);
        	step *= STEPDEC;
    	}
    	long totalparams = 0L;
    	for(TIntDoubleHashMap i : hm) {
    		if(i != null)
    			totalparams += i.size();
    	}
    	hm = aggregateCounts(-1, zs, _c._docs);
    	TIntDoubleHashMap [] hm2 = Arrays.copyOf(hm, hm.length + 1);
    	hm2[hm.length] = aggregateCounts(-2, zs, _c._docs)[0];
    	hm = hm2;
    	long transParams = 0L;
    	for(TIntDoubleHashMap i : hm) {
    		if(i != null)
    			transParams += i.size();
    	}
    	System.out.println("forward:");
    	step = STEPSIZE;
    	for(int i=0; i<NUMSTEPS; i++) {
        	gradientStep(_forwardDelta, hm, step, false);
        	step *= STEPDEC;
    	}
    	hm = aggregateCounts(1, zs, _c._docs);
    	hm2 = Arrays.copyOf(hm, hm.length + 1);
    	hm2[hm.length] = aggregateCounts(2, zs, _c._docs)[0];
    	hm = hm2;
    	System.out.println("backward:");
    	step = STEPSIZE;
    	for(int i=0; i<NUMSTEPS; i++) {
        	gradientStep(_backwardDelta, hm, step, false);
        	step *= STEPDEC;
    	}
    	System.out.println("full word deltas: " + Arrays.toString(_wordDelta));
    	System.out.println("full forward deltas: " + Arrays.toString(_forwardDelta));
    	System.out.println("full backward deltas: " + Arrays.toString(_backwardDelta));
    	System.out.println("total nonzero word-state params: " + totalparams);
    	System.out.println("total forward trans params: " + transParams);
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
    		updateModel(_c._z);
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
     * @param  outFile   writes model to here before each expansion
     */
	public void trainModel(int iterations, int updateInterval, String outFile) throws Exception {
		boolean done = false;
		int j=0;
		int maxVal = _struct.numLeaves();
		while(!done) {
			for(int i=1; i<=iterations; i++) {
				System.out.print(i + ": ");
				initZ(_c._docs, maxVal, false, true); //init scratch array to zero
				gibbsPass();
				if((i % updateInterval)==0) { 
					optimizeParameters(_c._scratchZ);
				}
				updateModel(_c._scratchZ);
				_c._z = _c._scratchZ;
			}
			if(this._USEEXPANSION == 1) {
				writeModel(outFile + "." + j);
				j++;
				done = !expandModel(j);
			}
			else 
				done = true;
		}
	}
	
	public void writeModel(String outFile) throws Exception {
		SparseBackoffTree [] sbtW = this._wordToState;
		SparseBackoffTree [] sbtF = this._forward;
		SparseBackoffTree [] sbtB = this._backward;
		SparseBackoffTree sbtS = this._startState;
		SparseBackoffTree sbtE = this._endState;
		SparseBackoffTreeStructure struct = this._struct;
		_struct = null;
		_wordToState = null;
		_forward = null;
		_backward = null;
		_startState = null;
		_endState = null;
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outFile));
		oos.writeObject(this);
		oos.close();
		_wordToState = sbtW;
		_forward = sbtF;
		_backward = sbtB;
		_struct = struct;
		_startState = sbtS;
		_endState = sbtE;
	}
	
	public static SBTSequenceModel readModel(String inFile) throws Exception {
		ObjectInputStream ois = new ObjectInputStream(new FileInputStream(inFile));
		SBTSequenceModel out = (SBTSequenceModel) ois.readObject();
		ois.close();
		out._struct = new SparseBackoffTreeStructure(out._branchingFactors);
		out.updateModel(out._c._z);
		return out;
	}
	
	
	
	//returns model x word array of probabilities 
	public static double [][] getEnsembleProbsForDoc(int doc, double [][] wordTopicMarginal, SBTSequenceModel [] ms,
			double [][][] tms) {
		
		double [][] ps = new double[ms.length][ms[0]._c._docs[doc].size()]; //dims [model][token]
		for(int m=0; m<ms.length; m++) {
			int numStates = wordTopicMarginal[m].length;
			SBTSequenceModel mod = ms[m];
			TIntArrayList words = mod._c._docs[doc];
			double [][] dist = new double[words.size()+1][numStates];
			double [][] tm = tms[m];
			
			dist[0][0] = 1.0;
			for(int i=0; i<dist.length - 1; i++) {
				int w = words.get(i);
				//step from i to i+1:
				double [] pWordGivenState = new double[numStates];
				double [] pStateGivenPrev = new double[numStates];
//				for(int j=0; j<dist[i].length;j++) {
//					sbtForward.addWeighted(this._forward[j], dist[i][j]);
//				}
				//intersect:
				double p = 0.0;
				if(i==0) {
					for(int j=0; j<numStates;j++) {
						pStateGivenPrev[j] = ms[m]._startState.getSmoothed(j) / ms[m]._startState._totalMass;
					}
				}
				else {
					for(int j=0; j<numStates;j++) {
						for(int k=0; k<numStates;k++) {
							pStateGivenPrev[k] += tm[j][k]*dist[i][j];
						}
					}
				}
				for(int t=0;t<numStates; t++) {
					//dividing by wordTopicMarginal ensures we use a distribution P(w | t) that sums to one
					pWordGivenState[t] = mod._wordToState[w].getSmoothed(t) / wordTopicMarginal[m][t];
					//double numt = sbtForward.getSmoothed(t);
					//pStateGivenPrev[t] = numt / sbtForward._totalMass;
					dist[i+1][t] = (pWordGivenState[t]*pStateGivenPrev[t]);
					p += dist[i+1][t];
				}
				ps[m][i] = p;
				//renormalize:
				for(int t=0;t<numStates;t++) {
					dist[i+1][t] /= p;
				}
			}
		}
		return ps;
	}
	
	//side effect, fills psOut with probs per tok
	public static double [] testEnsembleOnDocExact(int doc, double [][] wordTopicMarginal, SBTSequenceModel [] ms,
			double [][][] tms, double [] ws, double [] psOut) {
		double LL = 0.0;
				
		double [][] ps = getEnsembleProbsForDoc(doc, wordTopicMarginal, ms, tms);
		
		double numWords = 0.0;
		double [] ensembleps = new double[ms[0]._c._docs[doc].size()];
		for(int m=0; m < ms.length; m++) {
			for(int i=0; i<ms[m]._c._docs[doc].size(); i++) {
				ensembleps[i] += ws[m]*ps[m][i];
			}
		}
		
		for(int i=0; i<ms[0]._c._docs[doc].size(); i++) {
			int w = ms[0]._c._docs[doc].get(i);
			if(ms[0]._c._pWord[w] > 0) {
				numWords++;
				LL += Math.log(ensembleps[i]);
				psOut[i] = ensembleps[i];
				if(Double.isNaN(LL))
					System.out.println("got NaN with " + ensembleps[i]);
			}
		}
		return new double [] {LL, numWords};		
	}
	
	public static double [] testEnsembleOnDocFullPpl(int doc, double [][] wordTopicMarginal, SBTSequenceModel [] ms) {
		double LL = 0.0;
		double [][] ps = new double[ms.length][ms[0]._c._docs[doc].size()]; //dims [model][token]
		for(int m=0; m<ms.length; m++) {
			for(int j=0; j<ms[m]._NUMPARTICLES; j++) {
				TIntArrayList words = ms[m]._c._docs[doc];
				for(int i=0; i<words.size(); i++) {
					ms[m]._c._z[doc].set(i, -1);
				}
				for(int i=0; i<words.size(); i++) {
					int w = words.get(i);
					double p = 0.0;						
					
					//resample everything before:
					for(int k=0; k<i; k++) {
						int r = ms[m].sampleZ(doc, k, true);
						ms[m]._c._z[doc].set(k, r);
					}
					
					if(ms[m]._c._pWord[w]>0.0) {
						//test on this word:
						double [] pWordGivenState = new double[wordTopicMarginal[m].length];
						double [] pStateGivenPrev = new double[wordTopicMarginal[m].length];
						int prev = 0;
						if(i > 0)
							prev = ms[m]._c._z[doc].get(i-1);
						for(int t=0;t<wordTopicMarginal[m].length; t++) {
							//dividing by wordTopicMarginal ensures we use a distribution P(w | t) that sums to one
							
							pWordGivenState[t] = ms[m]._wordToState[w].getSmoothed(t) / wordTopicMarginal[m][t];
							if(ms[m]._forward[prev]._totalMass > 0.0) {
								double numt = ms[m]._forward[prev].getSmoothed(t);
								pStateGivenPrev[t] = numt / ms[m]._forward[prev]._totalMass;
							}
							else {
								pStateGivenPrev[t] = 1.0 / (double)pStateGivenPrev.length;
							}
							p += (pWordGivenState[t]*pStateGivenPrev[t]);
							if(Double.isNaN(p))
								System.err.println("nan.");
							if(p < 0.0)
								System.err.println("negative.");
						}
						ps[m][i] += p;
					}
					//sample this word:
					int r = ms[m].sampleZ(doc, i, true);
					ms[m]._c._z[doc].set(i, r);
				}
			}
		}
		double numWords = 0.0;
		double [] ensembleps = new double[ms[0]._c._docs[doc].size()];
		for(int m=0; m < ms.length; m++) {
			for(int i=0; i<ms[m]._c._docs[doc].size(); i++) {
				ensembleps[i] += ps[m][i]/((double)ms[m]._NUMPARTICLES * (double)ms.length);
			}
		}
		
		for(int i=0; i<ms[0]._c._docs[doc].size(); i++) {
			int w = ms[0]._c._docs[doc].get(i);
			if(ms[0]._c._pWord[w] > 0) {
				numWords++;
				LL += Math.log(ensembleps[i]);
				if(Double.isNaN(LL))
					System.out.println("got NaN with " + ensembleps[i]);
			}
		}
		return new double [] {LL, numWords};		
	}
	
	//currently requires all test words found in vocab
	public double [] testOnDocFullPplExact(int doc, double [] wordTopicMarginal, double [][] tm) {
		double LL = 0.0;
		double [] ps = new double[_c._docs[doc].size()];
		
		TIntArrayList words = _c._docs[doc];
		int numStates = wordTopicMarginal.length;
		double [][] dist = new double[words.size()+1][numStates];
		
		
		dist[0][0] = 1.0;
		for(int i=0; i<dist.length - 1; i++) {
			int w = words.get(i);
			//step from i to i+1:
			double [] pWordGivenState = new double[numStates];
			double [] pStateGivenPrev = new double[numStates];
//			for(int j=0; j<dist[i].length;j++) {
//				sbtForward.addWeighted(this._forward[j], dist[i][j]);
//			}
			//intersect:
			double p = 0.0;
			double checksum = 0.0;
			if(i==0) {
				for(int j=0; j<numStates;j++) {
					pStateGivenPrev[j] = _startState.getSmoothed(j) / _startState._totalMass;
					checksum += pStateGivenPrev[j]; 
				}
			}
			else {
				for(int j=0; j<numStates;j++) {
					for(int k=0; k<numStates;k++) {
						pStateGivenPrev[k] += tm[j][k]*dist[i][j];
					}
				}
			}
			for(int t=0;t<numStates; t++) {
				//dividing by wordTopicMarginal ensures we use a distribution P(w | t) that sums to one
				pWordGivenState[t] = _wordToState[w].getSmoothed(t) / wordTopicMarginal[t];
				//double numt = sbtForward.getSmoothed(t);
				//pStateGivenPrev[t] = numt / sbtForward._totalMass;
				dist[i+1][t] = (pWordGivenState[t]*pStateGivenPrev[t]);
				p += dist[i+1][t];
			}
			ps[i] = p;
			//renormalize:
			for(int t=0;t<numStates;t++) {
				dist[i+1][t] /= p;
			}
		}
		
		for(int i=0; i<ps.length; i++) {
			LL += Math.log(ps[i]);
			if(Double.isNaN(LL)) {
				System.out.println("got NaN with " + ps[i]);
			}
		}
		return new double [] {LL, ps.length};
	}
	
	//tests on a single document using left-to-right method
	//word topic marginal (i.e. P(topic) computed from P(topic, word)) is supplied to ensure evaluated distribution sums to one
	public double [] testOnDocFullPpl(int doc, double [] wordTopicMarginal) {
		boolean CHECKSUMTOONE = false; //slow.  Can use to confirm that test code uses a valid probability distribution over words
		double LL = 0.0;
		double [] ps = new double[_c._docs[doc].size()];
		for(int j=0; j<_NUMPARTICLES; j++) {
			TIntArrayList words = _c._docs[doc];
			for(int i=0; i<words.size(); i++) {
				_c._z[doc].set(i, -1);
			}
			for(int i=0; i<words.size(); i++) {
				int w = words.get(i);
				double p = 0.0;						
				
				//resample everything before:
				for(int k=0; k<i; k++) {
					int r = sampleZ(doc, k, true);
					_c._z[doc].set(k, r);
				}
				
				if(_c._pWord[w]>0.0) {
					//test on this word:
					double [] pWordGivenState = new double[wordTopicMarginal.length];
					double [] pStateGivenPrev = new double[wordTopicMarginal.length];
					double checkSum = 0.0f;
					int prev = 0;
					if(i > 0)
						prev = _c._z[doc].get(i-1);
					for(int t=0;t<wordTopicMarginal.length; t++) {
						//dividing by wordTopicMarginal ensures we use a distribution P(w | t) that sums to one
						
						pWordGivenState[t] = _wordToState[w].getSmoothed(t) / wordTopicMarginal[t];
						if(this._forward[prev]._totalMass > 0.0) {
						double numt = this._forward[prev].getSmoothed(t);
						checkSum += numt;
						pStateGivenPrev[t] = numt / this._forward[prev]._totalMass;
						}
						else {
							pStateGivenPrev[t] = 1.0 / (double)pStateGivenPrev.length;
						}
						p += (pWordGivenState[t]*pStateGivenPrev[t]);
						if(Double.isNaN(p))
							System.err.println("nan.");
						if(p < 0.0)
							System.err.println("negative.");
					}
					if(j==0 && (i==0 || (i== _c._docs[doc].size()-1)) && (CHECKSUMTOONE)) {
						float sDoc = 0.0f;
						float sWord = 0.0f;
						float totalTopicMarginal = 0.0f;
						for(int t=0; t<pWordGivenState.length; t++) {
							sDoc += pStateGivenPrev[t];
							sWord += pWordGivenState[t]*wordTopicMarginal[t];
							totalTopicMarginal += wordTopicMarginal[t];
						}
						sWord /= totalTopicMarginal;
						System.out.println("Check SUM: " + i + "\t" + p + "\t" + sDoc + "\t" + sWord + "\t" + _c._pWord[w]);
					}
				}
				ps[i] += p;
				//sample this word:
				int r = sampleZ(doc, i, true);
				_c._z[doc].set(i, r);
			}
		}
		double numWords = 0.0;
		for(int i=0; i<_c._docs[doc].size(); i++) {
			int w = _c._docs[doc].get(i);
			if(_c._pWord[w] > 0) {
				numWords++;
				LL += Math.log(ps[i]/(double)_NUMPARTICLES);
				if(Double.isNaN(LL))
					System.out.println("got NaN with " + ps[i]);
			}
		}
		return new double [] {LL, numWords};
	}

	public static void train(String inputFile, String outputFile, String configFile, String probsFile) throws Exception {
		SBTSequenceModel sbtsm = new SBTSequenceModel(configFile);
		sbtsm.initializeForCorpus(inputFile, sbtsm._struct.numLeaves());
		if(probsFile != null) {
			sbtsm.baseProbs = sbtsm.readBaseProbs(probsFile);
		}
		sbtsm.trainModel(sbtsm._NUMITERATIONS, sbtsm._OPTIMIZEINTERVAL, outputFile);
		sbtsm.writeModel(outputFile);
	}
	
	public static void train(String inputFile, String outputFile, String configFile) throws Exception {
		train(inputFile, outputFile, configFile, null);
	}

	//ps has dims docs x models x tokens
	public static double [] optimizeEnsembleWeights(double [][][] ps) {
		double [] ws = new double[ps[0].length];
		Arrays.fill(ws, 1.0/(double)ws.length);
		double STEPSIZE = 0.0000001;
		for(int iter = 0; iter<10000; iter++) {
			double [] grad = new double[ws.length];
			double LL = 0.0;
			double [] wordGrad = new double[ws.length];
			for(int i=0; i<ps.length; i++) {
				//TODO: fix this non-sequential array traverse
				for(int j=0; j<ps[i][0].length; j++) {
					double denom = 0.0;
					for(int m=0; m<ps[i].length; m++) {
						wordGrad[m] = ps[i][m][j];
						denom += ps[i][m][j]*ws[m];
					}
					for(int m=0; m<ps[i].length; m++) {
						grad[m] += wordGrad[m]/denom;
					}
					LL += Math.log(denom);
				}
			}
			//HACK: make gradient sum to zero:
			double gradSum = 0.0;
			for(int m=0; m<grad.length; m++) {
				gradSum += grad[m];
			}
			for(int m=0; m<grad.length; m++) {
				grad[m] -= gradSum / (double)grad.length;
			}
			double normalizer = 0.0;
			for(int i=0; i<ws.length;i++) {
				ws[i] += STEPSIZE*grad[i];
				normalizer += ws[i]; //not actually necessary if gradient sums to zero
			}
			for(int i=0; i<ws.length;i++) {
				ws[i] /= normalizer;
			}
			if(iter % 100 == 0) {
				System.out.println(iter + " LL is : " + LL + " Grad is " + Arrays.toString(grad));
				System.out.println("new weights: " + Arrays.toString(ws));
			}
		}
		return ws;
	}
	
	public static double [] trainEnsemble(SBTSequenceModel [] sbtsm, double [][] wordStateNormalizer, double [][][] tm, String validationFile, int numDocs) throws Exception {
	  System.out.println("train ens num docs: " + numDocs);
		for(int i=0; i<sbtsm.length; i++) {
			sbtsm[i]._c.reInitForCorpus(validationFile, numDocs);
		}
		double LL = 0.0;
		double numWords = 0.0;
		TestGetProbsDoer [] tds = new TestGetProbsDoer[numDocs];
		ExecutorService e = Executors.newFixedThreadPool(sbtsm[0]._NUMTHREADS);
		double [] eqW = new double[sbtsm.length];
		Arrays.fill(eqW,  1.0/(double)eqW.length);
		//ExecutorService e = Executors.newFixedThreadPool(1);
		for(int i=0; i<tds.length;i++) {
			if(sbtsm[0]._MAXTESTDOCSIZE < 0 || sbtsm[0]._c._docs[i].size() > sbtsm[0]._MAXTESTDOCSIZE) {
				System.out.println("skipping " + i + " size " + sbtsm[0]._c._docs[i].size());
				continue;
			}
			tds[i] = new TestGetProbsDoer(i, wordStateNormalizer, sbtsm, tm, eqW);
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
		double [][][] ps = new double[numDocs][][];
		for(int i=0; i<numDocs;i++) {
			ps[i] = tds[i]._probs;
		}
		double [] outW = optimizeEnsembleWeights(ps);
		return outW;		
	}
	
	//computes log likelihood of model using left-to-right method
	//returns {ppl, number of tested words}
	//numDocs is number in test file
	//maxDocs is number actually tested (starting from the beginning of the file)
	public static double [] testEnsemble(File [] modelFile, String testFile, int numDocs, int maxDocs, String configFile,
			String validationFile, int validationDocs, String probsOutfile) throws Exception {
		SBTSequenceModel [] sbtsm = new SBTSequenceModel[modelFile.length];
		for(int i=0; i< sbtsm.length; i++) {
			sbtsm[i] = readModel(modelFile[i].getAbsolutePath());
			sbtsm[i].readConfigFile(configFile);
		}
		
		double [][] wordStateNormalizer = new double[modelFile.length][];
		double [][][] tm = new double[modelFile.length][][];
		int [] numStates = new int[sbtsm.length];
		for(int i=0; i< sbtsm.length; i++) {
			wordStateNormalizer[i] = getCheckMarginal(sbtsm[i]._wordToState, sbtsm[i]._struct);
			numStates[i] = wordStateNormalizer[i].length;
			tm[i] = sbtsm[i].transModel();
		}
		
		double [] ws = new double[sbtsm.length];
		if(validationFile != null) {
			ws = trainEnsemble(sbtsm, wordStateNormalizer, tm, validationFile, validationDocs);
		}
		else {
			Arrays.fill(ws, 1.0/(double)ws.length);
		}
		for(int i=0; i<sbtsm.length; i++) {
			sbtsm[i]._c.reInitForCorpus(testFile, numDocs);
		}
		double LL = 0.0;
		double numWords = 0.0;
		TestEnsembleExactDoer [] tds = new TestEnsembleExactDoer[maxDocs];
		ExecutorService e = Executors.newFixedThreadPool(sbtsm[0]._NUMTHREADS);
		//ExecutorService e = Executors.newFixedThreadPool(1);
		for(int i=0; i<tds.length;i++) {
			if(sbtsm[0]._MAXTESTDOCSIZE < 0 || sbtsm[0]._c._docs[i].size() > sbtsm[0]._MAXTESTDOCSIZE) {
				System.out.println("skipping " + i + " size " + sbtsm[0]._c._docs[i].size());
				continue;
			}
			tds[i] = new TestEnsembleExactDoer(i, wordStateNormalizer, sbtsm, tm, ws);
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
		double [][] psOut = new double[sbtsm[0]._c._docs.length][];
		for(int i=0; i<tds.length; i++) {
			if(tds[i]==null)
				continue;
			LL += tds[i]._res[0];
			numWords += tds[i]._res[1];
			psOut[i] = tds[i]._ps;
		}
		System.out.println("Test LL: " + LL);
		System.out.println("Words: " + numWords);
		System.out.println("ppl: " + Math.exp(-LL/numWords));
		if(probsOutfile != null) {
			writeBaseProbs(probsOutfile, psOut);
		}
		return new double [] {LL, numWords};
	}
	
	//computes log likelihood of model using left-to-right method
	//returns {ppl, number of tested words}
	//numDocs is number in test file
	//maxDocs is number actually tested (starting from the beginning of the file)
	public static double [] testModel(String modelFile, String testFile, int numDocs, int maxDocs, String configFile) throws Exception {
		boolean CHECKWORDDISTRIBUTION = false;
		
		SBTSequenceModel sbtsm = readModel(modelFile);
		sbtsm.readConfigFile(configFile);
		double [] wordStateNormalizer = getCheckMarginal(sbtsm._wordToState, sbtsm._struct); 
		int numStates = wordStateNormalizer.length;
		double [][] tm = sbtsm.transModel();

		float [] pWordGivenState = new float[wordStateNormalizer.length] ; 
		if(CHECKWORDDISTRIBUTION) {
			for(int w=0; w<sbtsm._wordToState.length; w++) {
				for(int t=0;t<wordStateNormalizer.length; t++) {
					pWordGivenState[t] += sbtsm._wordToState[w].getSmoothed(t) / wordStateNormalizer[t];
				}
			}
			float sum = 0.0f;
			for(int i=0; i<pWordGivenState.length; i++) {
				sum += pWordGivenState[i];
			}
			
			System.out.println("Word sum: " + sum + " should be about " + numStates);
		}
		
		sbtsm._c.reInitForCorpus(testFile, numDocs);
		sbtsm.setThreadBreaks();
		double LL = 0.0;
		double numWords = 0.0;
		//TestDoer [] tds = new TestDoer[maxDocs];
		TestExactDoer [] tds = new TestExactDoer[maxDocs];
		ExecutorService e = Executors.newFixedThreadPool(sbtsm._NUMTHREADS);
		//ExecutorService e = Executors.newFixedThreadPool(1);
		for(int i=0; i<tds.length;i++) {
			if(sbtsm._MAXTESTDOCSIZE < 0 || sbtsm._c._docs[i].size() > sbtsm._MAXTESTDOCSIZE) {
				System.out.println("skipping " + i + " size " + sbtsm._c._docs[i].size());
				continue;
			}
			tds[i] = sbtsm.new TestExactDoer(i, wordStateNormalizer, tm);
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
		System.out.println("Test LL: " + Arrays.toString(ll));
		System.out.println("ppl: " + Math.exp(-ll[0]/ll[1]));
	}
	
	/**
	 * Outputs the sparse word-to-topic vector proportional to P(t | w)/P(t) for words and topics with positive counts
	 * @param modelFile	Contains the model
	 * @param outputFile	Where to write the output
	 * @throws Exception
	 */
	public static void outputWordToTopicFile(String modelFile, String outputFile) throws Exception {
		BufferedWriter bwOut = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), "UTF8"));
		SBTSequenceModel sbtsm = readModel(modelFile);
		for(int i=0; i<sbtsm._c._pWord.length; i++) {
			if(sbtsm._c._pWord[i] > 0.0) {
				bwOut.write(i + "\t");
				TIntDoubleHashMap hm = sbtsm._wordToState[i].getLeafCounts();
				TIntDoubleIterator it = hm.iterator();
				while(it.hasNext()) {
					it.advance();
					int wId = it.key();
					double val = it.value();
					bwOut.write(wId + ":" + val + " ");
				}
				bwOut.write("\r\n");
			}
		}
		bwOut.close();
	}
	
	public void testAggregateCounts() {
		TIntArrayList [] oneZs = new TIntArrayList[1];
		oneZs[0] = new TIntArrayList();
		oneZs[0].add(1);
		oneZs[0].add(2);
		oneZs[0].add(3);
		oneZs[0].add(4);
		oneZs[0].add(5);
		TIntArrayList [] oneWs = new TIntArrayList[1];
		oneWs[0] = new TIntArrayList();
		oneWs[0].add(3);
		oneWs[0].add(3);
		oneWs[0].add(3);
		oneWs[0].add(3);
		oneWs[0].add(3);
		TIntDoubleHashMap [] hm = aggregateCounts(-1, oneZs, oneWs);
		System.out.println("hm is " + hm.length);
		hm = aggregateCounts(1, oneZs, oneWs);
		System.out.println("hm is " + hm.length);
		hm = aggregateCounts(0, oneZs, oneWs);
		System.out.println("hm is " + hm.length);
	}
	
	public static void main(String[] args) throws Exception {
		if(args.length > 0 && args[0].equalsIgnoreCase("train")) {
			if(args.length != 4) {
				System.err.println("Usage: train <input_file> <model_output_file> <configuration_file>");
				return;
			}
			train(args[1], args[2], args[3]);
		}
		else if(args.length > 0 && args[0].equalsIgnoreCase("trainAug")) { //augment existing probs
			if(args.length != 5) {
				System.err.println("Usage: train <input_file> <model_output_file> <configuration_file> <probs_file>");
				return;
			}
			train(args[1], args[2], args[3], args[4]);
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
		else if(args.length > 0 && args[0].equalsIgnoreCase("testEns")) {
			if(args.length != 6 && args.length != 8) {
				System.err.println("Usage: test <model_dir> <test_file> <configuration_file> <num_docs_in_test_file> <num_docs_to_test>"
						+ " [<validation_file> <num_validation_docs>]");
			}
			File f = new File(args[1]);
			if(args.length==6)
				testEnsemble(f.listFiles(), args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3], null, 0, null);
			else
				testEnsemble(f.listFiles(), args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3], args[6], 
						Integer.parseInt(args[7]), null);
		}
		else if(args.length > 0 && args[0].equalsIgnoreCase("outputEnsProbs")) {
			if(args.length != 7 && args.length != 9) {
				System.err.println("Usage: test <model_dir> <test_file> <configuration_file> <num_docs_in_test_file> <num_docs_to_test>"
						+ " <probs_out_file> [<validation_file> <num_validation_docs>]");
			}
			File f = new File(args[1]);
			if(args.length==7)
				testEnsemble(f.listFiles(), args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3], null, 0, args[6]);
			else
				testEnsemble(f.listFiles(), args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3], args[7], 
						Integer.parseInt(args[8]), args[6]);			
		}
		else {
			System.err.println("Usage: <train|trainAug|test|wordreps|testEns|outputEnsProbs> <args...>");
		}
	}

	
}
