package edu.northwestern.cs.websail.sbt;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntDoubleHashMap;

public class SBTDocumentModel implements Serializable {

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private class GibbsDoer implements Runnable {
    	int _start;
    	int _end;
    	int _changes = 0;
    	
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
	
	//corpus info: (TODO: factor out as separate corpus class)
	private int _VOCABSIZE = -1;  
	private int _NUMDOCS = -1;
	private int _NUMTOKENS = -1;
	TIntArrayList [] _docs; //the words
	TIntArrayList [] _z; //the topic variable assignments
	double [] _pWord; //marginal word prob

	//Model information:
	int [] _branchingFactors;
	SparseBackoffTreeStructure _struct = null;
	//SparseBackoffTree [] _topicGivenDoc = null;
	SparseBackoffTree [] _topicGivenWord = null; 
	double [] _topicMarginal;
	double [] _docsDeltaEndPts = new double [] {0.1, 0.1};
	double [] _wordsDeltaEndPts = new double [] {0.1, 0.1};
	
	//threading:
	private int _NUMTHREADS = -1;
	private int[] _THREADBREAKS = null; //inclusive doc indices where threads *end* (initial thread starts
							//at 0 implicitly)
	
	//training config:
	private int _NUMITERATIONS = -1;
	private int _OPTIMIZEINTERVAL = -1;
	
	public int[] get_branchingFactors() {
		return _branchingFactors;
	}


	public void set_branchingFactors(int[] _branchingFactors) {
		this._branchingFactors = _branchingFactors;
	}


	public SBTDocumentModel(String configFile) throws Exception {
		readConfigFile(configFile);
		_struct = new SparseBackoffTreeStructure(_branchingFactors);
	}
	
	
	public int get_VOCABSIZE() {
		return _VOCABSIZE;
	}

	public void set_VOCABSIZE(int _VOCABSIZE) {
		this._VOCABSIZE = _VOCABSIZE;
	}

	public int get_NUMDOCS() {
		return _NUMDOCS;
	}

	public void set_NUMDOCS(int _NUMDOCS) {
		this._NUMDOCS = _NUMDOCS;
	}

	public int get_NUMTOKENS() {
		return _NUMTOKENS;
	}

	public void set_NUMTOKENS(int _NUMTOKENS) {
		this._NUMTOKENS = _NUMTOKENS;
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
		int approxToks = _NUMTOKENS / _NUMTHREADS;
		int ct = 0;
		int thread = 0;
		for(int i=0; i<_NUMDOCS; i++ ) {
			ct += _docs[i].size();
			if(ct > approxToks) {
				_THREADBREAKS[thread++] = i;
				ct = 0;
			}
		}
		//extra goes in last thread:
		_THREADBREAKS[_NUMTHREADS - 1] = _NUMDOCS - 1;
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
	
	//scans doc, resampling each variable.  Not used at test time.
	//returns number of changes
	private int sampleDoc(int docId) {
		int changes = 0;
		TIntArrayList zs = _z[docId];
		TIntArrayList doc = _docs[docId];
		TIntDoubleHashMap docCts = aggregateCounts(zs);
		SparseBackoffTree sbtDoc = new SparseBackoffTree(_struct);
		//System.out.println("adding " + docCts.toString());
		double [] docDs = interpolateEndPoints(this._docsDeltaEndPts, this._branchingFactors.length);
		sbtDoc.addAllMass(docCts, docDs);
		for(int i=0; i<doc.size(); i++) {
			int newZ = sampleZ(docId, i, false, sbtDoc);
			if(newZ != zs.get(i))
				changes++;
			zs.set(i, newZ);
		}
		return changes;
	}
	
	public int sampleZ(int doc, int pos, boolean testing, SparseBackoffTree topicGivenDoc) {
		int w = _docs[doc].get(pos);
		int curZ = _z[doc].get(pos);
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
	
    private int sampleRange(int start, int end) {
    	
    	int chg = 0;
    	
		for(int i=start; i<=end; i++) {
			chg += sampleDoc(i);
		}

		return chg;
    }
	
	public static TIntArrayList lineToList(String sLine) {
		String [] idFeats = sLine.split("\t");
		String feats = idFeats[1];
		String [] featVals = feats.split(" ");
		TIntArrayList al = new TIntArrayList(featVals.length);
		for(int i=0; i<featVals.length; i++) {
			String [] featVal = featVals[i].split(":");
			int feat = Integer.parseInt(featVal[0]);
			int val = Integer.parseInt(featVal[1]);
			for(int j=0; j<val; j++) {
				al.add(feat);
			}
		}
		return al;
	}
	
	
	//returns number of docs
	public int readCorpusDat(String inFile, boolean updatePWord) throws Exception {
		BufferedReader brIn = new BufferedReader(
				new InputStreamReader(new FileInputStream(inFile), "UTF8" ));
		String sLine;
		int i=0;
		int toks = 0;
		ArrayList<TIntArrayList> aldocs = new ArrayList<TIntArrayList>();
		if (updatePWord)
			_pWord = new double[_VOCABSIZE];
		while((sLine = brIn.readLine())!=null) {
			TIntArrayList ll = lineToList(sLine);
			aldocs.add(ll);
			if (updatePWord) {
				TIntIterator it = ll.iterator();
				while(it.hasNext()) {
					_pWord[it.next()]++;
				}
			}
			toks += ll.size();
		}
		brIn.close();
		System.out.println("Tokens: " + toks);
		_NUMTOKENS = toks;
		double dubToks = (double)toks;
		if(updatePWord)
			for(int j=0; j<_pWord.length; j++)
				if(_pWord[j]>0.0)
					_pWord[j] /= dubToks;
		_NUMDOCS = aldocs.size();
		_docs = aldocs.toArray(new TIntArrayList[aldocs.size()]);
		setThreadBreaks();
		return i;
	}
	
	//keeps word distributions, but re-sets docs, counts, and z's for new corpus
	//takes in number of docs in new corpus
	//keeps topicMarginal and _pWord fixed
	public void reInitForCorpus(String testFile, int numDocs) throws Exception {
		_docs = new TIntArrayList[numDocs];
		_NUMDOCS = numDocs;
		this.readCorpusDat(testFile, false);
		_z = new TIntArrayList[numDocs];
		for(int i=0; i<_docs.length; i++) {
			_z[i] = new TIntArrayList(_docs[i].size());
			for(int j=0; j<_docs[i].size(); j++) {
				int w = _docs[i].get(j);
				int r = -1;
				_z[i].add(r);
			}
		}
		setThreadBreaks();
	}
	
	public void initZRandom(TIntArrayList [] ds, int maxVal) {
		_z = new TIntArrayList[ds.length];
		for(int i=0; i<ds.length; i++) {
			_z[i] = new TIntArrayList(ds[i].size());
			for(int j=0; j<ds[i].size(); j++)
				_z[i].add(_r.nextInt(maxVal));
		}
	}
	
	public void updateWordParamsFromZs(double [] dsWords) {
		_topicGivenWord = new SparseBackoffTree[_VOCABSIZE];
		for(int i=0; i<_VOCABSIZE; i++) {
			_topicGivenWord[i] = new SparseBackoffTree(_struct);
		}
		//aggregate:
		TIntDoubleHashMap [] hm = new TIntDoubleHashMap[_VOCABSIZE];
		for(int i=0; i<_z.length; i++) {
			for(int j=0; j<_z[i].size(); j++) {
				int wordID = _docs[i].get(j);
				int z = _z[i].get(j);
				if(hm[wordID] == null)
					hm[wordID] = new TIntDoubleHashMap();
				hm[wordID].adjustOrPutValue(z, 1.0, 1.0);
			}
		}
		//add:
		for(int i=0; i<_VOCABSIZE; i++) {
			if(hm[i] == null)
				continue;
			TIntDoubleIterator it = hm[i].iterator();
			while(it.hasNext()) {
				it.advance();
				int z = it.key();
				double val = it.value();
				_topicGivenWord[i].smoothAndAddMass(z, val, dsWords);
			}
		}
	}
	
	public static double [] interpolateEndPoints(double [] endPts, int length) {
		double [] out = new double[length];
		double delt = (endPts[1] - endPts[0])/(length - 1);
		for(int i=0; i<length; i++) {
			out[i] = delt*i + endPts[0];
		}
		return out;
	}
	
	//reads the corpus and initializes zs and model
	public int initializeForCorpus(String inFile, int maxVal) throws Exception {
		int toks = readCorpusDat(inFile, true);
		initZRandom(_docs, maxVal);
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
	//such that if we divide leaf counts by amt_i, we get "proper smoothing"
    //TODO: improve this comment
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
		double maxSmooth = 0.0f;
		double sumCount = 0.0f;
		double sumSmoother = 0.0f;
		for(int i=0; i<numStates; i++) {
			double [] smoothAndCount = shdAgg.getSmoothAndCount(i);
			out[i] += smoothAndCount[0] + smoothAndCount[1];
		}
		return out;
	}
	
    public void updateModel() {
		System.out.println("first doc samples: " + _z[0].toString());
    	this.updateWordParamsFromZs(this.interpolateEndPoints(_wordsDeltaEndPts, _branchingFactors.length));
    	
		_topicMarginal = getNormalizers(_topicGivenWord, _struct);
		for(int i=0; i<_topicGivenWord.length; i++) 
			_topicGivenWord[i].divideCountsBy(_topicMarginal);
		System.out.println("\tdiscounts doc: " + Arrays.toString(this._docsDeltaEndPts) + "\tword: " + Arrays.toString(this._wordsDeltaEndPts));
    }
    
    public void optimizeParameters() {
    	System.out.println("here is where we would optimize parameters.");
    	this._docsDeltaEndPts[0] *= 0.7;
    	this._docsDeltaEndPts[1] *= 0.95;
    	this._wordsDeltaEndPts[0] *= 0.7;
    	this._wordsDeltaEndPts[1] *= 0.95;
    }
    
	public void trainModel(int iterations, int updateInterval) {
		
		for(int i=1; i<=iterations; i++) {
			System.out.print(i + ": ");
			gibbsPass();
			if((i % updateInterval)==0) { 
				optimizeParameters();
			}
			updateModel();
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
	
	public double [] testOnDocFullPpl(int doc, double [] ds, double [] wordTopicMarginal) {
		boolean CHECKSUMTOONE = false;
		double LL = 0.0;
		double numWords = 0.0;
		final int PARTICLES = 20;
		ArrayList<Integer> wordOrder = new ArrayList<Integer>();
		for(int i=0; i<_docs[doc].size(); i++) {
			wordOrder.add(i);
		}
		Collections.shuffle(wordOrder);
		double [] ps = new double[wordOrder.size()];
		for(int j=0; j<PARTICLES; j++) {
			SparseBackoffTree topicGivenDoc = new SparseBackoffTree(_struct);
			numWords = 0.0;
			for(int i=0; i<_docs[doc].size(); i++) {
				_z[doc].set(i, -1);
			}
			for(int i=0; i<_docs[doc].size(); i++) {
				int widx = wordOrder.get(i);
				int w = _docs[doc].get(widx);
				double p = 0.0;
				//only compute for seen words:
				if(_pWord[w] > 0) {
					//init topicGivenDoc with all samples:
					topicGivenDoc = new SparseBackoffTree(_struct);
					for(int k =0; k<_z[doc].size(); k++) {
						int r = _z[doc].get(k);
						if(r != -1)
							topicGivenDoc.smoothAndAddMass(r, 1.0, ds);
					}
						//resample everything before:
					for(int k=0; k<i; k++) {
						int wkidx = wordOrder.get(k);
						if(_z[doc].get(wkidx) != -1) { //seen word
							int r = sampleZ(doc, wkidx, true, topicGivenDoc);
							_z[doc].set(wkidx, r);
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
					if((i==0 || (i== _docs[doc].size()-1)) && (CHECKSUMTOONE)) {
						float sDoc = 0.0f;
						float sWord = 0.0f;
						float totalTopicMarginal = 0.0f;
						for(int t=0; t<pWordGivenState.length; t++) {
							sDoc += pStateGivenDoc[t];
							sWord += pWordGivenState[t]*wordTopicMarginal[t];
							totalTopicMarginal += wordTopicMarginal[t];
						}
						sWord /= totalTopicMarginal;
						System.out.println("Check SUM: " + i + "\t" + p + "\t" + sDoc + "\t" + sWord + "\t" + _pWord[w]);
					}
					ps[widx] += p;
					//sample this word:
					int r = -1;
					if(numWords>0)
						r = sampleZ(doc, widx, true, topicGivenDoc);
					else
						r = _topicGivenWord[w].sample(_r);
					_z[doc].set(widx, r);
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
		for(int i=0; i<_docs[doc].size(); i++) {
			int widx = wordOrder.get(i);
			int w = _docs[doc].get(widx);
			if(_pWord[w] > 0) {
				LL += Math.log(ps[widx]/(double)PARTICLES);
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
	public static double [] testModel(String modelFile, String testFile, int numDocs, int maxDocs) throws Exception {
		boolean CHECKWORDDISTRIBUTION = true;
		
		SBTDocumentModel sbtdm = readModel(modelFile);
		double [] ds = interpolateEndPoints(sbtdm._docsDeltaEndPts, sbtdm._branchingFactors.length);
		double numWordTrainToks = (double)sbtdm._NUMTOKENS;
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

		sbtdm.reInitForCorpus(testFile, numDocs);
		double LL = 0.0;
		double numWords = 0.0;
		double origWords = 0.0;
		TestDoer [] tds = new TestDoer[maxDocs];
		ExecutorService e = Executors.newFixedThreadPool(sbtdm._NUMTHREADS);
		//ExecutorService e = Executors.newFixedThreadPool(1);
		for(int i=0; i<tds.length;i++) {
			if(sbtdm._docs[i].size() > 2000) {
				System.out.println("skipping " + i);
				continue;
			}
			tds[i] = sbtdm.new TestDoer(i, ds, wordTopicNormalizer);
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
		System.out.println("total orig words: " + origWords);
		return new double [] {LL, numWords};
	}
	
	public static void test(String modelFile, String inputFile, int numDocsInFile, int numDocsToTest) throws Exception {
		double [] ll = testModel(modelFile, inputFile, numDocsInFile, numDocsToTest);
		System.out.println("Test ppl: " + Arrays.toString(ll));
	}
	
	public static void main(String[] args) throws Exception {
		if(args.length > 0 && args[0].equalsIgnoreCase("train")) {
			if(args.length != 4) {
				System.err.println("Usage: SBTDocumentModel train <input_file> <model_output_file> <configuration_file>");
				return;
			}
			train(args[1], args[2], args[3]);
		}
		else if(args.length > 0 && args[0].equalsIgnoreCase("test")) {
			if(args.length != 5) {
				System.err.println("Usage: SBTDocumentModel test <model_file> <test_file> <num_docs_in_file> <num_docs_to_test>");
			}
			test(args[1], args[2], Integer.parseInt(args[3]), Integer.parseInt(args[4]));
		}
		else {
			System.err.println("Usage: SBTDocumentModel <train|test> <args...>");
		}
	}

}
