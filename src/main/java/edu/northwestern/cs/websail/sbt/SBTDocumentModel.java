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
		TIntArrayList doc = _docs[docId];
		TIntDoubleHashMap docCts = aggregateCounts(doc);
		SparseBackoffTree sbtDoc = new SparseBackoffTree(_struct);
		sbtDoc.addAllMass(docCts);
		for(int i=0; i<doc.size(); i++) {
			int newZ = sampleZ(docId, i, false, sbtDoc);
			if(newZ != doc.get(i))
				changes++;
			doc.set(i, newZ);
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
	public int readCorpusDat(String inFile) throws Exception {
		BufferedReader brIn = new BufferedReader(
				new InputStreamReader(new FileInputStream(inFile), "UTF8" ));
		String sLine;
		int i=0;
		int toks = 0;
		ArrayList<TIntArrayList> aldocs = new ArrayList<TIntArrayList>();
		while((sLine = brIn.readLine())!=null) {
			TIntArrayList ll = lineToList(sLine);
			aldocs.add(ll);
			toks += ll.size();
		}
		brIn.close();
		System.out.println("Tokens: " + toks);
		_NUMTOKENS = toks;
		_NUMDOCS = aldocs.size();
		_docs = aldocs.toArray(new TIntArrayList[aldocs.size()]);
		setThreadBreaks();
		return i;
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
		int toks = readCorpusDat(inFile);
		initZRandom(_docs, maxVal);
		updateWordParamsFromZs(interpolateEndPoints(this._wordsDeltaEndPts, _branchingFactors.length));
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
    	System.out.print("\ttime: " + stTime + "\tchanges: " + changes);
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
			smoothing[i] = smoothAndCount[i];
			if(smoothing[i] > maxSmooth) {
				maxSmooth = smoothing[i];
			}
			count[i] = smoothAndCount[i];
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
	
    public void updateModel() {
		
    	this.updateWordParamsFromZs(this.interpolateEndPoints(_wordsDeltaEndPts, _branchingFactors.length));
    	
		_topicMarginal = getNormalizers(_topicGivenWord, _struct);
		for(int i=0; i<_topicGivenWord.length; i++) 
			_topicGivenWord[i].divideCountsBy(_topicMarginal);
		System.out.println("\tdiscounts doc: " + Arrays.toString(this._docsDeltaEndPts) + "\tword: " + Arrays.toString(this._wordsDeltaEndPts));
    }
    
    public void optimizeParameters() {
    	System.out.println("here is where we would optimize parameters.");
    }
    
	public void trainModel(int iterations, int updateInterval) {
		
		for(int i=1; i<=iterations; i++) {
			gibbsPass();
			if((iterations % updateInterval)==0) { 
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
		out.updateWordParamsFromZs(interpolateEndPoints(out._wordsDeltaEndPts, out._branchingFactors.length));
		return out;
	}
	
	public static void train(String inputFile, String outputFile, String configFile) throws Exception {
		SBTDocumentModel sbtdm = new SBTDocumentModel(configFile);
		sbtdm.initializeForCorpus(inputFile, sbtdm._struct.numLeaves());
		sbtdm.trainModel(sbtdm._NUMITERATIONS, sbtdm._OPTIMIZEINTERVAL);
		sbtdm.writeModel(outputFile);
	}
	
	public static double test(String modelFile, String inputFile) {
		return -1.0;
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
			if(args.length != 3) {
				System.err.println("Usage: SBTDocumentModel test <model_file> <test_file>");
			}
		}
		else {
			System.err.println("Usage: SBTDocumentModel <train|test> <args...>");
		}
	}

}
