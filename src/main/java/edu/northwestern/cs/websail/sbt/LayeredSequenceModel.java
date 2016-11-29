package edu.northwestern.cs.websail.sbt;
import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import junit.framework.Assert;

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
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class LayeredSequenceModel implements Serializable {

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
      double [] _res;
      
      public TestDoer(int i) {
        _doc = i;
      }
      
      public void run() {
        _res = testOnDocFullPpl(_doc);
        //_res = testOnDocFullPplExact(_doc, _wordTopicMarginal);
        System.out.println(_doc + "\t" + _res[0] + "\t" + _res[1]);
      }
    }
    
  //Model information:
  //state index zero is special start-of-sentence state
  int [] _branchingFactors; //dims: depth
  
  //hard-code these parameters:
  int _NUMLAYERS = 2; //number of layers of latent states in the sequence model
  boolean [][] _transitionTemplate = new boolean[][] {{true, false},{false, true}}; //dims: MUST BE _NUMLAYERS x _NUMLAYERS
  boolean [] _productionTemplate = new boolean[] {false, true}; //dims: MUST BE _NUMLAYERS
//  int _NUMLAYERS = 1; //number of layers of latent states in the sequence model
//  boolean [][] _transitionTemplate = new boolean[][] {{true}}; //dims: MUST BE _NUMLAYERS x _NUMLAYERS
//  boolean [] _productionTemplate = new boolean[] {true}; //dims: MUST BE _NUMLAYERS
  
  
  
  //assumed internal links are always from layer i to layer i+1
  
  SparseBackoffTree [][][] _forward; //dims: LAYERi-1 x LAYERi x state
  SparseBackoffTree []  _startState; //dims: LAYER
  SparseBackoffTree [][][] _backward; //dims: LAYERi+1 x LAYERi x state
  SparseBackoffTree [][] _wordToState; //dims: LAYER x state
  SparseBackoffTree [][] _down; //dims: LAYER-1 x state
  SparseBackoffTree [][] _up; //dims: LAYER-1 x state
  
  
  double [][] _wordStateMarginal; //dims: LAYER x state
  double [][] _startStateMarginal; //dims: LAYER x state
  double [][][] _backwardStateMarginal; //dims: LAYERi+1 x LAYERi x state
  double [][][] _forwardStateMarginal; //dims: LAYERi-1 x LAYERi x state
  double [][] _upStateMarginal; //dims: LAYER x state
  double [][] _downStateMarginal; //dims: LAYER x state
  
  
  SparseBackoffTreeStructure [] _struct; //dims: LAYER
  double [][] _forwardDelta = null; //dims: LAYER x depth
  double [][] _backwardDelta = null; //dims: LAYER x depth
  double [][] _wordDelta = null; //dims: LAYER x depth
  double [][] _upDelta = null; //dims: LAYER x depth
  double [][] _downDelta = null; //dims: LAYER x depth
    

  
  LayeredCorpus _c;
  
  Random _r = new Random();
  
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
  
  public LayeredSequenceModel(String configFile) throws Exception {
    _c = new LayeredCorpus(_NUMLAYERS);
    readConfigFile(configFile);
    if(_USEEXPANSION == 1) {
      if(_EXPANSIONSCHEDULE == null) {
        _EXPANSIONSCHEDULE = defaultExpansion(_branchingFactors.length);
      }
      _expansionBranchFactors = Arrays.copyOf(_branchingFactors, _branchingFactors.length);
      _branchingFactors = Arrays.copyOf(_branchingFactors, _EXPANSIONSCHEDULE[0]+1);
    }
    _struct = new SparseBackoffTreeStructure[_NUMLAYERS];
    for(int i=0; i<_struct.length; i++)
      _struct[i] = new SparseBackoffTreeStructure(_branchingFactors);
    initDeltas();
  }
  
  public void initDeltas() {
    _forwardDelta = new double[_NUMLAYERS][_branchingFactors.length];
    _backwardDelta = new double[_NUMLAYERS][_branchingFactors.length];
    _wordDelta =new double[_NUMLAYERS][_branchingFactors.length];
    _upDelta =new double[_NUMLAYERS][_branchingFactors.length];
    _downDelta =new double[_NUMLAYERS][_branchingFactors.length];
    double initDelta = 0.9 / (double)_branchingFactors.length;
    for(int i=0; i<_NUMLAYERS; i++) {
      Arrays.fill(_wordDelta[i], initDelta);
      Arrays.fill(_forwardDelta[i], initDelta);
      Arrays.fill(_backwardDelta[i], initDelta);
      Arrays.fill(_upDelta[i], initDelta);
      Arrays.fill(_downDelta[i], initDelta);
    }
  }
  
  public int get_NUMLAYERS() {
    return _NUMLAYERS;
  }
  
  public void set_NUMLAYERS(int _NUMLAYERS) {
    this._NUMLAYERS = _NUMLAYERS;
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
  
  //returns offset
  private int setThreadBreaks(int numPartitions, int partition) {
    
    _THREADBREAKS = new int[_NUMTHREADS];
    long approxToks = _c._NUMTOKENS / (_NUMTHREADS * numPartitions);
    long skipToks = partition * approxToks * _NUMTHREADS;
    long ct = 0;
    int thread = 0;
    int j=0;
    while(ct < skipToks)
      ct += _c._docs[j++].size();
    ct = 0;
    for(int i=j; (i<_c._NUMDOCS)&&thread<_NUMTHREADS; i++ ) {
      ct += _c._docs[i].size();
      if(ct > approxToks) {
        _THREADBREAKS[thread++] = i;
        ct = 0;
      }
    }
    //double-check that extra goes in last thread:
    if(partition ==numPartitions - 1)
      _THREADBREAKS[_NUMTHREADS - 1] = _c._NUMDOCS - 1;
    return j;
  }
  
  private void setThreadBreaks() {
    setThreadBreaks(1, 0);
//    _THREADBREAKS = new int[_NUMTHREADS];
//    long approxToks = _c._NUMTOKENS / _NUMTHREADS;
//    long ct = 0;
//    int thread = 0;
//    for(int i=0; i<_c._NUMDOCS; i++ ) {
//      ct += _c._docs[i].size();
//      if(ct > approxToks) {
//        _THREADBREAKS[thread++] = i;
//        ct = 0;
//      }
//    }
//    //extra goes in last thread:
//    _THREADBREAKS[_NUMTHREADS - 1] = _c._NUMDOCS - 1;
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
    TIntArrayList [] zs = _c._z[docId];
    TIntArrayList [] scratchZs = _c._scratchZ[docId];
    TIntArrayList doc = _c._docs[docId];
    for(int layer=0; layer<_NUMLAYERS; layer++) {
      for(int i=0; i<doc.size(); i++) {
        double chg = _c._changeFactor[docId][layer].get(i);
        if(_r.nextDouble() > chg) { //skip it
          scratchZs[layer].set(i, zs[layer].get(i));
          continue;       
        }
        int newZ = sampleZ(docId, layer, i, true, true);
        chg = _c.lambda * chg;
        if(newZ != zs[layer].get(i)) {
          changes++;
          chg = 1.0;
        }
        scratchZs[layer].set(i, newZ);
        _c._changeFactor[docId][layer].set(i, chg);
      }
    }
    return changes;
  }
  
  private static class MarkovBlanketContainer {
    SparseBackoffTree [] sbts;
    double [] subs;
  }
  
  //Returns the SBTs and amounts to subtract for the given latent variable
  public MarkovBlanketContainer getMarkovBlanket(int doc, int layer, int pos, boolean useWord, boolean sub) {
    int w = _c._docs[doc].get(pos);
    int curZ = _c._z[doc][layer].get(pos);
    ArrayList<SparseBackoffTree> sbts = new ArrayList<>();
    ArrayList<Double> doubs = new ArrayList<>();
    //forward:
    if(pos==0) {
      sbts.add(this._startState[layer]);
      if(sub)
        doubs.add(this._startStateMarginal[layer][curZ]);
    }
    else {
      for(int i=0; i<_NUMLAYERS; i++) {
        if(this._transitionTemplate[i][layer]) {
          if(_c._z[doc][i].get(pos-1) >= 0) {
            SparseBackoffTree sbt = this._forward[i][layer][_c._z[doc][i].get(pos-1)];
            sbts.add(sbt);
            if(sub)
              doubs.add(this._forwardStateMarginal[i][layer][curZ]);
          }
        }
      }
    }
    //backward:
    if(pos!=_c._docs[doc].size() - 1) {
      for(int i=0; i<_NUMLAYERS; i++) {
        if(this._transitionTemplate[layer][i]) {
          if(_c._z[doc][i].get(pos+1) >= 0) {
            SparseBackoffTree sbt = this._backward[i][layer][_c._z[doc][i].get(pos+1)];
            sbts.add(sbt);
            if(sub)
              doubs.add(this._backwardStateMarginal[i][layer][curZ]);
          }
        }
      }
    }
    if(useWord) {
      //word:
      if(_productionTemplate[layer]) {
        SparseBackoffTree sbtWord =this._wordToState[layer][w]; 
        if(sbtWord._totalMass > 0.0) {
          sbts.add(this._wordToState[layer][w]);
          if(sub)
            doubs.add(this._wordStateMarginal[layer][curZ]);
        }
      }
    }
      //layers:
    if(layer!=_NUMLAYERS-1) {
      if(_c._z[doc][layer+1].get(pos) >= 0) {
        sbts.add(this._up[layer+1][_c._z[doc][layer+1].get(pos)]);
        if(sub)
          doubs.add(this._upStateMarginal[layer+1][curZ]);
      }
    }
    if(layer!=0) {
      if(_c._z[doc][layer-1].get(pos) >= 0) {
        sbts.add(this._down[layer-1][_c._z[doc][layer-1].get(pos)]);
        if(sub)
          doubs.add(this._downStateMarginal[layer-1][curZ]);
      }
    }
    
    MarkovBlanketContainer out = new MarkovBlanketContainer();
    out.sbts = new SparseBackoffTree[sbts.size()];
    out.subs = new double[sbts.size()];
    for(int i=0; i<sbts.size(); i++) {
      out.sbts[i] = sbts.get(i);
      if(sub)
        out.subs[i] = doubs.get(i);
    }
    return out;
  }
  
  
  public int sampleZ(int doc, int layer, int pos, boolean useWord, boolean sub) {
    MarkovBlanketContainer mbc = getMarkovBlanket(doc, layer, pos, useWord, sub);
    SparseBackoffTreeIntersection sbti;
    int curZ = _c._z[doc][layer].get(pos);

    if(sub)
      sbti = new SparseBackoffTreeIntersection(mbc.sbts, curZ, mbc.subs, true);
    else
      sbti = new SparseBackoffTreeIntersection(mbc.sbts, true);
    
    int sample = sbti.sample(_r);
    return sample;
  }
  
  private long sampleRange(int start, int end) {
    
    long chg = 0;
    for(int i=start; i<=end; i++) {
      chg += (long)sampleDoc(i);
    }

    return chg;
  }
  
  //if scratch, copies z into scratch
  public void initZ(TIntArrayList [] ds, int [] maxVal, boolean scratch) {
    TIntArrayList [][] zs;
    if(scratch) {
      _c._scratchZ = new TIntArrayList[ds.length][_NUMLAYERS];
      zs = _c._scratchZ;
    }
    else {
      _c._z = new TIntArrayList[ds.length][_NUMLAYERS];
      zs = _c._z;
    }
    for(int i=0; i<ds.length; i++) {
      for(int j=0; j<_NUMLAYERS;j++) {
        zs[i][j] = new TIntArrayList(ds[i].size());
        for(int k=0; k<ds[i].size(); k++)
          zs[i][j].add((!scratch) ? _r.nextInt(maxVal[j]) : _c._z[i][j].get(k));
      }
    }
  }

  private static class AggregatedCounts {
    TIntDoubleHashMap [][][] forward; //LAYERi-1 x LAYERi x NUMSTATES
    TIntDoubleHashMap [][][] backward; //LAYERi+1 x LAYERi x NUMSTATES
    TIntDoubleHashMap [][] up; //LAYER x NUMSTATES...[i][] is counts *from* layer i
    TIntDoubleHashMap [][] down; //LAYER x NUMSTATES...[i][] is counts *from* layer i
    TIntDoubleHashMap [] start; //LAYER
    TIntDoubleHashMap [] end; //LAYER
    TIntDoubleHashMap [][] word; //LAYER x VOCAB
  }
  
  //doesn't create new hashmaps, only arrays
  private AggregatedCounts getInitializedAggregator() {
    AggregatedCounts ac = new AggregatedCounts();
    
    ac.forward = new TIntDoubleHashMap[_NUMLAYERS][_NUMLAYERS][];
    ac.backward = new TIntDoubleHashMap[_NUMLAYERS][_NUMLAYERS][];
    ac.up = new TIntDoubleHashMap[_NUMLAYERS][];
    ac.down = new TIntDoubleHashMap[_NUMLAYERS][];
    ac.start = new TIntDoubleHashMap[_NUMLAYERS];
    ac.end = new TIntDoubleHashMap[_NUMLAYERS];
    ac.word = new TIntDoubleHashMap[_NUMLAYERS][this.get_VOCABSIZE()];
    //init:
    for(int i=0; i<_NUMLAYERS; i++) {
      if(i!=_NUMLAYERS-1) {
        ac.down[i] = new TIntDoubleHashMap[this._struct[i].numLeaves()];
      }
      if(i > 0) {
        ac.up[i] = new TIntDoubleHashMap[this._struct[i].numLeaves()];
      }
      for(int j=0; j<_NUMLAYERS; j++) {
        ac.forward[i][j] = new TIntDoubleHashMap[this._struct[i].numLeaves()];
        ac.backward[i][j] = new TIntDoubleHashMap[this._struct[i].numLeaves()];
      }
    }
    return ac;
  }
  
  public void createAdjPut(TIntDoubleHashMap [] a, int idx, int j, double val) {
    if(a[idx]==null)
      a[idx] = new TIntDoubleHashMap();
    a[idx].adjustOrPutValue(j, val, val);
  }
  
  public AggregatedCounts aggregateCounts(TIntArrayList [][] zs, TIntArrayList [] ws) {
    AggregatedCounts ac = getInitializedAggregator();
    for(int doc=0; doc<zs.length; doc++) {
      TIntArrayList w= ws[doc];
      TIntArrayList [] z= zs[doc];
      //start
      for(int i=0; i<_NUMLAYERS; i++) {
        createAdjPut(ac.start, i, z[i].get(0), 1.0);
      }
      //transitions, layers, words:
      for(int i=0; i<w.size(); i++) {
        int word = w.get(i);
        if(i < w.size() - 1) {
          for(int j=0; j<_NUMLAYERS; j++) {
            for(int k=0; k<_NUMLAYERS; k++) {
              if(this._transitionTemplate[j][k]) {
                createAdjPut(ac.forward[j][k], z[j].get(i), z[k].get(i+1), 1.0); 
                createAdjPut(ac.backward[k][j], z[k].get(i+1), z[j].get(i), 1.0); 
              }
            }
          }
        }
        for(int j=0; j<_NUMLAYERS;j++) {
          if(j > 0) {
            createAdjPut(ac.up[j], z[j].get(i), z[j-1].get(i), 1.0);
          }
          if(j < _NUMLAYERS - 1) {
            createAdjPut(ac.down[j], z[j].get(i), z[j+1].get(i), 1.0);
          }
          if(this._productionTemplate[j]) {
            createAdjPut(ac.word[j], word, z[j].get(i), 1.0);
          }
        }
      }
      //end:
      for(int i=0; i<_NUMLAYERS; i++) {
        createAdjPut(ac.end, i, z[i].get(w.size()-1), 1.0);
      }
    }
    return ac;
  }
  
  
  public SparseBackoffTree [] getParamsFromHash(double [] ds, TIntDoubleHashMap [] hm, SparseBackoffTreeStructure struct) {
    SparseBackoffTree [] out = new SparseBackoffTree[hm.length];
    
    for(int i=0; i<hm.length; i++) {
      out[i] = new SparseBackoffTree(struct);
    }
    
    //add:
    for(int i=0; i<hm.length; i++) {
      if(hm[i] == null)
        continue;
      out[i] = getParamsFromHash(ds, hm[i], struct);
    }
    return out;   
  }
  
  public SparseBackoffTree getParamsFromHash(double [] ds, TIntDoubleHashMap hm, SparseBackoffTreeStructure struct) {
    SparseBackoffTree out = new SparseBackoffTree(struct);
    TIntDoubleIterator it = hm.iterator();
    while(it.hasNext()) {
      it.advance();
      int z = it.key();
      double val = it.value();
      out.smoothAndAddMass(z, val, ds);
    }
    return out;
  }
  
  public TIntDoubleHashMap expandHash(TIntDoubleHashMap hm, int multiplier) {
    TIntDoubleHashMap out = new TIntDoubleHashMap();
    TIntDoubleIterator it = out.iterator();
    while(it.hasNext()) {
      it.advance();
      int base = it.key()*multiplier;
      for(int i=base; i<base+multiplier; i++) {
        out.put(i, Math.max(0.0, it.value()/(double)multiplier));
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
  public int initializeForCorpus(String inFile, int [] maxVal) throws Exception {
    int toks = _c.readCorpusDat(inFile, true);
    setThreadBreaks();
    initZ(_c._docs, maxVal, false);
    updateModel(_c._z);
    return toks;
  }
  
  private int gibbsPass() {
    return gibbsPass(0);
  }
  
    private int gibbsPass(int offset) { 
      int changes= 0 ;
      GibbsDoer [] gds = new GibbsDoer[_NUMTHREADS];
      long stTime = System.currentTimeMillis();
        ExecutorService e = Executors.newFixedThreadPool(_NUMTHREADS);
        for(int i=0; i<_NUMTHREADS;i++) {
          gds[i] = new GibbsDoer();
          gds[i]._start = offset;
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
  
  private static double [] ones(int len) {
    double [] out = new double[len];
    Arrays.fill(out, 1.0);
    return out;
  }

  
  /**
   * Updates the model given the topic assignments (_z) and divides by marginal for next sampling pass
   */
    public void updateModel(TIntArrayList [][] zs) {
      System.out.println("arch: " + Arrays.toString(this._branchingFactors));
      System.out.println("second doc samples (layer 0): " + zs[1][0].toString());
      if(_NUMLAYERS > 1)
        System.out.println("second doc samples (layer 1): " + zs[1][1].toString());
      if(_c._changeFactor != null)
        System.out.println("second doc chg (layer 0): " + _c._changeFactor[1][0].toString());
      AggregatedCounts ac = this.aggregateCounts(zs, _c._docs);
      
      this._startState  = new SparseBackoffTree[_NUMLAYERS];
      this._wordToState  = new SparseBackoffTree[_NUMLAYERS][];
      this._forward  = new SparseBackoffTree[_NUMLAYERS][_NUMLAYERS][];
      this._backward = new SparseBackoffTree[_NUMLAYERS][_NUMLAYERS][];
      this._up  = new SparseBackoffTree[_NUMLAYERS][];
      this._down  = new SparseBackoffTree[_NUMLAYERS][];
      this._startStateMarginal = new double[_NUMLAYERS][];
      this._backwardStateMarginal = new double[_NUMLAYERS][_NUMLAYERS][];
      this._forwardStateMarginal =  new double[_NUMLAYERS][_NUMLAYERS][];
      this._downStateMarginal = new double[_NUMLAYERS][];
      this._upStateMarginal = new double[_NUMLAYERS][];
      this._wordStateMarginal = new double[_NUMLAYERS][];
      
      for(int i=0; i<_NUMLAYERS;i++) {
        this._startState[i] = getParamsFromHash(_forwardDelta[i], ac.start[i], this._struct[i]);
        for(int j=0; j<_NUMLAYERS;j++) {
          this._forward[i][j] = getParamsFromHash(_forwardDelta[j], ac.forward[i][j], this._struct[j]);
          this._backward[j][i] = getParamsFromHash(_backwardDelta[i], ac.backward[j][i], this._struct[i]);
        }
        this._wordToState[i] = getParamsFromHash(_wordDelta[i], ac.word[i], this._struct[i]);
        if(i > 0)
          this._down[i-1] = getParamsFromHash(_downDelta[i-1], ac.down[i-1], this._struct[i]);
        if(i < _NUMLAYERS - 1)
          this._up[i+1] = getParamsFromHash(_upDelta[i+1], ac.up[i+1], this._struct[i]);
      }

      for(int j=0; j<_NUMLAYERS;j++) {
        this._startStateMarginal[j] = ones(_struct[j].numLeaves());
        for(int k=0; k<_NUMLAYERS;k++) {
          boolean first = true;
          if(_transitionTemplate[j][k]) {
            //first forward always has normalizers of 1.0:
            if(first) {
              first=false;
              this._forwardStateMarginal[j][k] = ones(_struct[k].numLeaves());
            }
            else {
              this._forwardStateMarginal[j][k] = getNormalizers(_forward[j][k], _struct[k]);
            }
            this._backwardStateMarginal[j][k] = getNormalizers(_backward[j][k], _struct[k]);
            
            for(int i=0; i<_backward[j][k].length; i++) {
              _backward[j][k][i].divideCountsBy(_backwardStateMarginal[j][k]);
              _forward[j][k][i].divideCountsBy(_forwardStateMarginal[j][k]);
            }
          }
        }
        if(j < _NUMLAYERS - 1) {
          this._downStateMarginal[j] = getNormalizers(_down[j], _struct[j+1]);
          for(int i=0; i<_down[j].length; i++) {
            _down[j][i].divideCountsBy(_downStateMarginal[j]);
          }
        }
        if(j > 0) {
          this._upStateMarginal[j] = getNormalizers(_up[j], _struct[j-1]);
          for(int i=0; i<_up[j].length; i++) {
            _up[j][i].divideCountsBy(_upStateMarginal[j]);
          }
        }
        if(this._productionTemplate[j]) {
          this._wordStateMarginal[j] = getNormalizers(_wordToState[j], _struct[j]);
          for(int i=0; i<_wordToState[j].length; i++) {
            _wordToState[j][i].divideCountsBy(_wordStateMarginal[j]); 
          }
        }
      }
      System.out.println("\tdiscounts forward: " + Arrays.toString(_forwardDelta[0]) + 
          "\tbackward " + Arrays.toString(_backwardDelta[0]) +
              "\tword " + Arrays.toString(_wordDelta[0]));
    }
    
    /**
     * Adds the gradient of the log likelihood for the given observations
     * Assumes sum of deltas < 1
     * @param grad  the gradient is added here
     * @param hm  the observations
     * @param curDeltas the current hyperparameters
     * @param incremental if true, measure log likelihood when introducing observations one-by-one in random order;
     * otherwise, use leave-one-out cross validation
     * @return  log likelihood using current parameters
     */
    public double updateGradient(double [] grad, TIntDoubleHashMap hm, double [] curDeltas, boolean incremental,
        SparseBackoffTreeStructure struct) {
      SparseBackoffTree sbt = new SparseBackoffTree(struct);
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
          int [] localIdxTrace = struct.getLocalIdxTrace(z);
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
          int [] localIdxTrace = struct.getLocalIdxTrace(z);
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
     * @param gradient  gradient to normalize
     * @param curDeltas current hyperparameters to which gradient will be applied
     * @param stepSize  how far to move (in L1)
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
     * @param deltas  hyperparameters to start from
     * @param hm  observation counts
     * @param stepSize  how far to move in the step, in L1 norm
     * @param incremental if true, measure log likelihood when introducing observations one-by-one in random order;
     * otherwise, use leave-one-out cross validation
     * @return
     */
    public double [] gradientStep(double [] deltas, TIntDoubleHashMap [] hm, double stepSize, boolean incremental,
        SparseBackoffTreeStructure struct) {
      double [] gradient = new double[deltas.length];
      double LL = 0.0;
      for(int i=0; i<hm.length; i++) {
        if(hm[i] != null)
          LL += updateGradient(gradient, hm[i], deltas, incremental, struct);
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
    public void optimizeParameters(TIntArrayList [][] zs) {
      //TODO: implement
      return;
//      for(int i=0; i<_wordsDelta.length; i++) {
//        _wordsDelta[i] *= 0.9;
//        _docsDelta[i] *= 0.9;
//      }
      //TODO: multi-thread
//      double STEPSIZE = 0.01; //start stepping this far in L1 norm
//      double STEPDEC = 0.95; //decrease step size this much each step
//      int NUMSTEPS = 20; //number of steps to take
//      double STEPSIZE = 0.02; //start stepping this far in L1 norm
//      double STEPDEC = 0.8; //decrease step size this much each step
//      int NUMSTEPS = 10; //number of steps to take
//
//      TIntDoubleHashMap [] hm = aggregateCounts(0, zs, _c._docs);
//      System.out.println("words:");
//      double step = STEPSIZE;
//      for(int i=0; i<NUMSTEPS; i++) {
//          gradientStep(_wordDelta, hm, step, false);
//          step *= STEPDEC;
//      }
//      long totalparams = 0L;
//      for(TIntDoubleHashMap i : hm) {
//        if(i != null)
//          totalparams += i.size();
//      }
//      hm = aggregateCounts(-1, zs, _c._docs);
//      TIntDoubleHashMap [] hm2 = Arrays.copyOf(hm, hm.length + 1);
//      hm2[hm.length] = aggregateCounts(-2, zs, _c._docs)[0];
//      hm = hm2;
//      long transParams = 0L;
//      for(TIntDoubleHashMap i : hm) {
//        if(i != null)
//          transParams += i.size();
//      }
//      System.out.println("forward:");
//      step = STEPSIZE;
//      for(int i=0; i<NUMSTEPS; i++) {
//          gradientStep(_forwardDelta, hm, step, false);
//          step *= STEPDEC;
//      }
//      hm = aggregateCounts(1, zs, _c._docs);
//      hm2 = Arrays.copyOf(hm, hm.length + 1);
//      hm2[hm.length] = aggregateCounts(2, zs, _c._docs)[0];
//      hm = hm2;
//      System.out.println("backward:");
//      step = STEPSIZE;
//      for(int i=0; i<NUMSTEPS; i++) {
//          gradientStep(_backwardDelta, hm, step, false);
//          step *= STEPDEC;
//      }
//      System.out.println("full word deltas: " + Arrays.toString(_wordDelta));
//      System.out.println("full forward deltas: " + Arrays.toString(_forwardDelta));
//      System.out.println("full backward deltas: " + Arrays.toString(_backwardDelta));
//      System.out.println("total nonzero word-state params: " + totalparams);
//      System.out.println("total forward trans params: " + transParams);
    }
    
    /**
     * Expands the model
     * @param expansionIndex  The desired index in _EXPANSIONSCHEDULE to expand to.
     * @return  whether the model was expanded
     */
    public boolean expandModel(int expansionIndex) {
      if(expansionIndex < _EXPANSIONSCHEDULE.length) {
        for(int layer=0;layer<_NUMLAYERS;layer++) {
          int curLeaves = _struct[layer].numLeaves(); 
          //get new structure
          _branchingFactors = Arrays.copyOf(this._expansionBranchFactors, _EXPANSIONSCHEDULE[expansionIndex] + 1);
          _struct[layer] = new SparseBackoffTreeStructure(_branchingFactors);
          initDeltas();
          int multiplier = _struct[layer].numLeaves()/curLeaves;
          System.out.println("Expanding to " + Arrays.toString(_branchingFactors));
          for(int i=0; i<_c._z.length;i++) {
            for(int j=0; j<_c._z[i][layer].size(); j++) {
              int z = multiplier * _c._z[i][layer].get(j) + _r.nextInt(multiplier);
              _c._z[i][layer].set(j,  z);
              _c._changeFactor[i][layer].set(j, 1.0);
            }
          }
        }
        updateModel(_c._z);
        return true;
      }
      else
        return false;
    }
    
    public void trainModel(int iterations, int updateInterval, String outFile) throws Exception {
      trainModel(iterations, updateInterval, outFile, false);
    }
    
    public void translateStates(TIntArrayList [] zs, TIntIntHashMap trans) {
      for(int i=0; i<zs.length; i++) { 
        for(int j=0; j<zs[i].size(); j++) {
          zs[i].set(j, trans.get(zs[i].get(j)));
        }
      }
    }
    
    
    /**
     * Trains the model for a specified number of iterations.
     *
     * @param  iterations  the number of Gibbs passes to perform
     * @param  updateInterval number of iterations between hyperparameter updates
     * @param  outFile   writes model to here before each expansion
     */
  public void trainModel(int iterations, int updateInterval, String outFile, boolean expandOnEntry) throws Exception {
    boolean done = false;
    int j=0;
    int [] maxVal = new int[_NUMLAYERS];
    for(int i=0; i<_NUMLAYERS; i++) {
      maxVal[i] = _struct[i].numLeaves();
    }
    if(expandOnEntry) {
      done = !expandModel(j);
    }
    boolean PARTIALPASS = true;
    while(!done) {
      for(int i=1; i<=iterations; i++) {
        System.out.print(i + ": ");
        int offset = 0;
        if(PARTIALPASS) {
          offset = setThreadBreaks(10, (i-1) % 10);
          System.out.println(offset + "\t" + Arrays.toString(_THREADBREAKS));
        }
        initZ(_c._docs, maxVal, true); //init scratch array to zs
        gibbsPass(offset);
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
    SparseBackoffTree [][] sbtW = this._wordToState;
    SparseBackoffTree [][][] sbtF = this._forward;
    SparseBackoffTree [][][] sbtB = this._backward;
    SparseBackoffTree [] sbtS = this._startState;
    SparseBackoffTree [][] sbtU = this._up;
    SparseBackoffTree [][] sbtD = this._down;
    SparseBackoffTreeStructure [] struct = this._struct;
    _struct = null;
    _wordToState = null;
    _forward = null;
    _backward = null;
    _startState = null;
    _up = null;
    _down = null;
    ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outFile));
    oos.writeObject(this);
    oos.close();
    _wordToState = sbtW;
    _forward = sbtF;
    _backward = sbtB;
    _struct = struct;
    _startState = sbtS;
    _up = sbtU;
    _down = sbtD;
  }
  
  public static LayeredSequenceModel readModel(String inFile) throws Exception {
    ObjectInputStream ois = new ObjectInputStream(new FileInputStream(inFile));
    LayeredSequenceModel out = (LayeredSequenceModel) ois.readObject();
    ois.close();
    out._struct = new SparseBackoffTreeStructure[out._NUMLAYERS];
    for(int i=0; i<out._NUMLAYERS; i++)
      out._struct[i] = new SparseBackoffTreeStructure(out._branchingFactors);
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
//        for(int j=0; j<dist[i].length;j++) {
//          sbtForward.addWeighted(this._forward[j], dist[i][j]);
//        }
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
  
  
  //tests on a single document using left-to-right method
  //word topic marginal (i.e. P(topic) computed from P(topic, word)) is supplied to ensure evaluated distribution sums to one
  public double [] testOnDocFullPpl(int doc) {
    double LL = 0.0;
    double [] ps = new double[_c._docs[doc].size()];
    for(int j=0; j<_NUMPARTICLES; j++) {
      TIntArrayList words = _c._docs[doc];
      for(int layer=0;layer<_NUMLAYERS;layer++) {
        for(int i=0; i<words.size(); i++) {
          _c._z[doc][layer].set(i, -1);
        }
      }
      for(int i=0; i<words.size(); i++) {
        int w = words.get(i);
        double p = 0.0;           
        
        //sample this state forward:
        for(int layer=0;layer<_NUMLAYERS;layer++) {
          _c._z[doc][layer].set(i, this.sampleZ(doc, layer, i, false, false));
        }
        
        if(_c._pWord[w]>0.0) {
          
          //test on this word:
          double [] pOut = new double[this.get_VOCABSIZE()];
          double normalizer = 0.0;
          for(int k=0; k<pOut.length; k++) {
            if(_c._pWord[k] >0.0) {
              pOut[k] = 1.0;
              for(int layer=0; layer<_NUMLAYERS; layer++) {
                if(this._productionTemplate[layer])
                  pOut[k] *= this._wordToState[layer][k].getSmoothed(_c._z[doc][layer].get(i));
              }
              normalizer += pOut[k];
            }
          }
          for(int k=0; k<pOut.length; k++) {
            pOut[k] /= normalizer;
          }
          p = pOut[w];
        }
        ps[i] += p;
        //re-sample using this word:
        //backward:
//        for(int idx=i-1; idx>=0; idx--) {
//          for(int layer=0;layer<_NUMLAYERS;layer++) {
//            _c._z[doc][layer].set(idx, this.sampleZ(doc, layer, idx, true, false));
//          }
//        }
//        //forward:
//        for(int idx=0; idx<i; idx++) {
//          for(int layer=0;layer<_NUMLAYERS;layer++) {
//            _c._z[doc][layer].set(idx, this.sampleZ(doc, layer, idx, true, false));
//          }
//        }
        //current:
        for(int pass=0; pass<3; pass++) {
          for(int layer=0;layer<_NUMLAYERS;layer++) {
            _c._z[doc][layer].set(i, this.sampleZ(doc, layer, i, true, false));
          }
        }
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
  
  public static void train(String inputFile, String outputFile, String configFile) throws Exception {
    LayeredSequenceModel sbtsm = new LayeredSequenceModel(configFile);
    int [] maxVal = new int[sbtsm._NUMLAYERS];
    for(int i=0; i<maxVal.length; i++) {
      maxVal[i] = sbtsm._struct[i].numLeaves();
    }
    sbtsm.initializeForCorpus(inputFile, maxVal);
    sbtsm.trainModel(sbtsm._NUMITERATIONS, sbtsm._OPTIMIZEINTERVAL, outputFile);
    sbtsm.writeModel(outputFile);
  }
  
    
  
  //computes log likelihood of model using left-to-right method
  //returns {ppl, number of tested words}
  //numDocs is number in test file
  //maxDocs is number actually tested (starting from the beginning of the file)
  public static double [] testModel(String modelFile, String testFile, int numDocs, int maxDocs, String configFile) throws Exception {
    
    LayeredSequenceModel sbtsm = readModel(modelFile);
    sbtsm.readConfigFile(configFile);

    sbtsm._c.reInitForCorpus(testFile, numDocs);
    sbtsm.setThreadBreaks();
    double LL = 0.0;
    double numWords = 0.0;
    TestDoer [] tds = new TestDoer[maxDocs];
    ExecutorService e = Executors.newFixedThreadPool(sbtsm._NUMTHREADS);
    for(int i=0; i<tds.length;i++) {
      if(sbtsm._MAXTESTDOCSIZE < 0 || sbtsm._c._docs[i].size() > sbtsm._MAXTESTDOCSIZE) {
        System.out.println("skipping " + i + " size " + sbtsm._c._docs[i].size());
        continue;
      }
      tds[i] = sbtsm.new TestDoer(i);
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
        return;
      }
      test(args[1], args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3]);
    }
    else {
      System.err.println("Usage: <train|test> <args...>");
    }
  }
  
}
