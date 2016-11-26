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

public class LayeredSequenceModel {

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
    
  //Model information:
  //state index zero is special start-of-sentence state
  int [] _branchingFactors; //dims: depth
  
  //hard-code these parameters:
  int _NUMLAYERS = 2; //number of layers of latent states in the sequence model
  boolean [][] _transitionTemplate = new boolean[][] {{true, true},{true, true}}; //dims: MUST BE _NUMLAYERS x _NUMLAYERS
  boolean [] _productionTemplate = new boolean[] {false, true}; //dims: MUST BE _NUMLAYERS
  //assumed internal links are always from layer i to layer i+1
  
  SparseBackoffTree [][][] _forward; //dims: LAYERi x LAYERi-1 x state
  SparseBackoffTree []  _startState; //dims: LAYER
  SparseBackoffTree [] _endState; //dims: LAYER
  SparseBackoffTree [][][] _backward; //dims: LAYERi x LAYERi+1 x state
  SparseBackoffTree [][] _wordToState; //dims: LAYER x state
  SparseBackoffTree [] _down; //dims: LAYER-1
  SparseBackoffTree [] _up; //dims: LAYER-1
  
  
  double [][] _wordStateMarginal; //dims: LAYER x state
  double [][] _startStateMarginal; //dims: LAYER x state
  double [][] _endStateMarginal; //dims: LAYER x state
  double [][][] _backwardStateMarginal; //dims: LAYERi x LAYERi+1 x state
  double [][][] _forwardStateMarginal; //dims: LAYERi x LAYERi-1 x state
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
    _c = new LayeredCorpus();
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
    double initDelta = 0.9 / (double)_branchingFactors.length;
    for(int i=0; i<_NUMLAYERS; i++) {
      Arrays.fill(_wordDelta, initDelta);
      Arrays.fill(_forwardDelta, initDelta);
      Arrays.fill(_backwardDelta, initDelta);
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
        double chg = _c._changeFactor[layer][docId].get(i);
        if(_r.nextDouble() > chg) { //skip it
          scratchZs[layer].set(i, zs[layer].get(i));
          continue;       
        }
        int newZ = sampleZ(docId, layer, i, false);
        chg = _c.lambda * chg;
        if(newZ != zs[layer].get(i)) {
          changes++;
          chg = 1.0;
        }
        scratchZs[layer].set(i, newZ);
        _c._changeFactor[layer][docId].set(i, chg);
      }
    }
    return changes;
  }
  
  private static class MarkovBlanketContainer {
    SparseBackoffTree [] sbts;
    double [] subs;
  }
  
  //Returns the SBTs and amounts to subtract for the given latent variable
  public MarkovBlanketContainer getMarkovBlanket(int doc, int layer, int pos, boolean testing) {
    int w = _c._docs[doc].get(pos);
    int curZ = _c._z[doc][layer].get(pos);
    ArrayList<SparseBackoffTree> sbts = new ArrayList<>();
    ArrayList<Double> doubs = new ArrayList<>();
    //forward:
    if(pos==0) {
      sbts.add(this._startState[layer]);
      doubs.add(this._startStateMarginal[layer][curZ]);
    }
    else {
      for(int i=0; i<_NUMLAYERS; i++) {
        if(this._transitionTemplate[i][layer]) {
          SparseBackoffTree sbt = this._forward[layer][i][_c._z[doc][i].get(pos-1)];
          sbts.add(sbt);
          doubs.add(this._forwardStateMarginal[layer][i][curZ]);
        }
      }
    }
    //backward:
    if(!testing) {
      if(pos==_c._docs[doc].size() - 1) {
        sbts.add(this._endState[layer]);
        doubs.add(this._endStateMarginal[layer][curZ]);
      }
      else {
        for(int i=0; i<_NUMLAYERS; i++) {
          if(this._transitionTemplate[layer][i]) {
            SparseBackoffTree sbt = this._backward[layer][i][_c._z[doc][i].get(pos+1)];
            sbts.add(sbt);
            doubs.add(this._backwardStateMarginal[layer][i][curZ]);
          }
        }
      }
      //word:
      if(_productionTemplate[layer]) {
        SparseBackoffTree sbtWord =this._wordToState[layer][w]; 
        if(sbtWord._totalMass > 0.0) {
          sbts.add(this._wordToState[layer][w]);
          doubs.add(this._endStateMarginal[layer][curZ]);
        }
      }
    }
    //layers:
    if(!testing && layer!=_NUMLAYERS-1) {
      sbts.add(this._up[_c._z[doc][layer+1].get(pos)]);
      doubs.add(this._upStateMarginal[layer][curZ]);
    }
    if(layer!=0) {
      sbts.add(this._down[_c._z[doc][layer-1].get(pos)]);
      doubs.add(this._downStateMarginal[layer][curZ]);
    }
    
    MarkovBlanketContainer out = new MarkovBlanketContainer();
    out.sbts = (SparseBackoffTree []) sbts.toArray();
    out.subs = new double[sbts.size()];
    for(int i=0; i<sbts.size(); i++) {
      out.subs[i] = doubs.get(i);
    }
    return out;
  }
  
  
  public int sampleZ(int doc, int layer, int pos, boolean testing) {
    MarkovBlanketContainer mbc = getMarkovBlanket(doc, layer, pos, testing);
    SparseBackoffTreeIntersection sbti;
    int curZ = _c._z[layer][doc].get(pos);

    if(!testing)
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
  public void initZ(TIntArrayList [] ds, int maxVal, boolean scratch) {
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
          zs[i][j].add((!scratch) ? _r.nextInt(maxVal) : _c._z[i][j].get(k));
      }
    }
  }

  private static class AggregatedCounts {
    TIntDoubleHashMap [][][] forward; //LAYERi-1 x LAYERi x NUMSTATES
    TIntDoubleHashMap [][][] backward; //LAYERi+1 x LAYERi x NUMSTATES
    TIntDoubleHashMap [][] up; //LAYER x NUMSTATES
    TIntDoubleHashMap [][] down; //LAYER x NUMSTATES
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
        ac.up[i] = new TIntDoubleHashMap[this._struct[i].numLeaves()];
        ac.down[i] = new TIntDoubleHashMap[this._struct[i].numLeaves()];
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
        if(i < w.size()) {
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
            createAdjPut(ac.down[j+1], z[j].get(i), z[j+1].get(i), 1.0);
          }
          if(this._productionTemplate[j]) {
            createAdjPut(ac.word[j], z[j].get(i), word, 1.0);
          }
        }
      }
      //end:
      
    }
    return ac;
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
   * @param dsWords The discounts to use in the returned parameters
   * @param from  see above
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
  
  //TODO: remove redundancy with code directly above.
  public SparseBackoffTree [] getParamsFromHash(double [] dsWords, int from, TIntDoubleHashMap [] hm) {
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
  
  public TIntDoubleHashMap [] accToHash(int from) {
    TIntDoubleHashMap [] out;
    if(from==-2) { //start
      out = new TIntDoubleHashMap[1];
      out[0] = new TIntDoubleHashMap();
      for(int i=0; i<this._startAcc.length;i++) {
        out[0].put(i, _startAcc[i]);
      }
    }
    else if(from ==-1) { //forward
      out = new TIntDoubleHashMap[_forwardArrays.length];
      for(int i=0; i<out.length;i++) {
        out[i] = new TIntDoubleHashMap();
      }
      for(int j=0; j<this._transAcc.length;j++) {
        for(int k=0; k<this._transAcc[j].length; k++) {
            out[j].put(k, _transAcc[j][k]);
        }
      }
    }

    else if(from==0) { //word
      out = new TIntDoubleHashMap[_wordAcc.length];
      for(int i=0; i<out.length; i++) {
        out[i] = new TIntDoubleHashMap();
        for(int j=0; j<_wordAcc[i].length; j++) {
          out[i].put(j, _wordAcc[i][j]);
        }
      }
    }
     else if(from==1) { //backward
        out = new TIntDoubleHashMap[_backwardArrays.length];
        for(int i=0; i<out.length;i++) {
          out[i] = new TIntDoubleHashMap();
        }
        for(int j=0; j<this._transAcc.length;j++) {
          for(int k=0; k<this._transAcc[j].length; k++) {
              out[k].put(j, _transAcc[j][k]);
          }
        }
      }
     else { //end state
      out = new TIntDoubleHashMap[1];
      out[0] = new TIntDoubleHashMap();
      for(int i=0; i<this._endAcc.length;i++) {
        out[0].put(i, _endAcc[i]);
      }
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
    
  public SparseBackoffTree [] getParamsFromAccs(double [] ds, int from) {
    TIntDoubleHashMap [] hm = accToHash(from);
    return getParamsFromHash(ds, from, hm);
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
  
    //TODO: remove redundancy with gibbsPass
    private double emPass() { 
      double loglikelihood= 0 ;
      EMDoer [] eds = new EMDoer[_NUMTHREADS];
      long stTime = System.currentTimeMillis();
      ExecutorService e = Executors.newFixedThreadPool(_NUMTHREADS);
      int numStates = this._forwardArrays.length;
      int numWords = this.get_VOCABSIZE();
      _startAcc = new double [numStates];
      _transAcc = new double[numStates][numStates];
      _endAcc = new double[numStates];
      _wordAcc = new double[numWords][numStates];
      for(int i=0; i<_NUMTHREADS;i++) {
        eds[i] = new EMDoer();
        eds[i]._start = 0;
        if(i > 0)
          eds[i]._start = _THREADBREAKS[i-1] + 1;
        eds[i]._end = _THREADBREAKS[i];
        eds[i]._numStates = numStates;
        eds[i]._numWords = numWords;
        e.execute(eds[i]);
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
        loglikelihood += eds[i]._loglikelihood;
        reduce(eds[i]);
      }
      stTime = System.currentTimeMillis() - stTime;
      System.out.println("\ttime: " + stTime + "\tLL: " + loglikelihood + "\tskipWords: " + _numSkipWords);
      //System.out.println("state marginal: " + Arrays.toString(this._wordStateMarginal));
      _numSkipWords = 0;
      //System.out.println(Arrays.toString(_topicMarginal));
      return loglikelihood;
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
  
  public void reduce(EMDoer ed) {
    for(int i=0; i<_startAcc.length;i++) {
      _startAcc[i] += ed._startAcc[i];
      _endAcc[i] += ed._endAcc[i];
    }
    for(int i=0; i<_wordAcc.length;i++) {
      for(int j=0; j<_wordAcc[i].length; j++) {
        _wordAcc[i][j] += ed._wordAcc[i][j];
      }
    }
    for(int i=0; i<_transAcc.length;i++) {
      for(int j=0; j<_transAcc[i].length; j++) {
        _transAcc[i][j] += ed._transAcc[i][j];
      }
    }   
  }
  
  public void updateModelFromEM() {
    
    System.out.println("arch: " + Arrays.toString(this._branchingFactors));
    this._startState = getParamsFromAccs(_forwardDelta, -2)[0];
    this._forward = getParamsFromAccs(_forwardDelta, -1);
    this._backward = getParamsFromAccs(_backwardDelta, 1);
    this._endState = getParamsFromAccs(_backwardDelta, 2)[0];
    this._wordToState = getParamsFromAccs(_wordDelta, 0);
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
   * Updates the model given the topic assignments (_z) and divides by marginal for next sampling pass
   */
    public void updateModel(TIntArrayList [] zs) {
      System.out.println("arch: " + Arrays.toString(this._branchingFactors));
    System.out.println("second doc samples: " + zs[1].toString());
    if(_c._changeFactor != null)
      System.out.println("second doc chg: " + _c._changeFactor[1].toString());
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
     * @param grad  the gradient is added here
     * @param hm  the observations
     * @param curDeltas the current hyperparameters
     * @param incremental if true, measure log likelihood when introducing observations one-by-one in random order;
     * otherwise, use leave-one-out cross validation
     * @return  log likelihood using current parameters
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
//      for(int i=0; i<_wordsDelta.length; i++) {
//        _wordsDelta[i] *= 0.9;
//        _docsDelta[i] *= 0.9;
//      }
      //TODO: multi-thread
//      double STEPSIZE = 0.01; //start stepping this far in L1 norm
//      double STEPDEC = 0.95; //decrease step size this much each step
//      int NUMSTEPS = 20; //number of steps to take
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
     * @param expansionIndex  The desired index in _EXPANSIONSCHEDULE to expand to.
     * @return  whether the model was expanded
     */
    public boolean expandModel(int expansionIndex) {
      if(expansionIndex < _EXPANSIONSCHEDULE.length) {
        int curLeaves = _struct.numLeaves(); 
        if(!GIBBS) {
          System.out.println("executing special expansion gibbs for EM");
          for(int i=0; i<0; i++) {
            initZ(_c._docs, curLeaves, true); //init scratch array to zs
            gibbsPass();
            _c._z = _c._scratchZ;
          }
        }
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
            _c._changeFactor[i].set(j, 1.0);
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
    
    public void renumberStatesForCompression(TIntArrayList [] zs) { //changes states in zs to make more frequent states have lower ids
      TIntIntHashMap stateCount = new TIntIntHashMap();
      TIntIntHashMap stateTranslation = new TIntIntHashMap();
      for(int i=0; i<zs.length; i++) {
        for(int j=0; j<zs[i].size(); j++)
        stateCount.adjustOrPutValue(zs[i].get(j), 1, 1);
      }
      TreeMap<Double, Integer> sortedStates = new TreeMap<>(Collections.reverseOrder());
      TIntIntIterator it = stateCount.iterator();
      while(it.hasNext()) {
        it.advance();
        double freq = it.value();
        while(sortedStates.containsKey(freq)) { //break ties arbitrarily
          freq -= 0.0000001; //HACK (gross)
        }
        sortedStates.put(freq, it.key());
      }
      int ct = 0;
      for(double d : sortedStates.keySet()) {
        stateTranslation.put(sortedStates.get(d), ct++);
      }
      translateStates(zs, stateTranslation);
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
    int maxVal = _struct.numLeaves();
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
        if(GIBBS) {
          initZ(_c._docs, maxVal, true); //init scratch array to zs
          gibbsPass(offset);
          if((i % updateInterval)==0) { 
            optimizeParameters(_c._scratchZ);
//            if(j==0)
//              renumberStatesForCompression(_c._scratchZ);
          }
          updateModel(_c._scratchZ);
          _c._z = _c._scratchZ;
        }
        else {
          this.initArraysFromSBTs();
          emPass();
          if((i % updateInterval)==0) { 
            for(int k=0; k<this._forwardDelta.length;k++) {
              _forwardDelta[k] *= 0.995;
              _backwardDelta[k] *= 0.995;
              _wordDelta[k] *= 0.995;
            }
          }
          updateModelFromEM();
        }
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
    if(out._endArray != null) {
      System.out.println("arrays found, loading from arrays.");
      out.updateModelFromEM();
    }
    else
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
//      for(int j=0; j<dist[i].length;j++) {
//        sbtForward.addWeighted(this._forward[j], dist[i][j]);
//      }
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

  public static void extendTraining(String modelFile, String inputFile, String outputFile, String configFile) throws Exception {
    SBTSequenceModel sbtsm = readModel(modelFile);
    if(sbtsm._c._changeFactor == null) {
      sbtsm._c._changeFactor = new TDoubleArrayList[sbtsm._c._z.length];
      for(int i=0; i<sbtsm._c._changeFactor.length; i++) {
        sbtsm._c._changeFactor[i] = new TDoubleArrayList();
        for(int j=0; j<sbtsm._c._z[i].size();j++)
          sbtsm._c._changeFactor[i].add(1.0);
      }
      sbtsm._c.lambda = 0.99; //HACK
    }
    sbtsm.readConfigFile(configFile);
    sbtsm._expansionBranchFactors = sbtsm._branchingFactors;
    sbtsm.trainModel(sbtsm._NUMITERATIONS, sbtsm._OPTIMIZEINTERVAL, outputFile, true);
    sbtsm.writeModel(outputFile);
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
        ws[i] = Math.max(0,  ws[i]); //don't go below zero
        normalizer += ws[i];
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
   * @param modelFile Contains the model
   * @param outputFile  Where to write the output
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
  
  //assumes sorted
  public static double intersectTwo(TIntDoubleHashMap a, TIntDoubleHashMap b, TIntDoubleHashMap c) {
     double prod = 0.0;
     for(int k : a.keys()) {
       if(b.containsKey(k)) {
         double d =a.get(k)*b.get(k); 
         prod += d;
         if(c != null)
           c.put(k, d);
       }
     }
     return prod;
  }
  
  public static double intersectThree(TIntDoubleHashMap a, TIntDoubleHashMap b, TIntDoubleHashMap c) {
    TIntDoubleHashMap [] hms = new TIntDoubleHashMap [] {a, b, c};
    for(int i=0; i<hms.length;i++)
      if(hms[i]==null)
        hms[i] = new TIntDoubleHashMap();
    Arrays.sort(hms, (x, y) -> Integer.compare(x.size(), y.size()));
    TIntDoubleHashMap firstTwo = new TIntDoubleHashMap();
    double out = intersectTwo(hms[0], hms[1], firstTwo);
    out += intersectTwo(hms[0], hms[2], new TIntDoubleHashMap());
    out += intersectTwo(hms[1], hms[2], new TIntDoubleHashMap());
    out += intersectTwo(firstTwo, hms[2], new TIntDoubleHashMap());
    return out;
  }
  
  //returns 3-dim array with {time, count, totalFactor}
  public static long [] timePass(boolean useSBT, TIntDoubleHashMap [] hmF,
      TIntDoubleHashMap [] hmW, 
      TIntDoubleHashMap [] hmB,
      SBTSequenceModel sbtsm,
      int maxCount,
      int factorLow,
      int factorHigh) {
    long [] out = new long [3];
    double sum = 0.0;
    long start = System.currentTimeMillis();
    long totalFactor = 0L;
    for(int i=0; i<sbtsm._c._docs.length; i++) {
      TIntArrayList doc = sbtsm._c._docs[i];
      for(int j=1; j<doc.size()-1;j++) {
        int zPrev = sbtsm._c._z[i].get(j-1);
        int zNext = sbtsm._c._z[i].get(j+1);
        int wIdx = doc.get(j);
        TIntDoubleHashMap fw = hmF[zPrev];
        TIntDoubleHashMap bw = hmB[zNext];
        TIntDoubleHashMap word = hmW[wIdx];
        int [] sizes = new int [] {fw.size(), bw.size(), word.size()};
        Arrays.sort(sizes);
        int factor = 2 * sizes[0] + sizes[1];
        if(factor >= factorHigh || factor < factorLow)
          continue;
        if(useSBT) {
          SparseBackoffTree fwSBT = sbtsm._forward[zPrev];
          SparseBackoffTree bwSBT = sbtsm._backward[zNext];
          SparseBackoffTree wordSBT = sbtsm._wordToState[wIdx];
          sum += (new SparseBackoffTreeIntersection(new SparseBackoffTree [] {fwSBT, bwSBT, wordSBT}))._totalMass;
        }
        else
          sum += intersectThree(fw, bw, word);
        totalFactor += factor;
        if(++out[1] == maxCount)
          break;
      }
      if(out[1] == maxCount)
        break;
    }
    out[0] = System.currentTimeMillis() - start;
    out[2] = totalFactor;
    System.out.println(sum);
    return out;
  }
  
  public static void timeTrial(String modelFile) throws Exception {
    SBTSequenceModel sbtsm = readModel(modelFile);
    System.out.println("done reading.");
    TIntDoubleHashMap [] hmF = sbtsm.aggregateCounts(-1, sbtsm._c._z, sbtsm._c._docs);
    TIntDoubleHashMap [] hmW = sbtsm.aggregateCounts(0, sbtsm._c._z, sbtsm._c._docs);
    TIntDoubleHashMap [] hmB = sbtsm.aggregateCounts(1, sbtsm._c._z, sbtsm._c._docs);
    System.out.println("Done building hash maps.");
    System.out.println("Sampling 20000 words.");
    
    int [] buckets = new int [] {0, 62, 125, 250, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 6000, 8000}; 
    int [] cts = new int [buckets.length - 1];
    int [] times = new int[buckets.length - 1];
    int [] timesHash = new int[buckets.length - 1];
    int [] factors = new int[buckets.length - 1];
    
    for(int b=1; b<buckets.length;b++) {
      long [] timeCt = timePass(true, hmF, hmW, hmB, sbtsm, 20000, buckets[b-1], buckets[b]);
      times[b-1] = (int)timeCt[0];
      cts[b-1] = (int)timeCt[1];
      long [] timeCtHash = timePass(false, hmF, hmW, hmB, sbtsm, 20000, buckets[b-1], buckets[b]);
      timesHash[b-1] = (int)timeCtHash[0];
      factors[b-1] = (int)timeCtHash[2];
    }
    System.out.println("buckets: " + Arrays.toString(buckets));
    System.out.println("SBT times: " + Arrays.toString(times));
    System.out.println("Hash times: " + Arrays.toString(timesHash));
    System.out.println("Total factors: " + Arrays.toString(factors));
    System.out.println("cts: " + Arrays.toString(cts));
    for(int i=1; i<buckets.length;i++) {
      System.out.println(buckets[i] + "\t" + times[i-1] + "\t" + timesHash[i-1] + "\t" + cts[i-1] + "\t" + factors[i-1]);
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
    else if(args.length > 0 && args[0].equalsIgnoreCase("timeTrial")) {
      if(args.length != 2) {
        System.err.println("Usage: timeTrial <model_file>");
        return;
      }
      else {
        timeTrial(args[1]);
      }
    }
    else if(args.length > 0 && args[0].equalsIgnoreCase("trainFromModel")) { //expands existing model and continues
      if(args.length != 5) {
        System.err.println("Usage: trainFromModel <input_file> <model_output_file> <configuration_file> <model_input_file>");
        return;
      }
      extendTraining(args[4], args[1], args[2], args[3]);
    }
    else if(args.length > 0 && args[0].equalsIgnoreCase("trainAug")) { //augment existing probs
      if(args.length != 5) {
        System.err.println("Usage: trainAug <input_file> <model_output_file> <configuration_file> <probs_file>");
        return;
      }
      train(args[1], args[2], args[3], args[4]);
    }
    else if(args.length > 0 && args[0].equalsIgnoreCase("test")) {
      if(args.length != 6) {
        System.err.println("Usage: test <model_file> <test_file> <configuration_file> <num_docs_in_test_file> <num_docs_to_test>");
        return;
      }
      test(args[1], args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3]);
    }
    else if(args.length > 0 && args[0].equalsIgnoreCase("wordreps")) {
      if(args.length != 3) {
        System.err.println("Usage: wordreps <model_file> <output_file>");
        return;
      }
      outputWordToTopicFile(args[1], args[2]);
    }
    else if(args.length > 0 && args[0].equalsIgnoreCase("testEns")) {
      if(args.length != 6 && args.length != 8) {
        System.err.println("Usage: test <model_dir> <test_file> <configuration_file> <num_docs_in_test_file> <num_docs_to_test>"
            + " [<validation_file> <num_validation_docs>]");
        return;
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
        return;
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
