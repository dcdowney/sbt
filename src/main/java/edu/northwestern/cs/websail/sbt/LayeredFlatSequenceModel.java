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

public class LayeredFlatSequenceModel implements Serializable {

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
  int _NUMLAYERS = 3; //number of layers of latent states in the sequence model
  int [] _NUMSTATES = new int [] {4, 13, 41}; //number of latent state values at each layer



  boolean [][] _transitionTemplate = new boolean[][] {{true, false,false},{false, true, false},{false,false,true}}; //dims: MUST BE _NUMLAYERS x _NUMLAYERS
  //transitionTemplate[i][j] says whether state at layer i is connected to state at layer j in the next time step
  boolean [] _productionTemplate = new boolean[] {false, false, true}; //dims: MUST BE _NUMLAYERS
//  boolean [][] _transitionTemplate = new boolean[][] {{true}}; //dims: MUST BE _NUMLAYERS x _NUMLAYERS
//  boolean [] _productionTemplate = new boolean[] {true}; //dims: MUST BE _NUMLAYERS
  
  
  
  //assumed internal links are always from layer i to layer i+1
  
  TIntDoubleHashMap [][][] _forward; //dims: LAYERi-1 x LAYERi x state
  TIntDoubleHashMap []  _startState; //dims: LAYER
  TIntDoubleHashMap [][][] _backward; //dims: LAYERi+1 x LAYERi x state
  TIntDoubleHashMap [][] _wordToState; //dims: LAYER x state
  TIntDoubleHashMap [][] _down; //dims: LAYER-1 x state
  TIntDoubleHashMap [][] _up; //dims: LAYER-1 x state
  
  double [][] _marginal; //dims: LAYER x state, gives total counts in each state
  
  
  double _alpha = 10.0; //pseudo-counts to add to each state  

  LayeredCorpus _c;
  
  Random _r = new Random();
    
  //threading:
  private int _NUMTHREADS = -1;
  private int[] _THREADBREAKS = null; //inclusive doc indices where threads *end* (initial thread starts
              //at 0 implicitly)
  
  //training config:
  protected int _NUMITERATIONS = -1;
  
  //testing config:
  protected int _NUMPARTICLES = -1;
  private int _MAXTESTDOCSIZE = -1;
  
  private int _numSkipWords = 0;
  
  public LayeredFlatSequenceModel(String configFile) throws Exception {
    _c = new LayeredCorpus(_NUMLAYERS);
    readConfigFile(configFile);
  }
  
  
  public int get_NUMLAYERS() {
    return _NUMLAYERS;
  }
  
  public void set_NUMLAYERS(int _NUMLAYERS) {
    this._NUMLAYERS = _NUMLAYERS;
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
    for(int i=0; i<doc.size(); i++) {
      for(int layer=0; layer<_NUMLAYERS; layer++) {
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
    TIntDoubleHashMap [] sbts;
  }
  
  //Returns the SBTs and amounts to subtract for the given latent variable
  public MarkovBlanketContainer getMarkovBlanket(int doc, int layer, int pos, boolean useWord) {
    int w = _c._docs[doc].get(pos);
    ArrayList<TIntDoubleHashMap> sbts = new ArrayList<>();
    //forward:
    if(pos==0) {
      sbts.add(this._startState[layer]);
    }
    else {
      for(int i=0; i<_NUMLAYERS; i++) {
        if(this._transitionTemplate[i][layer]) {
          if(_c._z[doc][i].get(pos-1) >= 0) {
            TIntDoubleHashMap sbt = this._forward[i][layer][_c._z[doc][i].get(pos-1)];
            sbts.add(sbt);
          }
        }
      }
    }
    //backward:
    if(pos!=_c._docs[doc].size() - 1) {
      for(int i=0; i<_NUMLAYERS; i++) {
        if(this._transitionTemplate[layer][i]) {
          if(_c._z[doc][i].get(pos+1) >= 0) {
            TIntDoubleHashMap sbt = this._backward[i][layer][_c._z[doc][i].get(pos+1)];
            sbts.add(sbt);
          }
        }
      }
    }
    if(useWord) {
      //word:
      if(_productionTemplate[layer]) {
        TIntDoubleHashMap sbtWord =this._wordToState[layer][w]; 
          sbts.add(this._wordToState[layer][w]);
      }
    }
      //layers:
    if(layer!=_NUMLAYERS-1) {
      if(_c._z[doc][layer+1].get(pos) >= 0) {
        sbts.add(this._up[layer+1][_c._z[doc][layer+1].get(pos)]);
      }
    }
    if(layer!=0) {
      if(_c._z[doc][layer-1].get(pos) >= 0) {
        sbts.add(this._down[layer-1][_c._z[doc][layer-1].get(pos)]);
      }
    }
    
    MarkovBlanketContainer out = new MarkovBlanketContainer();
    out.sbts = new TIntDoubleHashMap[sbts.size()];
    for(int i=0; i<sbts.size(); i++) {
      out.sbts[i] = sbts.get(i);
    }
    return out;
  }
  
  
  public int sampleZ(int doc, int layer, int pos, boolean useWord, boolean sub) {
    MarkovBlanketContainer mbc = getMarkovBlanket(doc, layer, pos, useWord);
 
    int curZ = _c._z[doc][layer].get(pos);
    double total = 0.0;
    double [] distro = new double[_NUMSTATES[layer]];
    for(int i=0; i<_NUMSTATES[layer]; i++) {
      double adj = 0.0;
      if(sub && i==curZ)
        adj = -1.0;
      distro[i] = 1.0;
      boolean first = true;
      for(int j=0; j<mbc.sbts.length; j++) {
        if(mbc.sbts[j]!=null) { //ignore nulls
          if(!first)
            distro[i] *= (mbc.sbts[j].get(i) + _alpha/_NUMSTATES[layer] + adj)/(_marginal[layer][i] + _alpha  + adj);
          else {
            first = false;
            distro[i] *= (mbc.sbts[j].get(i) + _alpha/_NUMSTATES[layer] + adj); //no normalization on first one
          }
        }
      }
      total += distro[i];
    }
    double ch = total*_r.nextDouble();
    int sample = -1;
    for(int i=0; i<_NUMSTATES[layer];i++) {
      if(ch<distro[i]) {
        sample = i;
        break;
      }
      else
        ch -= distro[i];
    }
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
        ac.down[i] = new TIntDoubleHashMap[_NUMSTATES[i]];
      }
      if(i > 0) {
        ac.up[i] = new TIntDoubleHashMap[_NUMSTATES[i]];
      }
      for(int j=0; j<_NUMLAYERS; j++) {
        ac.forward[i][j] = new TIntDoubleHashMap[_NUMSTATES[i]];
        ac.backward[i][j] = new TIntDoubleHashMap[_NUMSTATES[i]];
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
  public static double [] getNormalizers(TIntDoubleHashMap [] shds, int numStates) {
    double [] out = new double[numStates];
    for(int i=0; i<numStates; i++) {
      for(int j=0; j<shds.length; j++)
        out[i] += shds[j].get(i);
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
      for(int i=0; i<_NUMLAYERS;i++) {
        System.out.println("second doc samples (layer " + i + "): " + zs[1][i].toString());
      }
//      if(_c._changeFactor != null)
//        System.out.println("second doc chg (layer 0): " + _c._changeFactor[1][0].toString());
      AggregatedCounts ac = this.aggregateCounts(zs, _c._docs);
      
      this._startState  = new TIntDoubleHashMap[_NUMLAYERS];
      this._wordToState  = new TIntDoubleHashMap[_NUMLAYERS][];
      this._forward  = new TIntDoubleHashMap[_NUMLAYERS][_NUMLAYERS][];
      this._backward = new TIntDoubleHashMap[_NUMLAYERS][_NUMLAYERS][];
      this._up  = new TIntDoubleHashMap[_NUMLAYERS][];
      this._down  = new TIntDoubleHashMap[_NUMLAYERS][];
      this._marginal = new double[_NUMLAYERS][];
      
      for(int i=0; i<_NUMLAYERS;i++) {
        _marginal[i] = new double[_NUMSTATES[i]];
        this._startState[i] = ac.start[i];
        this._wordToState[i] = ac.word[i];
        for(int j=0; j<_NUMLAYERS;j++) {
          this._forward[i][j] = ac.forward[i][j];
          this._backward[j][i] = ac.backward[j][i];
        }
        if(i > 0)
          this._down[i-1] = ac.down[i-1];
        if(i < _NUMLAYERS - 1)
          this._up[i+1] = ac.up[i+1];
      }
      //re-compute marginal:
      for(int i=0; i<zs.length;i++) {
        for(int j=0; j<zs[i].length; j++) {
          for(int k=0; k<zs[i][j].size(); k++) {
            _marginal[j][zs[i][j].get(k)]++;
          }
        }
      }
      for(int i=0; i<_marginal.length;i++)
        System.out.println("marginal[" + i + "]: " + Arrays.toString(_marginal[i]));
    }
    
    
    
    
    
    public void trainModel(int iterations, String outFile) throws Exception {
      trainModel(iterations, outFile, false);
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
  public void trainModel(int iterations, String outFile, boolean expandOnEntry) throws Exception {
    boolean done = false;
    int j=0;
    boolean PARTIALPASS = true;
    while(!done) {
      for(int i=1; i<=iterations; i++) {
        System.out.print(i + ": ");
        int offset = 0;
        if(PARTIALPASS) {
          offset = setThreadBreaks(5, (i-1) % 5);
          System.out.println(offset + "\t" + Arrays.toString(_THREADBREAKS));
        }
        initZ(_c._docs, _NUMSTATES, true); //init scratch array to zs
        gibbsPass(offset);
        updateModel(_c._scratchZ);
        _c._z = _c._scratchZ;
      }
      done = true;
    }
  }
  
  public void writeModel(String outFile) throws Exception {
    TIntDoubleHashMap [][] sbtW = this._wordToState;
    TIntDoubleHashMap [][][] sbtF = this._forward;
    TIntDoubleHashMap [][][] sbtB = this._backward;
    TIntDoubleHashMap [] sbtS = this._startState;
    TIntDoubleHashMap [][] sbtU = this._up;
    TIntDoubleHashMap [][] sbtD = this._down;
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
    _startState = sbtS;
    _up = sbtU;
    _down = sbtD;
  }
  
  public static LayeredFlatSequenceModel readModel(String inFile) throws Exception {
    ObjectInputStream ois = new ObjectInputStream(new FileInputStream(inFile));
    LayeredFlatSequenceModel out = (LayeredFlatSequenceModel) ois.readObject();
    ois.close();
    out.updateModel(out._c._z);
    return out;
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
        for(int pass=0; pass<3; pass++)
          for(int layer=0;layer<_NUMLAYERS;layer++) {
            _c._z[doc][layer].set(i, this.sampleZ(doc, layer, i, false, false));
          }
        
        if(_c._pWord[w]>0.0) {
          
          //test on this word:
          double [] pOut = new double[this.get_VOCABSIZE()];
          double normalizer = 0.0;
          for(int k=0; k<pOut.length; k++) {
            if(_c._pWord[k] >0.0) {
              int times = 0;
              pOut[k] = 1.0;
              for(int layer=0; layer<_NUMLAYERS; layer++) {
                if(this._productionTemplate[layer]) {
                  int zNow = _c._z[doc][layer].get(i);
                  pOut[k] *= (this._wordToState[layer][k].get(zNow)+_alpha/pOut.length)/(this._marginal[layer][zNow] + _alpha);
                  times++;
                }
              }
              //divide out redundant factors of p(word):
              for(;times>1;times--)
                pOut[k] /= _c._pWord[k];
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
//        //backward:
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
//        //current:
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
        if(ps[i]==0.0)
          System.out.println("problems.");
        LL += Math.log(ps[i]/(double)_NUMPARTICLES);
        if(Double.isNaN(LL))
          System.out.println("got NaN with " + ps[i]);
      }
    }
    return new double [] {LL, numWords};
  }
  
  public static void train(String inputFile, String outputFile, String configFile) throws Exception {
    LayeredFlatSequenceModel sbtsm = new LayeredFlatSequenceModel(configFile);
    sbtsm.initializeForCorpus(inputFile, sbtsm._NUMSTATES);
    sbtsm.trainModel(sbtsm._NUMITERATIONS, outputFile);
    sbtsm.writeModel(outputFile);
  }
  
    
  
  //computes log likelihood of model using left-to-right method
  //returns {ppl, number of tested words}
  //numDocs is number in test file
  //maxDocs is number actually tested (starting from the beginning of the file)
  public static double [] testModel(String modelFile, String testFile, int numDocs, int maxDocs, String configFile) throws Exception {
    
    LayeredFlatSequenceModel sbtsm = readModel(modelFile);
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
