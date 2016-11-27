package edu.northwestern.cs.websail.sbt;

import java.io.File;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import edu.berkeley.nlp.lm.ArrayEncodedProbBackoffLm;
import edu.berkeley.nlp.lm.ContextEncodedProbBackoffLm;
import edu.berkeley.nlp.lm.array.LongArray;
import edu.berkeley.nlp.lm.bits.BitList;
import edu.berkeley.nlp.lm.bits.CompressionUtils;
import edu.berkeley.nlp.lm.collections.TIntMap;
import edu.berkeley.nlp.lm.io.LmReaders;
import edu.berkeley.nlp.lm.util.Logger;
import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.map.hash.TIntIntHashMap;

public class SBTSequenceModelRunner {

  LongArray _wData; //word to topic data
  LongArray _tData; //transition data; encoding is 0-delimited records of the forward transition counts,
                    //where each record is <offset to next ID> <count> <offset to next ID> <count> 
  LongArray _tStartData; //transition data from start state
  
  SparseBackoffTreeStructure _struct;
  
  SBTSequenceModel _sbtsm; //temporary
  
  TIntMap<String> _dict; //top 26 bits is pointer to starting long in wData, bottom six bits are offset w/in that long
  TIntMap<String> _dictLen; //length of dictionary entry.  Shares keys with _dict
  
  TIntDoubleHashMap [] _wordToState; //temporary word-to-state until decompression implementation complete
  TIntDoubleHashMap [] _stateToState; //temporary state-to-state
  TIntDoubleHashMap _startState; //temporary start state
  
  TIntIntHashMap _translationTable; //maps stored state onto real state

  public int _singles = 0;
  
  double [] _wordMarginal; //holds marginal to divide by for words
      //TODO: remove this dependency?
  double [] _wordDs;
  double [] _forwardDs;
  
  double [] _wordDivider; //holds marginal to divide *counts* by for words
  
  public static final int RADIX = 3;
  
  
  public SBTSequenceModelRunner(String sbtFile) throws Exception {
    SBTSequenceModel sbtsm = SBTSequenceModel.readModel(sbtFile);
    _sbtsm = sbtsm;
    _struct = sbtsm._struct;
    _wordMarginal = SBTSequenceModel.getCheckMarginal(sbtsm._wordToState, sbtsm._struct);
    _translationTable = getTranslationTable(sbtsm);

    _tData = encodeTransData(sbtsm);
    _tStartData = encodeStartData(sbtsm);
    _dict = new TIntMap<String>();
    _dictLen = new TIntMap<String>();
    _forwardDs = sbtsm._forwardDelta;
    _wordDs = sbtsm._wordDelta;
    _wordDivider = sbtsm._wordStateMarginal;
    _wData = encodeWordData(sbtsm);
    System.out.println("run berkeley compression to get dictionary size.");
    System.out.println("transition model size: " + _tData.size()*64 + " bits.");
    System.out.println("start model size: " + _tStartData.size()*64 + " bits.");
    System.out.println("word model size: " + _wData.size()*64 + " bits.");
    System.out.println("total: " + (_tData.size()*64 +_tStartData.size()*64 + _wData.size()*64));
    System.out.println("NEED TO ADD WORD MARGINAL and discounts");
    System.out.println("singles: " + _singles);
  }
  
  private class TestDoer implements Runnable {
    TIntArrayList _doc;
    int _docID;
    double [] _res;
    
    public TestDoer(TIntArrayList doc, int docID) {
      _doc = doc;
      _docID = docID;
    }
    
    public void run() {
      _res = testModel(_doc, _wordDs, _forwardDs);
      //_res = testOnDocFullPplExact(_doc, _wordTopicMarginal);
      System.out.println(_docID + "\t" + _res[0] + "\t" + _res[1]);
    }
  }
  
  private static class TestEnsembleDoer implements Runnable {
    TIntArrayList _doc;
    int _docID;
    SBTSequenceModelRunner [] _runner;
    double [] _weights;
    double [] _res;
    
    public TestEnsembleDoer(TIntArrayList doc, int docID, SBTSequenceModelRunner [] runner, double [] weights) {
      _runner = runner;
      _doc = doc;
      _docID = docID;
      _weights = weights;
    }
    
    public void run() {
      _res = testEnsemble(_doc, _runner, _weights);
      System.out.println(_docID + "\t" + _res[0] + "\t" + _res[1]);
    }
  }
  
  /**
   * @param blockBits
   * @param array
   */
  public static void writeBlockToArray(final BitList blockBits, final LongArray array) {
    long curr = 0L;
    for (int i = 0; i <= Long.SIZE * (1 + blockBits.size() / Long.SIZE); ++i) {
      if (i % Long.SIZE == 0 && i > 0) {
        array.add(curr);
        curr = 0;
      }
      curr = (curr << 1) | ((i >= blockBits.size() || !blockBits.get(i)) ? 0 : 1);
    }
  }
  
  private static class IntDouble implements Comparable<IntDouble> {
    public int i;
    public double d;
    
    public IntDouble(int _i, double _d) {
      i = _i;
      d = _d;
    }
    
    public int compareTo(IntDouble id) {
      return Integer.compare(i,  id.i);
    }
  }
  
  public BitList bitsForHashMap(TIntDoubleHashMap hm) {
    BitList bl = new BitList();
    TIntDoubleIterator it = hm.iterator();
    ArrayList<IntDouble> al = new ArrayList<IntDouble>();
    
    while(it.hasNext()) {
      it.advance();
      int k = it.key();
//      k = _translationTable.get(k);
      int v = (int)Math.round(it.value());
      al.add(new IntDouble(k, v));
    }
    Collections.sort(al);
    int kLast = 0;
    for(int i=0; i<al.size(); i++) {
      int k = al.get(i).i;
      int v = (int)Math.round(al.get(i).d);
      if(v==1)
        _singles++;
      bl.addAll(CompressionUtils.variableCompress(k-kLast, RADIX));
      bl.addAll(CompressionUtils.variableCompress(v-1, RADIX));
      kLast = k;
    }
    return bl;
  }
  
  public LongArray encodeStartData(SBTSequenceModel sbtsm) {
    TIntDoubleHashMap [] ss = sbtsm.aggregateCounts(-2, sbtsm._c._z, sbtsm._c._docs);
    LongArray out = new LongArray(100L);
    writeBlockToArray(bitsForHashMap(ss[0]), out);
    _startState = ss[0];
    return out;
  }
  
  public TIntIntHashMap getTranslationTable(SBTSequenceModel sbtsm) {
    TIntDoubleHashMap [] fw = sbtsm.aggregateCounts(-1, sbtsm._c._z, sbtsm._c._docs); 
    TIntIntHashMap freq = new TIntIntHashMap();

    for(int i=0; i<fw.length; i++) {
      if(fw[i] == null)
        continue;
      TIntDoubleIterator it = fw[i].iterator();
      while(it.hasNext()) {
        it.advance();
        freq.adjustOrPutValue(it.key(), 1, 1);
      }
    }

    TIntDoubleHashMap [] wordSt = sbtsm.aggregateCounts(0, sbtsm._c._z, sbtsm._c._docs); 

    for(int i=0; i<wordSt.length; i++) {
      if(wordSt[i] == null)
        continue;
      TIntDoubleIterator it = wordSt[i].iterator();
      while(it.hasNext()) {
        it.advance();
        freq.adjustOrPutValue(it.key(), 1, 1);
      }
    }
    
    ArrayList<Map.Entry<Integer, Integer>> stateFreq = new ArrayList<>();
    TIntIntIterator it = freq.iterator();
    while(it.hasNext()) {
      it.advance();
      Map.Entry<Integer, Integer> e = new AbstractMap.SimpleImmutableEntry<Integer, Integer>(it.key(), it.value());
      stateFreq.add(e);
    }
    //sort by value desc:
    stateFreq.sort(new Comparator<Map.Entry<Integer, Integer>>(){
        @Override
        public int compare(Map.Entry<Integer, Integer> e1, Map.Entry<Integer, Integer> e2) {
          return -e1.getValue().compareTo(e2.getValue());
        }
    });
    
    freq = new TIntIntHashMap(); //re-use
    int ct = 0;
    for(Map.Entry<Integer, Integer> e : stateFreq) {
      freq.put(e.getKey(), ct++);
    }
    return freq;
  }
  
  public LongArray encodeTransData(SBTSequenceModel sbtsm) {
    LongArray out = new LongArray(100L);
    TIntDoubleHashMap [] fw = sbtsm.aggregateCounts(-1, sbtsm._c._z, sbtsm._c._docs); 
    TIntIntHashMap freq = new TIntIntHashMap();
        
    //temporary:
    _stateToState = fw;
    
    BitList bl = new BitList();
    int totalParams = 0;
    int zeroStates = 0;
    for(int i=0; i<fw.length; i++) {
      TIntDoubleHashMap hm = fw[i];
      if(hm!=null) {
        bl.addAll(bitsForHashMap(hm));
        totalParams += hm.size();
      }
      else {
        zeroStates++;
      }
      bl.addAll(CompressionUtils.variableCompress(0, RADIX));
    }
    writeBlockToArray(bl, out);
    System.out.println("total trans params: " + totalParams);
    System.out.println("total zero states: " + zeroStates);
    return out;
  }
  
  public LongArray encodeWordData(SBTSequenceModel sbtsm) {
    LongArray out = new LongArray(100L);
    TIntDoubleHashMap [] wc = sbtsm.aggregateCounts(0, sbtsm._c._z, sbtsm._c._docs);
    this._wordToState = wc; //temporary
    int longIdx = 0;
    int offset = 0;
    BitList bOut = new BitList();
    int wordParams = 0;
    boolean OUTPUTWORDS = false;
    for(int i=0; i<wc.length; i++) {
      TIntDoubleHashMap hm = wc[i];
      if(hm == null)
        continue;
      BitList bl = new BitList();
      bl.addAll(bitsForHashMap(hm));
      if(OUTPUTWORDS)
        System.out.println(i + "\t" + bl.size());
      bOut.addAll(bl);
      int loc = longIdx << 6;
      loc += offset;
      this._dict.put(Integer.toString(i), loc);
      longIdx += bl.size() / 64;
      offset += bl.size() % 64;
      if(offset >= 64)
        longIdx++;
      offset %= 64;
      wordParams += hm.size();
    }
    writeBlockToArray(bOut, out);
    System.out.println("total word params: " + wordParams);
    return out;
  }
  
  public static String filePathOfResource(String path) {
    return SBTSequenceModelRunner.class.getResource(path).getFile();
  }
  
  public SparseBackoffTree sbtForWord(int w) {
    SparseBackoffTree out = new SparseBackoffTree(_struct);
    out.addAllMass(this._wordToState[w], _wordDs);
    return out;
    
//    int loc = _dict.get(w, -1);
//    int offset = loc % 128;
//    int longIdx = loc >>> 6;
//    
//    if(i==_dicthigh.size() - 1) {
//      w2EndIdx = (int)_data.size() - 1;
//      w2EndLow = 63;    
//    }
//    else {
//      w2EndIdx = _dicthigh.values[i+1];
//      w2EndLow = _dictlow.values[i+1];
//    }
//    return out;    
    

  }
  
  public static void testBerkSize(String modelFile) throws Exception {
    Logger.setGlobalLogger(new Logger.SystemLogger(System.out, System.err));
    ArrayEncodedProbBackoffLm<String> lm = LmReaders.readArrayEncodedLmFromArpa(modelFile, true);
     double d = lm.getLogProb(Arrays.asList("he", "said", "that"));
     System.out.println("prob is:" + d);
  }
  
  public static void getEnsembleParams(String inDirectory) throws Exception {
    File dir = new File(inDirectory);
    File [] files = dir.listFiles();
    int total = 0;
    for(int i=0; i<files.length; i++) {
      String inFile = files[i].getAbsolutePath();
      //SBTSequenceModel sbtsm = SBTSequenceModel.readModel(inFile);
      SBTSequenceModelRunner sbtsmr = new SBTSequenceModelRunner(inFile);
      total += 64*(sbtsmr._tData.size()+sbtsmr._tStartData.size() + sbtsmr._wData.size());
    }
    System.out.println("total size: " + total);
  }
  
  public static void testWordEncoding() throws Exception {
    //String inFile = filePathOfResource("/ptb.model.1");
    String inFile = "e:\\data\\ptb\\ens_2\\8down.2";
    SBTSequenceModel sbtsm = SBTSequenceModel.readModel(inFile);
    SBTSequenceModelRunner sbtsmr = new SBTSequenceModelRunner(inFile);
    System.out.println("w2s length: " + sbtsm._wordToState.length + " should be 1 + dict size: " + sbtsmr._dict.size());
    
  }
  
  public static void incrementNextState(double weight, double [] nextPState, double [] nextPStateCt, TIntDoubleHashMap hm) {
    if(hm==null)
      return;
    TIntDoubleIterator it = hm.iterator();
    double sum = 0.0;
    while(it.hasNext()) {
      it.advance();
      sum += it.value();
    }
    it = hm.iterator();
    final double C = weight/sum;
    while(it.hasNext()) {
      it.advance();
      nextPState[it.key()] += it.value()*C;
      nextPStateCt[it.key()] += C;
    }
  }
  
  public static double [] weightedSmooth(double ct, double [] ds) {
    double [] out = new double[ds.length];
    for(int i=0; i<out.length; i++)
      out[i] = ds[i]*ct;
    return out;
  }
  
  public SparseBackoffTree decompressWordMass(SparseBackoffTreeStructure str, int i) {
    
    SparseBackoffTree out = new SparseBackoffTree(str);
    out.addAllMass(_wordToState[i], _wordDs);
    out.divideCountsBy(this._wordDivider);
    return out;
  }
  
  public static void normalize(double [] d) {
    double sum = 0.0;
    for(int i=0; i<d.length;i++)
      sum += d[i];
    if(sum > 0.0)
      for(int i=0; i<d.length;i++)
        d[i] /= sum;
  }
  
  public static double [] testEnsemble(String modelDir, String testFile, int maxDocs) throws Exception {
    
    File [] models = new File(modelDir).listFiles();
    SBTSequenceModelRunner [] runner = new SBTSequenceModelRunner[models.length];
    double [] weights = new double[runner.length];
    for(int i=0; i<models.length;i++) {
      runner[i] = new SBTSequenceModelRunner(models[i].getAbsolutePath());
      weights[i] = 1.0/weights.length;
    }
    Corpus c = new Corpus();
    c.readCorpusDat(testFile, false);
    double out = 0.0;
    double toks = 0.0;
    
    TestEnsembleDoer [] tds = new TestEnsembleDoer[Math.min(c._docs.length, maxDocs)];
    ExecutorService e = Executors.newFixedThreadPool(runner[0]._sbtsm.get_NUMTHREADS());
    //ExecutorService e = Executors.newFixedThreadPool(1);
    for(int i=0; i<tds.length;i++) {
      tds[i] = new TestEnsembleDoer(c._docs[i], i, runner, weights);
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
      out += tds[i]._res[0];
      toks += tds[i]._res[1];
    }
    System.out.println("final LL: " + out);
    System.out.println("Total words: " + toks);
    System.out.println("ppl: " + Math.exp(-out/toks));
    return new double [] {out, toks};
  }
  
  public static double [] testSingleModel(String modelFile, String testFile, int docsToTest) throws Exception {
    SBTSequenceModelRunner runner = new SBTSequenceModelRunner(modelFile);
    Corpus c = new Corpus();
    c.readCorpusDat(testFile, false);
    double out = 0.0;
    double toks = 0.0;
    
    TestDoer [] tds = new TestDoer[Math.min(docsToTest, c._docs.length)];
    ExecutorService e = Executors.newFixedThreadPool(runner._sbtsm.get_NUMTHREADS());
    //ExecutorService e = Executors.newFixedThreadPool(1);
    for(int i=0; i<tds.length;i++) {
      tds[i] = runner.new TestDoer(c._docs[i], i);
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
      out += tds[i]._res[0];
      toks += tds[i]._res[1];
    }
    System.out.println("final LL: " + out);
    System.out.println("Total words: " + toks);
    System.out.println("ppl: " + Math.exp(-out/toks));
    return new double [] {out, toks};
  }
  
  //TODO: remove redundancy with trans
  public static void fillStartState(double [] nextPState, double [] nextPStateCt, TIntDoubleHashMap hm) {
    TIntDoubleIterator it = hm.iterator(); 
    while(it.hasNext()) {
      it.advance();
      nextPState[it.key()] += it.value();
      nextPStateCt[it.key()] += 1.0;
    }
  }
  
  public static double [] testEnsemble(TIntArrayList ws, SBTSequenceModelRunner [] runner,
      double [] weights) {
    int M = runner.length;
    int [] numStates = new int [M];
    
    double [][] pState = new double[M][];
    double out = 0.0;
    double [][] nextPState = new double[M][];
    double [][] nextPStateCt = new double[M][];
    for(int i=0; i<M; i++) {
      numStates[i] = runner[i]._struct.numLeaves();
      pState[i] = new double[numStates[i]];
      nextPState[i] = new double[numStates[i]];
      nextPStateCt[i] = new double[numStates[i]];
    }
    for(int w=0; w<ws.size(); w++) {
      double likelihood = 0.0;
      for(int mod=0;mod<M;mod++) {
        int word = ws.get(w);
        //compute next state dist:
        if(w==0) {
          fillStartState(nextPState[mod], nextPStateCt[mod], runner[mod]._startState);
        }
        else
          for(int i=0; i<numStates[mod];i++) {
            if(pState[mod][i] > 0.0) {
              incrementNextState(pState[mod][i], nextPState[mod], nextPStateCt[mod], runner[mod]._stateToState[i]);
            }
          }
        SparseBackoffTree sbtNext = new SparseBackoffTree(runner[mod]._struct);
        for(int i=0; i<numStates[mod];i++) {
          sbtNext.smoothAndAddMass(i, nextPState[mod][i], weightedSmooth(nextPStateCt[mod][i], runner[mod]._forwardDs));
        }
        SparseBackoffTree sbtWord = runner[mod].decompressWordMass(runner[mod]._struct, word);
  //      SparseBackoffTreeIntersection sbti = new SparseBackoffTreeIntersection(new SparseBackoffTree [] {sbtNext, sbtWord});
        pState[mod] = sbtNext.toDoubleArray();
        normalize(pState[mod]);
        double [] pWord = sbtWord.toDoubleArray();
        for(int i=0; i<numStates[mod]; i++) {
          pState[mod][i] *= pWord[i]/runner[mod]._wordMarginal[i];
          likelihood += pState[mod][i]*weights[mod];
        }
        normalize(pState[mod]);
        Arrays.fill(nextPState[mod], 0.0);
        Arrays.fill(nextPStateCt[mod], 0.0);
      }
      out += Math.log(likelihood);
    }
    return new double [] {out, ws.size()};
  }
  
  public double [] testModel(TIntArrayList ws, double [] wordDs, double[] stateDs) {
    int numStates = _struct.numLeaves();
    double [] pState = new double[numStates];
    double out = 0.0;
    double [] nextPState = new double[numStates];
    double [] nextPStateCt = new double[numStates];
    for(int w=0; w<ws.size(); w++) {
      int word = ws.get(w);
      //compute next state dist:
      if(w==0) {
        fillStartState(nextPState, nextPStateCt, _startState);
      }
      else
        for(int i=0; i<numStates;i++) {
          if(pState[i] > 0.0) {
            incrementNextState(pState[i], nextPState, nextPStateCt, _stateToState[i]);
          }
        }
      SparseBackoffTree sbtNext = new SparseBackoffTree(_struct);
      for(int i=0; i<numStates;i++) {
        sbtNext.smoothAndAddMass(i, nextPState[i], weightedSmooth(nextPStateCt[i], stateDs));
      }
      SparseBackoffTree sbtWord = decompressWordMass(_struct, word);
//      SparseBackoffTreeIntersection sbti = new SparseBackoffTreeIntersection(new SparseBackoffTree [] {sbtNext, sbtWord});
      pState = sbtNext.toDoubleArray();
      normalize(pState);
      double [] pWord = sbtWord.toDoubleArray();
      double likelihood = 0.0;
      for(int i=0; i<numStates; i++) {
        pState[i] *= pWord[i]/_wordMarginal[i];
        likelihood += pState[i];
      }
      out += Math.log(likelihood);
      normalize(pState);
      Arrays.fill(nextPState, 0.0);
      Arrays.fill(nextPStateCt, 0.0);
    }
    return new double [] {out, ws.size()};
  }
  
	public static void main(String[] args) throws Exception {
	  //testWordEncoding();
//	  testEnsemble("e:\\data\\10m\\ens_adv", "e:\\data\\10m\\test.dat", 3080);
//	  testSingleModel("e:\\data\\ptb\\gibbslarge\\ptb.model", "e:\\data\\ptb\\penn_test.dat", 3761);
//	    testSingleModel("e:\\data\\10m\\gibbslarge6\\10m.model.6", "e:\\data\\10m\\test.dat", 3080);
	   testEnsemble("e:\\data\\10m\\gibbens", "e:\\data\\10m\\test.dat", 3080);
	   //2000 full passes:
	   //18: 538
	   //18 5: 350
	   //18 5 3: 270
	   //18 5 3 3: 220
	   
	   //2000 1/4 passes:
	  //18: 532
	   //18 5: 352
	   //18 5 3: 273
	   //18 5 3 3: 223
	   //18 5 3 3 3: 187
	   
	   //3000 1/10 passes:
	   //18: 551
	   //18 5: 357
	   //18 5 3: 275
	   //18 5 3 3: 220
	   //18 5 3 3 3: 185
     //18 5 3 3 3 3: 167
     //18 5 3 3 3 3 3: 159

     //3000 1/10 passes (3):
     //18: 540
     //18 5: 348
     //18 5 3: 272
     //18 5 3 3: 219 
     //18 5 3 3 3: 185
     //18 5 3 3 3 3: 167
     //18 5 3 3 3 3 3: 158

     //3000 1/10 passes (2):
     //18: 572
     //18 5: 356
     //18 5 3: 270
     //18 5 3 3: 220 
     //18 5 3 3 3: 186
     //18 5 3 3 3 3: 168
     //18 5 3 3 3 3 3: 159
     
	    //3000 1/10 passes (4):
	   //18: 
	   //18 5: 
	   //18 5 3: 
	   //18 5 3 3: 222
	   //18 5 3 3 3: 186
	   //18 5 3 3 3 3: 167
	   //18 5 3 3 3 3 3: 158
  	 //18 5 3 3 3 3 3 3: 157
 	   //18 5 3 3 3 3 3 3: 160
	    
	   //3000 1/10 passes:
	    //17:
	    //17 13: 290 
	    //17 13 11: 190
	    //17 13 11 7: 166
	    //17 13 11 7 5: 162
	    //17 13 11 7 5 5: 168
      //17 13 11 7 5 5 3: 181
	    //17 13 11 7 5 5 3 3: 194
      
	    //3000 /110 passes"
	    //17 11 7 5: 174
	    //17 11 7 5 3: 162
	    //17 11 7 5 3 2: 159
	    //17 11 7 5 3 2 2: 158
	    
	  //no translation table:
//	  testSingleModel("e:\\data\\10m\\ens_adv\\8down.7", "e:\\data\\10m\\test.dat", 3080);
//5028532 train after 526 for 15 10 with 0.99 every 5
	  
//	  getEnsembleParams("e:\\data\\ptb\\ens_cheery\\");
	  //testBerkSize("x:\\data\\10m\\trigramModelAll");
	}
}
