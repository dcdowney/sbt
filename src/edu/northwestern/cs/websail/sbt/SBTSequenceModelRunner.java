package edu.northwestern.cs.websail.sbt;
import mpi.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;


import edu.berkeley.nlp.lm.ArrayEncodedProbBackoffLm;
import edu.berkeley.nlp.lm.ContextEncodedProbBackoffLm;
import edu.berkeley.nlp.lm.array.LongArray;
import edu.berkeley.nlp.lm.bits.BitList;
import edu.berkeley.nlp.lm.bits.CompressionUtils;
import edu.berkeley.nlp.lm.collections.TIntMap;
import edu.berkeley.nlp.lm.io.LmReaders;
import edu.berkeley.nlp.lm.util.Logger;
import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.map.hash.TIntDoubleHashMap;

public class SBTSequenceModelRunner {

  LongArray _wData; //word to topic data
  LongArray _tData; //transition data; encoding is 0-delimited records of the forward transition counts,
                    //where each record is <offset to next ID> <count> <offset to next ID> <count> 
  LongArray _tStartData; //transition data from start state
  
  SparseBackoffTreeStructure _struct;
  
  TIntMap<String> _dict; //top 26 bits is pointer to starting long in wData, bottom six bits are offset w/in that long
  TIntMap<String> _dictLen; //length of dictionary entry.  Shares keys with _dict
  
  public static final int RADIX = 4;
  
  
  public SBTSequenceModelRunner(String sbtFile) throws Exception {
    SBTSequenceModel sbtsm = SBTSequenceModel.readModel(sbtFile);
    _struct = sbtsm._struct;
    _tData = encodeTransData(sbtsm);
    _tStartData = encodeStartData(sbtsm);
    _dict = new TIntMap<String>();
    _dictLen = new TIntMap<String>();
    
    _wData = encodeWordData(sbtsm);
    System.out.println("run berkeley compression to get dictionary size.");
    System.out.println("transition model size: " + _tData.size()*64 + " bits.");
    System.out.println("start model size: " + _tStartData.size()*64 + " bits.");
    System.out.println("word model size: " + _wData.size()*64 + " bits.");
    System.out.println("total: " + (_tData.size()*64 +_tStartData.size()*64 + _wData.size()*64));
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
    	return Integer.valueOf(i).compareTo(id.i);
    }
  }
  
  public BitList bitsForHashMap(TIntDoubleHashMap hm) {
    BitList bl = new BitList();
    TIntDoubleIterator it = hm.iterator();
    ArrayList<IntDouble> al = new ArrayList<IntDouble>();
    
    while(it.hasNext()) {
      it.advance();
      int k = it.key();
      int v = (int)Math.round(it.value());
      al.add(new IntDouble(k, v));
    }
    Collections.sort(al);
    int kLast = 0;
    for(int i=0; i<al.size(); i++) {
      int k = al.get(i).i;
      int v = (int)Math.round(al.get(i).d);
      bl.addAll(CompressionUtils.variableCompress(k-kLast, RADIX));
      bl.addAll(CompressionUtils.variableCompress(v, RADIX));
      kLast = k;
    }
    return bl;
  }
  
  public LongArray encodeStartData(SBTSequenceModel sbtsm) {
    TIntDoubleHashMap [] ss = sbtsm.aggregateCounts(-2, sbtsm._c._z, sbtsm._c._docs);
    LongArray out = new LongArray(100L);
    writeBlockToArray(bitsForHashMap(ss[0]), out);
    return out;
  }
  
  public LongArray encodeTransData(SBTSequenceModel sbtsm) {
    LongArray out = new LongArray(100L);
    TIntDoubleHashMap [] fw = sbtsm.aggregateCounts(-1, sbtsm._c._z, sbtsm._c._docs);
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
    int longIdx = 0;
    int offset = 0;
    BitList bOut = new BitList();
    int wordParams = 0;
    for(int i=0; i<wc.length; i++) {
      TIntDoubleHashMap hm = wc[i];
      if(hm == null)
        continue;
      BitList bl = new BitList();
      bl.addAll(bitsForHashMap(hm));
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
  
  public SparseBackoffTree sbtForWord(String w) {
    SparseBackoffTree out = new SparseBackoffTree(_struct);
    
    int loc = _dict.get(w, -1);
return out;    
    
//    int offset = loc % 128;
//    int longIdx = loc >>> 6;
//    if(i==_dicthigh.size() - 1) {
//      w2EndIdx = (int)_data.size() - 1;
//      w2EndLow = 63;    
//    }
//    else {
//      w2EndIdx = _dicthigh.values[i+1];
//      w2EndLow = _dictlow.values[i+1];
//    }
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
      SBTSequenceModel sbtsm = SBTSequenceModel.readModel(inFile);
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
  
	public static void main(String[] args) throws Exception {
	  //testWordEncoding();
	  getEnsembleParams("e:\\data\\10m\\ens_6\\");
	  //testBerkSize("x:\\data\\10m\\trigramModelAll");
	}

}
