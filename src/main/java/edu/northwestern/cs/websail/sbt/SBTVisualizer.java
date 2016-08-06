package edu.northwestern.cs.websail.sbt;

import java.io.*;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;

import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.map.hash.TIntDoubleHashMap;


//By Junhan Liu (@HERMANNITY), Zhaocheng Yu, Doug Downey

public class SBTVisualizer {

  //adapted from http://stackoverflow.com/questions/1265282/recommended-method-for-escaping-html-in-java
  public static String escapeHTML(String s) {
    StringBuilder out = new StringBuilder(Math.max(16, s.length()));
    for (int i = 0; i < s.length(); i++) {
        char c = s.charAt(i);
        //if (c > 127 || c == '"' || c == '<' || c == '>' || c == '&') {
        if (c == '"' || c == '\\') {
            out.append("\\" + c);
            
//            out.append((int) c);
//            out.append(';');
        } else {
            out.append(c);
        }
    }
    if(out.length() == 0)
      return "NULL";
    if(out.toString().equals("'"))
      return "SINGQ";
    return out.toString();
}
  
  public static HashMap<Integer, String> readDict(String inFile) throws Exception {
    BufferedReader brIn = new BufferedReader(
        new InputStreamReader(new FileInputStream(inFile), "UTF8" ));
    String sLine;
    HashMap<Integer, String> out = new HashMap<Integer, String>();
    int ct = 0;
    while((sLine = brIn.readLine())!=null) {
      String [] fields = sLine.split("\t");
      int id = Integer.parseInt(fields[1]);
      //int freq = Integer.parseInt(fields[2]);
      String w = fields[0];
      w = escapeHTML(w);
      out.put(id, w);
      if(++ct % 10000 == 0) {
        System.out.println(ct + " words read from dict");
      }
    }
    brIn.close();
    return out;
  }
  
  public static void outputWordToTopicFile(String modelFile, String outputFile, String modelType, String dictFile) throws Exception {
    
    SparseBackoffTreeStructure root;
    SparseBackoffTree [] wordSBTs;
    Corpus c;
    
    HashMap<Integer, String> dict = readDict(dictFile);
    
    //wouldn't be necessary if the doc model and sequence model inherited from a base class...
    if(modelType.equalsIgnoreCase("doc")) {
      SBTDocumentModel sbtdm = SBTDocumentModel.readModel(modelFile);
      root = sbtdm._struct;
      wordSBTs = sbtdm._topicGivenWord;
      c = sbtdm._c;
    }
    else {
      SBTSequenceModel sbtsm = SBTSequenceModel.readModel(modelFile);
      root = sbtsm._struct;
      wordSBTs = sbtsm._wordToState;
      c = sbtsm._c;
    }
    
    HashMap<Integer,TreeMap<Integer,Double>> topic =  word_distribute(c, wordSBTs); 

    
    BufferedWriter bwOut = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), "UTF8"));
    writingHelper(bwOut, root, 0, topic, dict);
    bwOut.close();
  }
  
  public static void showChars(String inFile, int targetChar) throws Exception {
    BufferedReader brIn = new BufferedReader(
        new InputStreamReader(new FileInputStream(inFile), "UTF8" ));
    String sLine;
    int ct = 0;
    while((sLine = brIn.readLine())!=null) {
      if(ct > targetChar - 100 && ct < targetChar + 100) {
        System.out.println("L: " + sLine);
      }
      if(ct <= targetChar && ct+sLine.length() >= targetChar) {
        System.out.println("problem line: ");
        System.out.println(sLine);
        System.out.println("problem char in center of: ");
        int lbound = Math.max(0, targetChar - (ct + 2));
        int rbound = Math.min(sLine.length(), targetChar - (ct - 2));
        System.out.println(sLine.substring(lbound, rbound));
      }
      ct += sLine.length() + 1;
    }
    brIn.close();
  }
  
//this function writes json
 public static void writingHelper(BufferedWriter bwOut, SparseBackoffTreeStructure curr, int indent, HashMap<Integer,TreeMap<Integer,Double>> topic,
     HashMap<Integer, String> dict) throws Exception {
   String ind = "";
   for (int i = 0; i < indent; i++){
     ind += "\t";
   }

   bwOut.write(ind + "{\"name\": " + "\""+ curr._minGlobalIndex + "-" + (curr._minGlobalIndex + curr.numLeaves() - 1) + "\"," + "\n");
   if(curr._children != null){
     bwOut.write(ind + "\"children\": " + "[\n");
     int lenOfChildren = curr._children.length;
     for(int i = 0; i < lenOfChildren - 1; i++){
       writingHelper(bwOut, curr._children[i], indent + 1, topic, dict);
       bwOut.write(",\n");
     }
     writingHelper(bwOut, curr._children[lenOfChildren - 1], indent + 1, topic, dict);
     bwOut.write("\n" + ind + "\t]");
   }
   // if we reach the second last layer
   else {
     bwOut.write(ind + "\"children\": " + "[\n");
     int len = curr._numLeaves.length;
     // output the topics
     int per = 3;
     for(int j = 0; j < len ; j++){
       int topNum = j + curr._minGlobalIndex;
       bwOut.write(ind + "\t{\"name\": " + "\""+ topNum + "\"");
       bwOut.write(",\n" + ind + "\t\"children\": " + "[\n");
       
       
       int idx = 0;
       if(topic.get(topNum) != null) {
         bwOut.write(ind + "\t\t{\"name\": \"");
         String vals = "";
         for(int k : topic.get(topNum).keySet()){
           if(idx > 0 && idx < per){
             bwOut.write(" " + dict.get(k));
             vals += topic.get(topNum).get(k) + " ";
           }
           else if(idx == per){
             break;
           }
           else{
             bwOut.write(dict.get(k));
             vals += topic.get(topNum).get(k) + " ";
           }
           idx ++;
         }
         bwOut.write("\"");
         bwOut.write(", \"value\": \"" + vals + "\"");
         bwOut.write("}\n" + ind);
       }
       bwOut.write("\t]");
       if(j < len - 1){
         bwOut.write("\n" + ind + "\t},\n");
       }
       else{
         bwOut.write("\n" + ind + "\t}\n");
       }
     }
           
     
     bwOut.write("\n" + ind + "\t]");
   }
   bwOut.write("\n");
   bwOut.write(ind + "}");
 }
 
 @SuppressWarnings("unchecked")
 public static  HashMap<Integer,TreeMap<Integer,Double>> word_distribute(Corpus c, SparseBackoffTree [] sbtToOutput) throws Exception{
 // Set up a hash table of the word distribution given topic
   @SuppressWarnings("rawtypes")
   HashMap<Integer,HashMap<Integer,Double>> word_givenTopic = new HashMap();
   for(int i = 0; i < c._pWord.length; i++) {
     if(c._pWord[i] > 0.0) {
       TIntDoubleHashMap hm = sbtToOutput[i].getLeafCounts();
       TIntDoubleIterator it = hm.iterator();
       while(it.hasNext()) {
         it.advance();
         int sId = it.key();
         double val = it.value()/c._pWord[i];
         if(word_givenTopic.containsKey(sId)){
                       HashMap<Integer,Double> sublist = word_givenTopic.get(sId);
                       sublist.put(i,val);
         }else{
           @SuppressWarnings("rawtypes")
           HashMap<Integer,Double> sublist = new HashMap();
           sublist.put(i,val);
           word_givenTopic.put(sId,sublist);
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
//    showChars(args[0], 337266);
    if(args.length != 4) {
      System.err.println("usage: SBTVisualizer <modelType=doc|seq> <model file> <output file> <dictionary file>");
      return;
    }
    else {
      outputWordToTopicFile(args[1], args[2], args[0], args[3]);
    }
  }

}
