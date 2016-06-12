package edu.northwestern.cs.websail.sbt;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.util.ArrayList;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TIntArrayList;

public class Corpus implements Serializable {
	
	private static final long serialVersionUID = 1L;
	
	public int _VOCABSIZE = -1;  
	public int _NUMDOCS = -1;
	public long _NUMTOKENS = -1;
	TIntArrayList [] _docs; //the words
	TIntArrayList [] _z; //the topic variable assignments  //TODO: think about factoring out
	TIntArrayList [] _scratchZ; //accumulator for new zs
	double [] _pWord; //marginal word prob

	//returns number of docs
	public int readCorpusDat(String inFile, boolean updatePWord) throws Exception {
		BufferedReader brIn = new BufferedReader(
				new InputStreamReader(new FileInputStream(inFile), "UTF8" ));
		String sLine;
		int i=0;
		long toks = 0L;
		ArrayList<TIntArrayList> aldocs = new ArrayList<TIntArrayList>();
		if (updatePWord)
			_pWord = new double[_VOCABSIZE];
		int ct = 0;
		while((sLine = brIn.readLine())!=null) {
			TIntArrayList ll = lineToList(sLine);
			aldocs.add(ll);
			if (updatePWord) {
				TIntIterator it = ll.iterator();
				while(it.hasNext()) {
					_pWord[it.next()]++;
				}
			}
			toks += (long)ll.size();
			ct++;
		}
		brIn.close();
		System.out.println("Tokens: " + toks);
		System.out.println("Lines processed: " + ct + " Lines saved: " + aldocs.size());
		_NUMTOKENS = toks;
		double dubToks = (double)toks;
		if(updatePWord)
			for(int j=0; j<_pWord.length; j++)
				if(_pWord[j]>0.0)
					_pWord[j] /= dubToks;
		_NUMDOCS = aldocs.size();
		_docs = aldocs.toArray(new TIntArrayList[aldocs.size()]);
		return i;
	}
	
	public static TIntArrayList lineToList(String sLine) {
		String [] idFeats = sLine.split("\t");
		String feats;
		if(idFeats.length > 1) 
			feats = idFeats[1];
		else
			feats = idFeats[0];
		String [] featVals = feats.split(" ");
		TIntArrayList al = new TIntArrayList(featVals.length);
		for(int i=0; i<featVals.length; i++) {
			String [] featVal = featVals[i].split(":");
			int val;
			if(featVal.length==1) { //implicit 1
				val = 1;
			}
			else{
				 val = Integer.parseInt(featVal[1]);
			}
			int feat = Integer.parseInt(featVal[0]);
			for(int j=0; j<val; j++) {
				al.add(feat);
			}
		}
		return al;
	}

	
	//HACK: keeps word distributions, but re-sets docs, counts, and z's for new corpus
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
		}
		

}
