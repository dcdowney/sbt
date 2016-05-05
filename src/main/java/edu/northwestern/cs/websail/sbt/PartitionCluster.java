package edu.northwestern.cs.websail.sbt;

import java.util.ArrayList;
import java.util.Set;
import java.util.TreeSet;
import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntDoubleHashMap;

public class PartitionCluster {
	
	public int id;
	public int token;
	public Set<Integer> wordSet;
	
	public PartitionCluster(int id)
	{
		this.id = id;
		wordSet = new TreeSet<Integer>();
		token = 0;
	}
	
	public void insert(TIntArrayList sentence){
		token += sentence.size();		
		TIntIterator it = sentence.iterator();
		while(it.hasNext()) {
			int i = it.next();
			wordSet.add(i);
		}
	}

}
