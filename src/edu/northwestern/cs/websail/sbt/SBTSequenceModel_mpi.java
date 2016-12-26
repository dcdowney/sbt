package edu.northwestern.cs.websail.sbt;

import mpi.*;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.LinkedList;
import java.util.Set;
import com.javamex.classmexer.MemoryUtil;
import com.javamex.classmexer.MemoryUtil.VisibilityFilter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.lang.reflect.Method;
import java.net.InetAddress;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntDoubleHashMap;

//TODO: factor out all the gross duplication between this class and SBTDocumentModel

public class SBTSequenceModel_mpi implements Serializable {

	private static final long serialVersionUID = 2L;

	private class GibbsDoer implements Runnable {
    	int _start;
    	int _end;
    	long _changes = 0;
    	
    	public void run() {
    		_changes = sampleRange(_start, _end);
    	}
    }
	
	private class DebugDoer implements Runnable {
    	int _start;
    	int _end;
    	long result = 0;
    	
    	public void run() {
    		System.out.println("Start: " + _start + " End : " + _end);
    		result = DebugRun(_start, _end);
    	}
    }
	
    public long DebugRun(int start, int end) {
    	
    	long chg = 0;
    	
    	long result = 0;
		for(int i=start; i<=end; i++) {
				for(int k = 0; k< 1500; k++)
					for(int l = 0;l<100;l++)
						result += (i+k+l);
		}

		return chg;
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

    public static class TestGetProbsDoer extends TestEnsembleExactDoer {
    	double [][] _probs;
    	
    	public TestGetProbsDoer(int i, double [][] wordTopicMarginal, SBTSequenceModel_mpi [] models, double [][][] tm, double [] ws) {
    		super(i, wordTopicMarginal, models, tm, ws);
    	}
    	
    	@Override
    	public void run() {
    		_probs = getEnsembleProbsForDoc(_doc, _wordTopicMarginal, _sbtsm, _tm);
    	}
    }
    
    private class TestExactDoer implements Runnable {
     	int _doc;
    	double [] _wordTopicMarginal;
    	double [] _res;
    	double [][] _tm;
    	
    	public TestExactDoer(int i, double [] wordTopicMarginal, double [][] transModel) {
    		_doc = i;
    		_wordTopicMarginal = wordTopicMarginal;
    		_tm = transModel;
    	}
    	
    	public void run() {
    		//_res = testOnDocFullPpl(_doc, _wordTopicMarginal);
    		_res = testOnDocFullPplExact(_doc, _wordTopicMarginal, _tm);
    		System.out.println(_doc + "\t" + _res[0] + "\t" + _res[1]);
    	} 	
    }
    
    public static class TestEnsembleDoer implements Runnable {
    	int _doc;
    	double [][] _wordTopicMarginal;
    	SBTSequenceModel_mpi [] _sbtsm;
    	double [] _res;
    	
    	public TestEnsembleDoer(int i, double [][] wordTopicMarginal, SBTSequenceModel_mpi [] models) {
    		_doc = i;
    		_wordTopicMarginal = wordTopicMarginal;
    		_sbtsm = models;
    	}
    	
    	public void run() {
    		_res = testEnsembleOnDocFullPpl(_doc, _wordTopicMarginal, _sbtsm);
    		System.out.println(_doc + "\t" + _res[0] + "\t" + _res[1]);
    	}
    }
	
    public static class TestEnsembleExactDoer implements Runnable {
    	int _doc;
    	double [][] _wordTopicMarginal;
    	double [][][] _tm;
    	SBTSequenceModel_mpi [] _sbtsm;
    	double [] _res;
    	double [] _ws;
    	double [] _ps;
    	
    	public TestEnsembleExactDoer(int i, double [][] wordTopicMarginal, SBTSequenceModel_mpi [] models, double [][][] tm,
    			double [] ws) {
    		_doc = i;
    		_wordTopicMarginal = wordTopicMarginal;
    		_sbtsm = models;
    		_tm = tm;
    		_ws = ws;
    	}
    	
    	public void run() {
    		_ps = new double[_sbtsm[0]._c._docs[_doc].size()];
    		_res = testEnsembleOnDocExact(_doc, _wordTopicMarginal, _sbtsm, _tm, _ws, _ps); //side effect, fills _ps with probs per tok
    		System.out.println(_doc + "\t" + _res[0] + "\t" + _res[1]);
    	}
    }
    
	//Model information:
    //state index zero is special start-of-sentence state
	int [] _branchingFactors;
	SparseBackoffTree [] _forward;
	SparseBackoffTree _startState;
	SparseBackoffTree _endState;
	SparseBackoffTree [] _backward;
	SparseBackoffTree [] _wordToState;
	double [] _wordStateMarginal;
	double [] _backwardStateMarginal;
	private LinkedList<Integer> _changesList;
	private int listLength = 40;
	private double listTrends = 0.40;		//increase percent  0.1 = 10%;
	private double listThreshold = 0.03;	//within,   5% percent (value)
	private int buffSize = 100;
	private static int _wordStoreFlag[];		// if -1: word doesn't exist in this mpi process; other wise, the rank of the mpi process aggregate the sbt
	private static int _wordReduceFlag[];		// the rank of the mpi process aggregate the sbt
	private double gibbsChunkSize = 0.05;	//percentage of gibbs partial update. 0.05 = update 5% each sampling 
	private double gibbsNextPointer = 0.0;	//next gibbs sampling start from gibbsNextPointer% data. consecutive sampling.
	
	SparseBackoffTreeStructure _struct;
	double [] _forwardDelta = null;
	double [] _backwardDelta = null;
	double [] _wordDelta = null;
	double MaxTokenThreshold = 0.08;
	
	Corpus _c;
	
	double [][] baseProbs = null; //holds existing model probabilities when training a layer of an ensemble
	
	Random _r = new Random(4);
	
	//TODO: think about wrapping all of threading/train config up into a trainer class
	
	//threading:
	private int _NUMTHREADS = -1;
	private int[] _THREADBREAKS = null; //inclusive doc indices where threads *end* (initial thread starts
							//at start of processbreak implicitly)
	private int[] _PROCESSBREAKS = null; //inclusive doc indices where mpi-process *end* (initial thread starts
	//at 0 implicitly)
	
	private int _MPIRANK = -1;
	private int _MPINUM = -1;
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
	
	public SBTSequenceModel_mpi(String configFile) throws Exception {
		_c = new Corpus();
		readConfigFile(configFile);
		if(_USEEXPANSION == 1) {
			if(_EXPANSIONSCHEDULE == null) {
				_EXPANSIONSCHEDULE = defaultExpansion(_branchingFactors.length);
			}
			_expansionBranchFactors = Arrays.copyOf(_branchingFactors, _branchingFactors.length);
			_branchingFactors = Arrays.copyOf(_branchingFactors, _EXPANSIONSCHEDULE[0]+1);
		}
		_struct = new SparseBackoffTreeStructure(_branchingFactors);
		initDeltas();
	}
	
	public void initDeltas() {
		_forwardDelta = new double[_branchingFactors.length];
		_backwardDelta = new double[_branchingFactors.length];
		_wordDelta =new double[_branchingFactors.length];
		double initDelta = 0.9 / (double)_branchingFactors.length;
		Arrays.fill(_wordDelta, initDelta);
		Arrays.fill(_forwardDelta, initDelta);
		Arrays.fill(_backwardDelta, initDelta);
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
	
	private void setThreadBreaks(int from) throws MPIException{
		
		if(from==0)
		{
			// random split
			/*
			TIntArrayList t;
			_THREADBREAKS = new int[_NUMTHREADS];
			for(int i=0;i<_c._NUMDOCS;i++)
			{
				int index = _r.nextInt(i + 1);
				t = _c._docs[i];
				_c._docs[i] = _c._docs[index];
				_c._docs[index] = t;
			}
			long approxToks = _c._NUMTOKENS / _NUMTHREADS;
			long ct = 0;
			int thread = 0;
			for(int i=0; i<_c._NUMDOCS; i++ ) {
				ct += _c._docs[i].size();
				if(ct > approxToks) {
					_THREADBREAKS[thread++] = i;
					ct = 0;
				}
			}
			//extra goes in last thread:
			_THREADBREAKS[_NUMTHREADS - 1] = _c._NUMDOCS - 1;
			*/

			//Jaccard split
			
			_THREADBREAKS = new int[_NUMTHREADS];
			int assign[] = new int[_c._NUMDOCS];
			
			int rank = MPI.COMM_WORLD.getRank();
			int size = MPI.COMM_WORLD.getSize();
			
			PartitionCluster clusters[] = new PartitionCluster[size];
			_PROCESSBREAKS = new int[size];
			
			if(rank==0){
				for(int i=0;i<_c._NUMDOCS;i++)
					assign[i] = -1;
				
				for(int i=0;i<size;i++)
					clusters[i] = new PartitionCluster(i);
				
				long tokenPerCluster = 0;
				
				for(int i=0;i<_c._NUMDOCS;i++)
					tokenPerCluster += _c._docs_bk[i].size();
				tokenPerCluster /= size;
				tokenPerCluster += 10;
				
				long MaxTokenLimit = (long) (tokenPerCluster * (1.0 + MaxTokenThreshold));
				
				for(int i=0;i<_c._NUMDOCS;i++)
				{
					// get sentence 
					TIntArrayList sentence = _c._docs_bk[i];
					int clusterId = -1;
					for(int j=0;j<size;j++)
					{
						double jaccardIndex = 0.0;
						double tempJaccard = 0.0;	
						long union = -1;
						
						if(clusters[j].token <= MaxTokenLimit)
						{
							//select a cluster to assign
							Set<Integer> wordSet = clusters[j].wordSet;
							union = 0;
							TIntIterator it = sentence.iterator();
							while(it.hasNext()) {
								int num = it.next();
								if(!wordSet.contains(num))
									union ++;
							}
							long intersection =  sentence.size() - union;
							tempJaccard = (double)intersection / union;
							if( tempJaccard >= jaccardIndex ){
								jaccardIndex = tempJaccard;
								clusterId = j;
							}
						}
					}
					
					if(clusterId==-1){
						System.out.println("why -1 ? strange");
						System.out.println("max token size = " + MaxTokenLimit);
						for(int j=0;j<size;j++)
							System.out.println("cluster : " + j + " has token " +  clusters[j].token);
					}
					
					clusters[clusterId].insert(sentence);
					assign[i] = clusterId;
				}
				
				int p = 0;
				for(int i=0;i<size;i++)
				{
					for(int j=0;j<_c._NUMDOCS;j++)
					{
						if(assign[j]==i)
						{		
							TIntArrayList t;
							t = _c._docs_bk[j];
							_c._docs_bk[j] = _c._docs_bk[p];
							_c._docs_bk[p] = t;
							
							int t_int;
							
							t_int = assign[j];
							assign[j] = assign[p];
							assign[p] = t_int;
							
							p++;
						}
					}
					_PROCESSBREAKS[i] = p;
				}
				_PROCESSBREAKS[size - 1] = _c._NUMDOCS;
			}
			
			MPI.COMM_WORLD.bcast(_PROCESSBREAKS, _PROCESSBREAKS.length, MPI.INT, 0);
			
			int local_total_length = 0;
			int local_arrays_lengths[] = null;
			int l = 0;
			if(rank==0)
				l = _PROCESSBREAKS[rank];
			else
				l = _PROCESSBREAKS[rank] - _PROCESSBREAKS[rank-1];
			local_arrays_lengths = new int[l];
			
			int buffSizes[] = new int[size];
			
			if(rank==0){
				for(int i=0;i<size;i++)
					buffSizes[i] = 0;
				for(int i=0;i<_c._NUMDOCS;i++){
					buffSizes[assign[i]] += _c._docs_bk[i].size();
				}
			}
			
			MPI.COMM_WORLD.bcast(buffSizes, size, MPI.LONG, 0);
			
			for(int i = 0;i<buffSizes.length;i++)
				assert(buffSizes[i]>=0);
			
			
			
			int recvBuff[] = new int[buffSizes[rank]];
			
			// set recev buff

			for(int i=0;i<size;i++){
				if(i==0){
					if(rank==0){
						int p=0;
						for(int j=0;j<_PROCESSBREAKS[rank];j++){
							int temp[] = _c._docs_bk[j].toArray();
							local_arrays_lengths[j] = temp.length;
							for(int k=0;k<temp.length;k++)
								recvBuff[p++] = temp[k];
						}
					}
				}
				else{
					
					if(rank==0){
						int send_local_arrays_lengths[] = new int[_PROCESSBREAKS[i] - _PROCESSBREAKS[i-1]];
						int sendBuff[] = new int[buffSizes[i]];
						int p = 0;
						for(int j=_PROCESSBREAKS[i-1];j<_PROCESSBREAKS[i];j++){
							send_local_arrays_lengths[j-_PROCESSBREAKS[i-1]] = _c._docs_bk[j].size();
							int temp[] = _c._docs_bk[j].toArray();
							for(int k=0;k<temp.length;k++)
								sendBuff[p++] = temp[k];
						}
						MPI.COMM_WORLD.send(sendBuff, sendBuff.length, MPI.INT, i, 0);
						MPI.COMM_WORLD.send(send_local_arrays_lengths, send_local_arrays_lengths.length, MPI.INT, i, 1);
					}
					else{
						if(i==rank){
							MPI.COMM_WORLD.recv(recvBuff, recvBuff.length, MPI.INT, 0, 0);
							MPI.COMM_WORLD.recv(local_arrays_lengths, local_arrays_lengths.length, MPI.INT, 0, 1);
						}
					}
				}
			}
			
			//split recvbuff to _doc & split threads
			
			int p = 0;
			int l_p = 0;
			ArrayList<TIntArrayList> aldocs = new ArrayList<TIntArrayList>();
			for(l=0;l<_c._NUMDOCS;l++){
				if(!within_doc_range(l)){
					aldocs.add(null);
					continue;
				}
				// why out of bound
				int temp[] = new int[local_arrays_lengths[l_p++]];
				for(int j=0;j<temp.length;j++)
					temp[j] = recvBuff[p++];
				TIntArrayList ll = new TIntArrayList();
				ll.addAll(temp);
				aldocs.add(ll);
			}
			
			_c._docs = aldocs.toArray(new TIntArrayList[aldocs.size()]);
			
			if(rank==0)
				_c._docs_bk = null;
			
			long approxToks = buffSizes[rank] / _NUMTHREADS;
			long ct = 0;
			int thread = 0;
			int start = -1;
			int end = -1;
			
			if(rank==0)
				start = 0;
			else
				start = _PROCESSBREAKS[rank-1];
			
			end = _PROCESSBREAKS[rank];
			
			for(int i=start; i<end; i++ ) {
				ct += _c._docs[i].size();
				if(ct > approxToks) {
					_THREADBREAKS[thread++] = i;
					ct = 0;
				}
			}
			//extra goes in last thread:
			_THREADBREAKS[_NUMTHREADS - 1] = _PROCESSBREAKS[rank];
			
			/*
			_THREADBREAKS = new int[_NUMTHREADS];
			long approxToks = _c._NUMTOKENS / _NUMTHREADS;
			long ct = 0;
			int thread = 0;
			for(int i=0; i<_c._NUMDOCS; i++ ) {
				ct += _c._docs[i].size();
				if(ct > approxToks) {
					_THREADBREAKS[thread++] = i;
					ct = 0;
				}
			}
			//extra goes in last thread:
			_THREADBREAKS[_NUMTHREADS - 1] = _c._NUMDOCS - 1;
			*/
		}
		else
		{
			_THREADBREAKS = new int[_NUMTHREADS];
			long approxToks = _c._NUMTOKENS / _NUMTHREADS;
			long ct = 0;
			int thread = 0;
			for(int i=0; i<_c._NUMDOCS; i++ ) {
				ct += _c._docs[i].size();
				if(ct > approxToks) {
					_THREADBREAKS[thread++] = i;
					ct = 0;
				}
			}
			//extra goes in last thread:
			_THREADBREAKS[_NUMTHREADS - 1] = _c._NUMDOCS - 1;
		}
		
		
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
		
		BufferedReader brIn = new BufferedReader(
				new InputStreamReader(new FileInputStream(inFile), "UTF8" ));
		String sLine;
		Method [] ms = this.getClass().getMethods();
		
		
		int rank = MPI.COMM_WORLD.getRank();
		
		MPI.COMM_WORLD.barrier();
		long stTime = System.currentTimeMillis();
		if(rank==0)
			System.out.println("reading config file.");
		
		while((sLine = brIn.readLine())!=null) {
			if(sLine.startsWith("#"))
				continue;
			if(rank==0)
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
				if(rank==0)
					System.out.println(m.getName() + "\t" + arg);
				m.invoke(this, arg);
			}
			else if(c.isArray() && c.getComponentType().getName().equals("int")) {
				String [] vals = fields[1].split(" ");
				int [] arg = new int[vals.length];
				for(int i=0; i<arg.length; i++) {
					arg[i] = Integer.parseInt(vals[i]);
				}
				if(rank==0)
					System.out.println(m.getName() + "\t" + Arrays.toString(arg));
				m.invoke(this, arg);
			}
		}
		brIn.close();
		
		MPI.COMM_WORLD.barrier();
		stTime = System.currentTimeMillis() - stTime;
		if(rank==0)
			System.out.println("finish reading -------reading time = " + stTime);
	}

	private double [][] readBaseProbs(String inFile) throws Exception {
		double [][] out = new double[_c._docs.length][];
		BufferedReader brIn = new BufferedReader(
				new InputStreamReader(new FileInputStream(inFile), "UTF8" ));
		String sLine;

		for(int i=0; i<out.length; i++) {
			sLine = brIn.readLine();
			String [] sps = sLine.split(" ");
			double [] ps = new double[_c._docs[i].size()];
			for(int j=0; j<ps.length; j++) {
				ps[j] = Double.parseDouble(sps[j]);
			}
			out[i] = ps;
		}
		brIn.close();
		return out;
	}
	
	private static void writeBaseProbs(String outFile, double [][] ps) throws Exception {
		BufferedWriter bwOut = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "UTF8"));
		for(int i=0; i<ps.length; i++) {
			for(int j=0; j<ps[i].length;j++) {
				if(j>0)
					bwOut.write(" ");
				bwOut.write(ps[i][j] + "");
			}
			bwOut.write("\r\n");
		}
		bwOut.close();
	}
	
	private double [][] transModel() {
		int numStates = this._wordStateMarginal.length;
		double [][] out = new double[numStates][numStates];
		double diff = 0.0;
		double total = 0.0;
		double [] nextStateMarginal = new double[numStates];
		double [] invBackwardMarginal = new double[numStates];
		for(int i=0; i<invBackwardMarginal.length; i++) {
			invBackwardMarginal[i] = 1.0 / _backwardStateMarginal[i];
		}
		for(int i=0; i<_backward.length; i++) {
			_backward[i].divideCountsBy(invBackwardMarginal);
		}
		//first compute estimates based on forward only:
		for(int i=0; i<numStates; i++) {
			for(int j=0; j<numStates;j++) {
				double forwardEst = 0.0;
				double countMass = _forward[i].getSmoothed(j); 
				if(_forward[i]._totalMass > 0.0)
					forwardEst = countMass /_forward[i]._totalMass;
				else
					forwardEst = 1.0/(double)numStates;
				out[i][j] = forwardEst;
				nextStateMarginal[j] += countMass;
				total += countMass;
			}
		}
		total /= (double)numStates;
		for(int i=0; i<numStates; i++) {
			nextStateMarginal[i] /= total;
		}
		
		double [][] be = new double [numStates][numStates];
		for(int i=0; i<numStates; i++) {
			double jSum = 0.0;
			for(int j=0; j<numStates;j++) {
				double backwardEst = 0.0;
				if(_backward[j]._totalMass > 0.0)
					//backwardEst = nextStateMarginal[j] * _backward[j].getSmoothed(i) / _backward[j]._totalMass;
					backwardEst = _backward[j].getSmoothed(i);
				else
					backwardEst = _backwardStateMarginal[i] / (double)numStates; //not sure about this
				be[i][j] = backwardEst;
				jSum += backwardEst;
			}
			for(int j=0; j<numStates;j++) {
				be[i][j] /= jSum;
			}
		}
		total = 0.0;
		for(int i=0; i<numStates; i++) {
			for(int j=0; j<numStates;j++) {
				double forwardEst = out[i][j];
				double backwardEst = be[i][j];
				diff += Math.abs(forwardEst - backwardEst);
				if(Math.abs(forwardEst - backwardEst) > 0.001) {
					System.out.println("diff on " +i + "," + j + " is " + Math.abs(forwardEst - backwardEst));
				}
				backwardEst = forwardEst;
				out[i][j] = (forwardEst + backwardEst)/2.0;
				total += out[i][j];
			}
		}
		System.out.println("trans model diff: " + diff);
		System.out.println("trans model total: " + total);
		return out;
	}
	

	protected TIntDoubleHashMap aggregateCounts(TIntArrayList occurs) {
		System.out.println("call aggregate counts !!!!  delete this output imme");
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
		TIntArrayList zs = _c._z[docId];
		TIntArrayList scratchZs = _c._scratchZ[docId];
		TIntArrayList doc = _c._docs[docId];
		//System.out.println("adding " + docCts.toString());

		for(int i=0; i<doc.size(); i++) {
			int newZ = sampleZ(docId, i, false);
			if(newZ != zs.get(i))
				changes++;
			scratchZs.set(i, newZ);
		}
		return changes;
	}
	
	public int sampleZ(int doc, int pos, boolean testing) {
		int w = _c._docs[doc].get(pos);
		int curZ = _c._z[doc].get(pos);
		SparseBackoffTree sbtForward = null;
		if(pos==0) {
			sbtForward = this._startState;
		}
		else {
			if(!testing && _SAMPLEONSAMPLES) { //HACK:we must sample left-to-right for this to work
				sbtForward = this._forward[_c._scratchZ[doc].get(pos-1)];
			}
			else 
				sbtForward = this._forward[_c._z[doc].get(pos-1)];
		}
		SparseBackoffTree sbtBackward = null;
		boolean forward = (testing && (pos==_c._docs[doc].size() - 1 || _c._z[doc].get(pos+1)==-1)); 
		if(!forward) { //don't use backward distribution if doing forward prediction
			if(pos==_c._docs[doc].size() - 1) {  
				sbtBackward = this._endState;
			}
			else {
				int nextZ = _c._z[doc].get(pos+1);
				if(nextZ >= 0) {
					sbtBackward = this._backward[nextZ];
				}
			}
		}
		SparseBackoffTree sbtWord = this._wordToState[w];
		boolean unseen = false;
		//if unseen, using simple smoothing distribution for the purpose of sampling:
		if(sbtWord==null)
			System.out.println(" null word = " +w);
		if(sbtWord._totalMass == 0.0) {
			unseen = true;
		}
		double bMarginalInv = 1.0f;
		double wMarginalInv = 1.0f;
		if(curZ >= 0) {
			bMarginalInv = 1.0 / this._backwardStateMarginal[curZ];
			wMarginalInv = 1.0 / this._wordStateMarginal[curZ];
		}
		
		SparseBackoffTree [] sbts;
		double [] subs;
		if(!forward) {
			if(unseen) {
				sbts = new SparseBackoffTree [] {sbtForward, sbtBackward};
				subs = new double [] {1.0, wMarginalInv};
			}
			else {
				sbts = new SparseBackoffTree [] {sbtWord, sbtForward, sbtBackward};
				subs = new double [] {wMarginalInv, 1.0, bMarginalInv};
			}
		}
		else {
			if(unseen) {
				sbts = new SparseBackoffTree [] {sbtForward};
				subs = new double [] {1.0};
			}
			else {
				sbts = new SparseBackoffTree [] {sbtForward, sbtWord};
				subs = new double [] {1.0, wMarginalInv};
			}
		}
		if(testing) //don't subtract at test time:
			subs = null;
		if(!testing && _SAMPLEONSAMPLES && pos > 0) {
			if(_c._z[doc].get(pos-1)==_c._scratchZ[doc].get(pos-1))
				subs[0] = 0.0;
		}
		SparseBackoffTreeIntersection sbti;
		if(subs!=null) 
			sbti = new SparseBackoffTreeIntersection(sbts, curZ, subs, true);
		else
			sbti = new SparseBackoffTreeIntersection(sbts, true);
		
		if(doc==2000000 && pos == 2){
    		System.out.println("doc: " + doc  + " pos: " + pos + " word: " + w + " curZ: " + curZ);
    		for(int i=0;i<sbts.length;i++){
    			System.out.println("SBT " + i + "---------");
    			System.out.println(Arrays.toString(sbts[i].toDoubleArray()));
    		}
		}
		
		
		int sample = 0;
		if(baseProbs !=null && !unseen && !forward) {
			double wordOutProb = sbti._totalMass / sbti._siblingSbtis[0]._totalMass;
			wordOutProb *= (double)this._struct.numLeaves();
			wordOutProb /= this.get_NUMTOKENS();
			double rng = wordOutProb + baseProbs[doc][pos];
			double choice = _r.nextDouble() * rng;
			if(choice < wordOutProb)
				sample = sbti.sample(_r);
			else {
				sample = sbti._siblingSbtis[0].sample(_r); //note -- dangerous, depends on this siblingSbti data structure
					//containing the intersection of forward and backward distributions
				_numSkipWords++;
			}
		}
		else
			sample = sbti.sample(_r);
		return sample;
	}
	
    private long sampleRange(int start, int end) {
    	
    	long chg = 0;

		for(int i=start; i<end; i++) {
			
			if((double)(i-start)/(double)(end-start)>=gibbsNextPointer && (double)(i-start)/(double)(end-start)<(gibbsNextPointer+gibbsChunkSize)){
				chg += (long)sampleDoc(i);
			}
			else{
				TIntArrayList scratchZs = _c._scratchZ[i];
				int tArray[] = _c._z[i].toArray();
				for(int j=0;j<tArray.length;j++)
					scratchZs.set(j, tArray[j]);
			}
				
		}

		return chg;
    }	
	
    //Zheng added: incFile
	public void initZ(TIntArrayList [] ds, int maxVal, boolean random, boolean scratch, String incFile) throws Exception {
		TIntArrayList [] zs;
		if(scratch) {
			_c._scratchZ = new TIntArrayList[ds.length];
			zs = _c._scratchZ;
		}
		else {
			_c._z = new TIntArrayList[ds.length];
			zs = _c._z;
		}
		
		if(incFile==null){
			for(int i=0; i<ds.length; i++) {
				if(!within_doc_range(i))
				{
					continue;
				}
				zs[i] = new TIntArrayList(ds[i].size());
				for(int j=0; j<ds[i].size(); j++)
				{
					zs[i].add(random ? _r.nextInt(maxVal) : 0); //don't use zero, reserved for start/end
//					zs[i].add(random ? 2 : 0); //don't use zero, reserved for start/end
				}
			}
		}
		else{
			SBTSequenceModel_mpi t_sbtsm = readModel(incFile);
			// copy explicitly
			
			System.out.println("in file the length of z is " + t_sbtsm._c._z.length);
			System.out.println("in restart the length of z is " + ds.length);
			
			for(int i=0; i<ds.length; i++) {
				if(!within_doc_range(i))
				{
					continue;
				}
				
				zs[i] = new TIntArrayList(ds[i].size());
				for(int j=0; j<ds[i].size(); j++)
					zs[i].add(t_sbtsm._c._z[i].get(0));
			}
		}
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
	public TIntDoubleHashMap [] aggregateCounts(int from, TIntArrayList [] zs, TIntArrayList [] ws, int partialUpdate) throws MPIException{
		TIntDoubleHashMap [] hm;
		
		//handle start/end as special case
		if(from==2 || from==-2) { 
			hm = new TIntDoubleHashMap[1];
			hm[0] = new TIntDoubleHashMap();
			for(int i=0; i<zs.length; i++) {
				if(!within_doc_range(i) && (partialUpdate==1))
				{
					continue;
				}
				if(from==-2) {
					hm[0].adjustOrPutValue(zs[i].get(0), 1.0, 1.0);
				}
				else {
					hm[0].adjustOrPutValue(zs[i].get(zs[i].size()-1), 1.0, 1.0);
				}
			}
		}
		else{
			if(from==0)
				hm = new TIntDoubleHashMap[_c._VOCABSIZE];
			else
				hm = new TIntDoubleHashMap[_struct.numLeaves()];
			
			//forward steps from word j to j+1
			//backward steps from j+1 to j
			//word steps to j
			for(int i=0; i<zs.length; i++) {
				if(!within_doc_range(i) && (partialUpdate==1))
				{
					continue;
				}
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
					else if(from==0){// modified structure !! zheng
						xID = ws[i].get(j);
						z = zs[i].get(j);
					}
					else
					{
						xID = ws[i].get(j);
						z = zs[i].get(j);
					}
					if(hm[xID] == null)
						hm[xID] = new TIntDoubleHashMap();
					hm[xID].adjustOrPutValue(z, 1.0, 1.0);
				}
			}
		}
		
		
		/*
		// reduce hm
		int totalNum = _struct.numLeaves() * hm.length;
		
		double hmArray[] = new double[totalNum];
		
		int pos = 0;
		for(int i=0;i<hm.length;i++){
			if(hm[i]==null){
				for(int j=0;j<_struct.numLeaves();j++){
					hmArray[pos++] = 0.0;
				}
			}
			else{
				for(int j=0;j<_struct.numLeaves();j++){
					hmArray[pos++] = hm[i].get(j);
				}
			}
		}
		
		MPI.COMM_WORLD.allReduce(hmArray, hmArray.length, MPI.DOUBLE, MPI.SUM);

		TIntDoubleHashMap [] newHm;
		
		newHm = new TIntDoubleHashMap[hm.length];
		
		pos = 0;
		
		for(int i=0;i<hm.length;i++){
			
			for(int j=0;j<_struct.numLeaves();j++){
				if(newHm[i]==null && hmArray[pos]!=0){
					newHm[i] = new TIntDoubleHashMap();
				}
				if(hmArray[pos]!=0)
					newHm[i].put(j, hmArray[pos]);
				pos++;
			}
		}
		*/
		
		double hmArray[] = new double[_struct.numLeaves() * buffSize];
		TIntDoubleHashMap [] newHm;	
		newHm = new TIntDoubleHashMap[hm.length];
		
		for(int i = 0;i< hm.length;){
			int len = buffSize;
			if(i+buffSize>hm.length)
				len = hm.length - i;
			
			int pos = 0;
			for(int j=0;j<len;j++){
				if(hm[i+j]==null){
					for(int k=0;k<_struct.numLeaves();k++){
						hmArray[pos++] = 0.0;
					}
				}
				else{
					for(int k=0;k<_struct.numLeaves();k++){
						hmArray[pos++] = hm[i+j].get(k);
					}
				}
			}
			
			MPI.COMM_WORLD.allReduce(hmArray, len*_struct.numLeaves(), MPI.DOUBLE, MPI.SUM);
			
			pos = 0;
			
			for(int j=0;j<len;j++){
				if(from==0 && partialUpdate==1 && _wordStoreFlag[i+j]==-1){
					pos += _struct.numLeaves();
				}
				else{
					for(int k=0;k<_struct.numLeaves();k++){
						if(newHm[i+j]==null && hmArray[pos]!=0.0){
							newHm[i+j] = new TIntDoubleHashMap();
						}
						if(hmArray[pos]!=0.0)
							newHm[i+j].put(k, hmArray[pos]);
						pos++;
					}
				}
			}
			
			
			i+=len;
		}
		
		return newHm;
	}
	
	private boolean within_doc_range(int a) throws MPIException {
		int rank=MPI.COMM_WORLD.getRank();
		int size=MPI.COMM_WORLD.getSize();
		int start = -1;
		int end = -1;
		if(rank==0)
			start = 0;
		else
			start = this._PROCESSBREAKS[rank-1];
		
		end = this._PROCESSBREAKS[rank];
		
		if(start<=a&&end>a)
			return true;
		else
			return false;
	}
	
	/**
	 * Returns sbts with given discounts of given type based on current sampler state
	 * if from = -2 then returns start-state parameteres
	 * if from = -1 then returns forward parameters
	 * if from = 0 then returns word-to-state parameters
	 * if from = 1 then backward parameters
	 * if from = 2 then returns end-state parameters
	 * @param dsWords	The discounts to use in the returned parameters
	 * @param from	see above
	 * @return
	 */
	public SparseBackoffTree [] getParamsFromZs(double [] dsWords, int from, TIntArrayList [] zs,
			TIntArrayList [] ws, int partialUpdate) throws MPIException{
		//aggregate:
		// if from ==0 make a word map to reduce the parameter space
		
		TIntDoubleHashMap [] hm = aggregateCounts(from, zs, ws, partialUpdate);
		
		// null hm entrires to reduce parameter space
		SparseBackoffTree [] out = new SparseBackoffTree[hm.length];
		
		for(int i=0; i<hm.length; i++) {			
			if(hm[i]==null && i!=0 && from==0)
				out[i] = null;
			else
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
	
	/**
	 * reads the corpus and initializes zs and model
	 * @param inFile
	 * @param maxVal
	 * @param incFile
	 * @return
	 * @throws Exception
	 */
	//Zheng added: incFile: null, if init; not null if restart
	public int initializeForCorpus(String inFile, int maxVal, String incFile) throws Exception {
		
		int rank = MPI.COMM_WORLD.getRank();
		MPI.COMM_WORLD.barrier();
		if(rank==0)
			System.out.println("before reading corpus----");
		int toks = _c.readCorpusDat(inFile, true);
		
		setThreadBreaks(0);
		
		this._wordStoreFlag = new int[_c._VOCABSIZE];
		this._wordReduceFlag = new int[_c._VOCABSIZE];
		
		for(int i=0;i<_c._VOCABSIZE;i++){
			_wordStoreFlag[i] = -1;
			_wordReduceFlag[i] = -1;
		}
		
		if(rank==0)
			System.out.println("doc length = " + _c._docs.length);
		
		for(int i=0;i<_c._docs.length;i++){
			if(within_doc_range(i)){
				for(int j=0;j<_c._docs[i].size();j++){
					_wordStoreFlag[_c._docs[i].get(j)] = 0;
					_wordReduceFlag[_c._docs[i].get(j)] = rank;
				}
			}
		}
		
		MPI.COMM_WORLD.allReduce(_wordReduceFlag, _wordReduceFlag.length, MPI.INT, MPI.MAX);
		if(rank==0)
			System.out.println("before InitZ");
		initZ(_c._docs, maxVal, true, false, incFile);
		updateModel(_c._z,1);
		getWordToStateSize();
		
		return toks;
	}
	
	private double[] getWordToStateSize() throws Exception{
		int rank = MPI.COMM_WORLD.getRank();
		int size = MPI.COMM_WORLD.getSize();
		long wordToStateSize = 0;
		long tLong[] = {0};
		double r[] = new double[]{-1.0,-1.0};
		
		for(int k=0;k<_wordToState.length;k++)
			wordToStateSize += SparseBackoffTree.SparseBackoffTreeSize(_wordToState[k]);
		tLong[0] = wordToStateSize;
		MPI.COMM_WORLD.allReduce(tLong, 1, MPI.LONG, MPI.SUM);
		r[0] = (double)(tLong[0]/1024/1024)/1024.0/(double)size;
//		if(rank==0)
//			System.out.println("ave wordToStateSBT size is" + r[0] + " GB ");
		
		tLong[0] = wordToStateSize;
		MPI.COMM_WORLD.allReduce(tLong, 1, MPI.LONG, MPI.MAX);
		r[1] = (double)(tLong[0]/1024/1024)/1024.0;
//		if(rank==0)
//			System.out.println("largest wordToStateSBT size is" + r[1] + " GB ");
		
		return r;
	}
	
	private double[] getUsedMemSize()throws Exception{
		int rank = MPI.COMM_WORLD.getRank();
		int size = MPI.COMM_WORLD.getSize();
		double[] r = new double[]{-1.0,-1.0};
		long tLong[] = {-1};
		long ttLong[] = {-1};
		
		Runtime rt = Runtime.getRuntime();
		rt.gc();
		
		tLong[0] = rt.totalMemory()-rt.freeMemory();
		
		MPI.COMM_WORLD.allReduce(tLong, ttLong, 1, MPI.LONG, MPI.SUM);
		r[0] = (double)(ttLong[0]/1024/1024)/1024.0/(double)size;
//		if(rank==0)
//			System.out.println("ave memory usage is " + r[0] + " GB per MPI_process");
		
		MPI.COMM_WORLD.allReduce(tLong, ttLong, 1, MPI.LONG, MPI.MAX);
		r[1] = (double)(ttLong[0]/1024/1024)/1024.0;
//		if(rank==0)
//			System.out.println("max memory usage is " + r[1] + " GB");
		
		return r;
	}
	
    private long gibbsPass() throws MPIException{
    	long changes= 0 ;
	    GibbsDoer [] gds = new GibbsDoer[_NUMTHREADS];
	    long stTime = System.currentTimeMillis();
	    //MPI part
		int rank = MPI.COMM_WORLD.getRank();
		int size = MPI.COMM_WORLD.getSize();
		for(int i=0; i<_NUMTHREADS;i++) {
			gds[i] = new GibbsDoer();
			if(rank==0)
			{
				if(i==0)
					gds[i]._start = 0;
				else
					gds[i]._start = _THREADBREAKS[i-1];
			}
			else
			{
				if(i==0)
					gds[i]._start = _PROCESSBREAKS[rank-1];
				else
					gds[i]._start = _THREADBREAKS[i-1];
			}
			gds[i]._end = _THREADBREAKS[i];
		}
		
		ExecutorService e = Executors.newFixedThreadPool(_NUMTHREADS);
		for(int i=0;i<_NUMTHREADS;i++)
		{
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
//		System.out.println("rank " + rank + " finished job, start = " + gds[0]._start + " , end = " + gds[0]._end);
		
		for(int i=0; i<_NUMTHREADS; i++) {
			changes += gds[i]._changes;
		}
		
		long g_changes[] = {0};
		long t_changes[] = {0};
		t_changes[0] = changes;
		// reduce changes
		
		MPI.COMM_WORLD.allReduce(t_changes, g_changes, 1, MPI.LONG, MPI.SUM);
		
		changes = g_changes[0];
	    
    	stTime = System.currentTimeMillis() - stTime;
    	if(rank==0)
    		System.out.println("All MPI job finished \ttime: " + stTime + "\tchanges: " + changes + "\tskipWords: " + _numSkipWords);
    	_numSkipWords = 0;
    	//System.out.println(Arrays.toString(_topicMarginal));
    	return changes;
    }
	
    
    
	//returns array of amt_i
	//such that if we divide leaf counts by amt_i, we get smoothing with marginal P(z)
    //(see paper)
	public static double [] getNormalizers(SparseBackoffTree [] shds, SparseBackoffTreeStructure struct, int partialUpdate, int from) throws MPIException{
		int numStates = struct.numLeaves();
		double [] count = new double[numStates];
		double [] smoothing = new double[numStates];
		
		int rank = MPI.COMM_WORLD.getRank();
		SparseBackoffTree shdAgg;
		
		if(from==0 && partialUpdate==1){
			// add only one time for duplicated word
			
			shdAgg = SparseBackoffTree.sum(shds, struct, _wordReduceFlag, rank);
			
		}
		else{
			shdAgg = SparseBackoffTree.sum(shds, struct);
		}
		
		/*
		double t_countsum = 0.0;
		double t_totalmass = 0.0;
		for(int j=0;j<shds.length;j++){
			shdAgg = shds[j];
			t_totalmass = shds[j]._totalMass;
			for(int i=0;i<numStates;i++){
				double [] smoothAndCount = shdAgg.getSmoothAndCount(i);
				t_countsum += smoothAndCount[1];
			}
		}
		
		System.out.println("local count sum = " + t_countsum);
		System.out.println("local mass = " + t_totalmass);
		*/
		
		
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
		
		
		if(partialUpdate==1 && from==0){	
			
			MPI.COMM_WORLD.allReduce(count, count.length, MPI.DOUBLE, MPI.SUM);
			MPI.COMM_WORLD.allReduce(smoothing, smoothing.length, MPI.DOUBLE, MPI.SUM);
			
			maxSmooth = 0.0;
			sumSmoother = 0.0;
			sumCount = 0.0;
			
			for(int i=0; i<numStates; i++) {
				sumSmoother += smoothing[i];
				sumCount += count[i];
				if(smoothing[i] > maxSmooth) {
					maxSmooth = smoothing[i];
				}
			}
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
	
	/**
	 * Updates the model given the topic assignments (_z) and divides by marginal for next sampling pass
	 */
	// Zheng added 'partialUpdate' parameter: 1, partial update wordToState    0, update all
	
    public void updateModel(TIntArrayList [] zs,int partialUpdate) throws MPIException{
    	int rank = MPI.COMM_WORLD.getRank();
    	int size = MPI.COMM_WORLD.getSize();
    	
    	if(rank==0){
    		System.out.println("arch: " + Arrays.toString(this._branchingFactors));
//    		System.out.println("second doc samples: " + zs[1].toString());
    	}
    	
		this._startState = getParamsFromZs(_forwardDelta, -2, zs, _c._docs, partialUpdate)[0];
		this._forward = getParamsFromZs(_forwardDelta, -1, zs, _c._docs, partialUpdate);
		this._backward = getParamsFromZs(_backwardDelta, 1, zs, _c._docs, partialUpdate);
		this._endState = getParamsFromZs(_backwardDelta, 2, zs, _c._docs, partialUpdate)[0];
		this._wordToState = getParamsFromZs(_wordDelta, 0, zs, _c._docs, partialUpdate);
				
		SparseBackoffTree [] comb = Arrays.copyOf(_backward, _backward.length + 1);
		comb[_backward.length] = _endState;
		
		//Zheng add:  partial get Normalizer or global get Normalizer : 1 partial update 0: update all
		
		this._backwardStateMarginal = getNormalizers(comb, _struct, partialUpdate, 1);
		
		this._wordStateMarginal = getNormalizers(_wordToState, _struct, partialUpdate, 0);
		
		for(int i=0; i<_backward.length; i++) {
			_backward[i].divideCountsBy(_backwardStateMarginal);
		}
		this._endState.divideCountsBy(_backwardStateMarginal);
		for(int i=0; i<_wordToState.length; i++) {
			if(_wordToState[i]!=null)
				_wordToState[i].divideCountsBy(_wordStateMarginal);	
		}
		
//		if(rank==0){
//			System.out.println("\tdiscounts forward: " + Arrays.toString(_forwardDelta) + 
//					"\tbackward " + Arrays.toString(_backwardDelta) +
//						"\tword " + Arrays.toString(_wordDelta));
//		}
    }
    
    /**
     * Adds the gradient of the log likelihood for the given observations
     * Assumes sum of deltas < 1
     * @param grad	the gradient is added here
     * @param hm	the observations
     * @param curDeltas	the current hyperparameters
     * @param incremental	if true, measure log likelihood when introducing observations one-by-one in random order;
     * otherwise, use leave-one-out cross validation
     * @return	log likelihood using current parameters
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
     * @param gradient	gradient to normalize
     * @param curDeltas	current hyperparameters to which gradient will be applied
     * @param stepSize	how far to move (in L1)
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
     * @param deltas	hyperparameters to start from
     * @param hm	observation counts
     * @param stepSize	how far to move in the step, in L1 norm
     * @param incremental	if true, measure log likelihood when introducing observations one-by-one in random order;
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
    public void optimizeParameters(TIntArrayList [] zs) throws MPIException{
//    	for(int i=0; i<_wordsDelta.length; i++) {
//    		_wordsDelta[i] *= 0.9;
//    		_docsDelta[i] *= 0.9;
//    	}
    	//TODO: multi-thread
//    	double STEPSIZE = 0.01; //start stepping this far in L1 norm
//    	double STEPDEC = 0.95; //decrease step size this much each step
//    	int NUMSTEPS = 20; //number of steps to take
    	double STEPSIZE = 0.02; //start stepping this far in L1 norm
    	double STEPDEC = 0.8; //decrease step size this much each step
    	int NUMSTEPS = 10; //number of steps to take

    	TIntDoubleHashMap [] hm = aggregateCounts(0, zs, _c._docs,1);
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
    	hm = aggregateCounts(-1, zs, _c._docs, 1);
    	TIntDoubleHashMap [] hm2 = Arrays.copyOf(hm, hm.length + 1);
    	hm2[hm.length] = aggregateCounts(-2, zs, _c._docs, 1)[0];
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
    	hm = aggregateCounts(1, zs, _c._docs, 1);
    	hm2 = Arrays.copyOf(hm, hm.length + 1);
    	hm2[hm.length] = aggregateCounts(2, zs, _c._docs, 1)[0];
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
     * @param expansionIndex	The desired index in _EXPANSIONSCHEDULE to expand to.
     * @return	whether the model was expanded
     */
    public boolean expandModel(int expansionIndex) throws MPIException{
    	if(expansionIndex < _EXPANSIONSCHEDULE.length) {
    		int curLeaves = _struct.numLeaves(); 
    		//get new structure
    		_branchingFactors = Arrays.copyOf(this._expansionBranchFactors, _EXPANSIONSCHEDULE[expansionIndex] + 1);
    		_struct = new SparseBackoffTreeStructure(_branchingFactors);
    		initDeltas();
    		int multiplier = _struct.numLeaves()/curLeaves;
    		int rank = MPI.COMM_WORLD.getRank();
    		if(rank==0)
    			System.out.println("Expanding to " + Arrays.toString(_branchingFactors));
    		for(int i=0; i<_c._z.length;i++) {
    			for(int j=0; j<_c._z[i].size(); j++) {
    				int z = multiplier * _c._z[i].get(j) + _r.nextInt(multiplier);
//    				int z = multiplier * _c._z[i].get(j) + 1;
    				_c._z[i].set(j,  z);
    			}
    		}
    		updateModel(_c._z, 1);
    		return true;
    	}
    	else
    		return false;
    }
    
    /**
     * Trains the model for a specified number of iterations.
     *
     * @param  iterations  the number of Gibbs passes to perform
     * @param  updateInterval	number of iterations between hyperparameter updates
     * @param  outFile   writes model to here before each expansion
     */
	public void trainModel(int iterations, int updateInterval, String outFile, int start_j) throws Exception {
		boolean done = false;
		int j=start_j;
		int maxVal = _struct.numLeaves();
		int rank = MPI.COMM_WORLD.getRank();
		int size = MPI.COMM_WORLD.getSize();
		int start = -1;
		int end = -1;
		
		FileWriter fw = null;
		BufferedWriter bw = null;
		
		FileWriter fw_1 = null;
		BufferedWriter bw_1 = null;
		
		long totalTime = 0;
		if(rank==0){

			fw = new FileWriter("./time.log");
			bw = new BufferedWriter(fw);
			bw.write("MPI_num = " + size);
			bw.newLine();
			bw.write("Thread per MPI = " + _NUMTHREADS);
			bw.newLine();
			bw.write("total thread = " + size*_NUMTHREADS);
			bw.newLine();
			bw.flush();
			
			fw_1 = new FileWriter("./changes.log");
			bw_1 = new BufferedWriter(fw_1);
		}
		
		_changesList = new LinkedList<Integer>();
		
		long totalChanges = 0;
		long initTimer = 0;
		long gibbsTimer = 0;
		long updateTimer = 0;
		long optTimer = 0;
		
		while(!done) {
			for(int i=1; i<=iterations*(int)(1.0/gibbsChunkSize); i++) {
				
				MPI.COMM_WORLD.barrier();
				long stTime = System.currentTimeMillis();
				long tTime = stTime;
				if(rank==0)
					System.out.println("before initZ: "+timeTag());
				
				initZ(_c._docs, maxVal, false, true, null); //init scratch array to zero
				
				MPI.COMM_WORLD.barrier();
				if(rank==0)
					System.out.println("after initZ & before gibbs pass"+timeTag());
				tTime = System.currentTimeMillis() - tTime;
				initTimer += tTime;
				
				tTime = System.currentTimeMillis();
				
		    	for(int x=2000000;x<2000005;x++){
		    		if(_c._z[x]!=null){
		    			System.out.println(x + " doc-samles: " + _c._z[x].toString() + " doc-content: " + _c._docs[x].toString());
		    		}
		    	}
		    	
				long changes = gibbsPass();
				
				MPI.COMM_WORLD.barrier();
				tTime = System.currentTimeMillis() - tTime;
				gibbsTimer += tTime;
				
				gibbsNextPointer += gibbsChunkSize;
				
				if(gibbsNextPointer>=1.0)
					gibbsNextPointer = 0.0;
				
				if(rank==0)
					System.out.println("after gibbs pass" + timeTag());
				
				totalChanges += changes;
				
				MPI.COMM_WORLD.barrier();
				tTime = System.currentTimeMillis();
				if((i % (updateInterval*(int)(1.0/gibbsChunkSize)))==0) {
					if(rank==0)
						System.out.println("before optimize" + timeTag());
					optimizeParameters(_c._scratchZ);
					if(rank==0)
						System.out.println("after optimize" + timeTag());
				}
				
				MPI.COMM_WORLD.barrier();
				tTime = System.currentTimeMillis() - tTime;
				optTimer += tTime;
				
				if(rank==0)
					System.out.println("before update Model" + timeTag());
				tTime = System.currentTimeMillis();
				
				updateModel(_c._scratchZ, 1);
				
				MPI.COMM_WORLD.barrier();
				tTime = System.currentTimeMillis() - tTime;
				updateTimer += tTime;
				
				_c._z = _c._scratchZ;
				MPI.COMM_WORLD.barrier();
				if(rank==0)
					System.out.println("passed the barrier after updateModel" + timeTag());
				stTime = System.currentTimeMillis() - stTime;
				totalTime += stTime;
				
				double memSize[] = getUsedMemSize();
				
				if(rank==0){
					
					bw_1.write("iteration" + i + "\n");
					bw_1.write("ave mem usage = " + memSize[0] + "GB" + "\n");
					bw_1.write("max mem usage = " + memSize[1] + "GB" + "\n");
					bw_1.write("changes= " + Long.toString(totalChanges) + "\n");
					bw_1.write("wordStateMarginal= " );
					for(int k=0;k<_wordStateMarginal.length;k++)
						bw_1.write(Double.toString(_wordStateMarginal[k])+", ");
					bw_1.newLine();
					bw_1.write("_forward[0]._totalMass= " + Double.toString(_forward[0]._totalMass) + "\n");
					bw_1.write("discounts forward: " + Arrays.toString(_forwardDelta) + 
							"\tbackward " + Arrays.toString(_backwardDelta) +
								"\tword " + Arrays.toString(_wordDelta)+"\n\n");
					bw_1.flush();
					totalChanges = 0;
					bw.write("iter " + i + " initZ_time= " + initTimer + "\n");
					bw.write("iter " + i + " gibbs_time= " + gibbsTimer + "\n");
					bw.write("iter " + i + " update_time= " + updateTimer + "\n");
					bw.write("iter " + i + " optimization_time= " + optTimer + "\n");
					bw.write("iter " + i + " total_time= " + stTime + "\n");
					bw.newLine();
					bw.flush();
					initTimer = 0;
					gibbsTimer = 0;
					updateTimer = 0;
					optTimer = 0;
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
		
		if(rank==0){
			bw.write("total time = " + totalTime);
			bw.newLine();
			bw.close();
			System.out.println("total training time = " + totalTime + timeTag());
			
			bw_1.close();
		}
	}
	
	private String timeTag(){
		String timeStamp = new SimpleDateFormat("HH:mm:ss").format(Calendar.getInstance().getTime());
		return " (time = " + timeStamp +")";
	}
	
	private boolean changesLevelOff(int changes){
		
		Integer changesArray[];
		double sum = 0.0;
		
		if(_changesList.size()>=listLength){
			_changesList.remove(0);
			_changesList.add(changes);
			changesArray = _changesList.toArray(new Integer[listLength]);
			
			for(int i=0;i<listLength-1;i++){
				if(changesArray[i]<=changesArray[i+1])
					sum += 1.0/(double)(listLength-1);
			}
			
			double diff;
			
			if(sum > listTrends){
				for(int i=0;i<listLength-1;i++){
					diff = (double)(Math.abs(changesArray[i] - changesArray[listLength-1]))/(double)(changesArray[listLength-1]);
					if(diff>listThreshold)
						return false;
				}
				return true;
			}
			else{
				return false;
			}
		}
		else
		{
			_changesList.add(changes);
			return false;
		}
	}
	
	public void writeModel(String outFile) throws Exception {
		SparseBackoffTree [] sbtW = this._wordToState;
		SparseBackoffTree [] sbtF = this._forward;
		SparseBackoffTree [] sbtB = this._backward;
		SparseBackoffTree sbtS = this._startState;
		SparseBackoffTree sbtE = this._endState;
		SparseBackoffTreeStructure struct = this._struct;
		int wordStoreFlag[] = this._wordStoreFlag;
		int wordReduceFlag[] = this._wordReduceFlag;
		
		_struct = null;
		_wordToState = null;
		_forward = null;
		_backward = null;
		_startState = null;
		_endState = null;
		_wordStoreFlag = null;
		_wordReduceFlag = null;
		
		//collect _c._z[]
		

		
		// Use MPI_Allgetherv instead, maybe
		int rank = MPI.COMM_WORLD.getRank();
		int size = MPI.COMM_WORLD.getSize();
		int start,end;
		
		int buff[] = null;
		int buff_array_size[] = null;
		
		
		// copy _c._z
		
		for(int k=1;k<size;k++)
		{
			int root = k;
			if(k==0)
				start = 0;
			else
				start = _PROCESSBREAKS[k-1];
			
			end = _PROCESSBREAKS[k];
			
			int buff_size = end-start;
			buff_array_size = new int[buff_size];
			
			int total_count = 0;
			
			if(rank==root){
				for(int l=start;l<end;l++){
					buff_array_size[l-start] = _c._z[l].size();
					total_count += _c._z[l].size();
				}
				MPI.COMM_WORLD.send(buff_array_size, buff_array_size.length, MPI.INT, 0, 0);
				buff = new int[total_count];
				
				int p=0;
				for(int l=start;l<end;l++){
					int ori_array[] = _c._z[l].toArray();
					for(int i=0;i<ori_array.length;i++)
						buff[p++] = ori_array[i];
				}
				MPI.COMM_WORLD.send(buff, buff.length, MPI.INT, 0, 1);
			}
			else if(rank==0){
				MPI.COMM_WORLD.recv(buff_array_size, buff_array_size.length, MPI.INT, root, 0);
				for(int i=0;i<buff_array_size.length;i++)
					total_count += buff_array_size[i];
				buff = new int[total_count];
				MPI.COMM_WORLD.recv(buff, buff.length, MPI.INT, root, 1);
				int p = 0;
				for(int i=0;i<buff_array_size.length;i++){
					int ori_array[] = new int[buff_array_size[i]];
					for(int j=0;j<buff_array_size[i];j++)
						ori_array[j] = buff[p++];
					TIntArrayList t = new TIntArrayList();
					t.reset();
					t.add(ori_array);
					_c._z[start+i] = t;
				}
			}
		}
		
		if(rank==0)
		{
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outFile+".part_1"));
			oos.writeObject(this);
			oos.close();
			
			for(int i=0;i<_c._z.length;i++){
				if(!within_doc_range(i))
					_c._z[i]=null;
			}
		}	
		
		// copy _c._doc
		
		for(int k=1;k<size;k++)
		{
			int root = k;
			if(k==0)
				start = 0;
			else
				start = _PROCESSBREAKS[k-1];
			
			end = _PROCESSBREAKS[k];
			
			int buff_size = end-start;
			buff_array_size = new int[buff_size];
			
			int total_count = 0;
			
			if(rank==root){
				for(int l=start;l<end;l++){
					buff_array_size[l-start] = _c._docs[l].size();
					total_count += _c._docs[l].size();
				}
				MPI.COMM_WORLD.send(buff_array_size, buff_array_size.length, MPI.INT, 0, 0);
				buff = new int[total_count];
				
				int p=0;
				for(int l=start;l<end;l++){
					int ori_array[] = _c._docs[l].toArray();
					for(int i=0;i<ori_array.length;i++)
						buff[p++] = ori_array[i];
				}
				MPI.COMM_WORLD.send(buff, buff.length, MPI.INT, 0, 1);
			}
			else if(rank==0){
				MPI.COMM_WORLD.recv(buff_array_size, buff_array_size.length, MPI.INT, root, 0);
				for(int i=0;i<buff_array_size.length;i++)
					total_count += buff_array_size[i];
				buff = new int[total_count];
				MPI.COMM_WORLD.recv(buff, buff.length, MPI.INT, root, 1);
				int p = 0;
				for(int i=0;i<buff_array_size.length;i++){
					int ori_array[] = new int[buff_array_size[i]];
					for(int j=0;j<buff_array_size[i];j++)
						ori_array[j] = buff[p++];
					TIntArrayList t = new TIntArrayList();
					t.reset();
					t.add(ori_array);
					_c._docs[start+i] = t;
				}
			}
		}

		if(rank==0)
		{
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outFile+".part_2"));
			oos.writeObject(this);
			oos.close();
			for(int i=0;i<_c._docs.length;i++){
				if(!within_doc_range(i))
					_c._docs[i]=null;
			}
		}
		_wordToState = sbtW;
		_forward = sbtF;
		_backward = sbtB;
		_struct = struct;
		_startState = sbtS;
		_endState = sbtE;
		_wordStoreFlag = wordStoreFlag;
		_wordReduceFlag = wordReduceFlag;
		
		/*
		int blank_count = 0;
		for(int i=0;i<_wordToState.length;i++)
		{
			if(_wordToState[i]==null)
				blank_count++;
		}
		
		System.out.println("rank :" + rank + "  "+"word to state size = " + _wordToState.length + "blank count" + blank_count + "num of leaf = " + _wordToState[0]._struct.numLeaves());
		System.out.println("rank :" + rank + "  "+"forward size = " + _forward.length+ "num of leaf = " + _forward[0]._struct.numLeaves());
		System.out.println("rank :" + rank + "  "+"backward size = " + _backward.length + "num of leaf = " + _backward[0]._struct.numLeaves());
		*/

	}
	
	public static SBTSequenceModel_mpi readModel(String inFile) throws Exception {
		ObjectInputStream ois = new ObjectInputStream(new FileInputStream(inFile+".part_1"));
		ObjectInputStream ois_2 = new ObjectInputStream(new FileInputStream(inFile+".part_2"));
		
		SBTSequenceModel_mpi out = (SBTSequenceModel_mpi) ois.readObject();
		SBTSequenceModel_mpi out_2 = (SBTSequenceModel_mpi) ois_2.readObject();
		ois.close();
		ois_2.close();
		
		out._c._docs = out_2._c._docs;
		
		out_2 = null;
		
		out._struct = new SparseBackoffTreeStructure(out._branchingFactors);
		System.out.println("before update");
		out.updateModel(out._c._z, 0);
		System.out.println("after update");
		return out;
	}
	
	
	
	//returns model x word array of probabilities 
	public static double [][] getEnsembleProbsForDoc(int doc, double [][] wordTopicMarginal, SBTSequenceModel_mpi [] ms,
			double [][][] tms) {
		
		double [][] ps = new double[ms.length][ms[0]._c._docs[doc].size()]; //dims [model][token]
		for(int m=0; m<ms.length; m++) {
			int numStates = wordTopicMarginal[m].length;
			SBTSequenceModel_mpi mod = ms[m];
			TIntArrayList words = mod._c._docs[doc];
			double [][] dist = new double[words.size()+1][numStates];
			double [][] tm = tms[m];
			
			dist[0][0] = 1.0;
			for(int i=0; i<dist.length - 1; i++) {
				int w = words.get(i);
				//step from i to i+1:
				double [] pWordGivenState = new double[numStates];
				double [] pStateGivenPrev = new double[numStates];
//				for(int j=0; j<dist[i].length;j++) {
//					sbtForward.addWeighted(this._forward[j], dist[i][j]);
//				}
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
	public static double [] testEnsembleOnDocExact(int doc, double [][] wordTopicMarginal, SBTSequenceModel_mpi [] ms,
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
	
	public static double [] testEnsembleOnDocFullPpl(int doc, double [][] wordTopicMarginal, SBTSequenceModel_mpi [] ms) {
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
//			for(int j=0; j<dist[i].length;j++) {
//				sbtForward.addWeighted(this._forward[j], dist[i][j]);
//			}
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
				if(_wordToState[w]==null){
					System.out.println("w = " + w + "  is empty -> replaced with unk wordID: 832");
					_wordToState[w] = _wordToState[832];
				}
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
	//Zheng added: incFile paramter: null if init, not null if restart from previous Z
	public static void train(String inputFile, String outputFile, String configFile, String probsFile, String incFile) throws Exception {
		int rank = MPI.COMM_WORLD.getRank();
		SBTSequenceModel_mpi sbtsm = new SBTSequenceModel_mpi(configFile);
		if(rank==0)
			System.out.println("----construction function----");
		sbtsm.initializeForCorpus(inputFile, sbtsm._struct.numLeaves(),incFile);
		
		if(probsFile != null) {
			sbtsm.baseProbs = sbtsm.readBaseProbs(probsFile);
		}
		if(incFile==null)
			sbtsm.trainModel(sbtsm._NUMITERATIONS, sbtsm._OPTIMIZEINTERVAL, outputFile,0);
		else
		{
			String[] aa = incFile.split("\\.");
			int start_j = Integer.parseInt(aa[aa.length-1]);
			System.out.println("before train mode : " + start_j);
			sbtsm.trainModel(sbtsm._NUMITERATIONS, sbtsm._OPTIMIZEINTERVAL, outputFile,start_j+1);
		}
		
		sbtsm.writeModel(outputFile);
	}
	
	public static void train(String inputFile, String outputFile, String configFile) throws Exception {
		train(inputFile, outputFile, configFile, null,null);
	}
	
	public long debugfun(int numthreads) throws Exception{
		
		int glen = 100000;
		
		int rank = MPI.COMM_WORLD.getRank();
		int size = MPI.COMM_WORLD.getSize();
		
		if(glen%size!=0)
			System.out.println("glen not dividable !!!!");
		//num of summation for each MPI process
		int llen = glen/size;
		
		if(llen%numthreads!=0)
			System.out.println("llen not dividable !!!");
		//num of summation for each threads
		int slen = llen/numthreads;
		
		int lstart = rank*llen;
		DebugDoer [] dds = new DebugDoer[numthreads];
		
		for(int i=0;i<numthreads;i++){
			dds[i] = new DebugDoer();
			dds[i]._start = lstart + slen*i;
			dds[i]._end = lstart + slen*(i+1)-1;
		}
		ExecutorService e = Executors.newFixedThreadPool(numthreads);
		MPI.COMM_WORLD.barrier();
		long stTime = System.currentTimeMillis();
		
		for(int i=0;i<numthreads;i++){
			e.execute(dds[i]);
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
		
		MPI.COMM_WORLD.barrier();
		stTime = System.currentTimeMillis() - stTime;
    	if(rank==0)
    		System.out.println("\ttime: " + stTime);
    	
		return 0;
	}
	
	public static void debug(int numthreads) throws Exception{
		SBTSequenceModel_mpi sbtsm;
		
		String configFile = "./SBT.config";
		
		sbtsm = new SBTSequenceModel_mpi(configFile);
		
		sbtsm.debugfun(numthreads);
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
				normalizer += ws[i]; //not actually necessary if gradient sums to zero
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
	
	public static double [] trainEnsemble(SBTSequenceModel_mpi [] sbtsm, double [][] wordStateNormalizer, double [][][] tm, String validationFile, int numDocs) throws Exception {
	  System.out.println("train ens num docs: " + numDocs);
		for(int i=0; i<sbtsm.length; i++) {
			sbtsm[i]._c.reInitForCorpus(validationFile, numDocs);
		}
		double LL = 0.0;
		double numWords = 0.0;
		
		
		TestGetProbsDoer [] tds = new TestGetProbsDoer[numDocs];
		ExecutorService e = Executors.newFixedThreadPool(1);
		double [] eqW = new double[sbtsm.length];
		Arrays.fill(eqW,  1.0/(double)eqW.length);
		//ExecutorService e = Executors.newFixedThreadPool(1);
		// mpi-incomplete here zheng
		for(int i=0; i<tds.length;i++) {
			if(sbtsm[0]._MAXTESTDOCSIZE < 0 || sbtsm[0]._c._docs[i].size() > sbtsm[0]._MAXTESTDOCSIZE) {
				System.out.println("skipping " + i + " size " + sbtsm[0]._c._docs[i].size());
				continue;
			}
			tds[i] = new TestGetProbsDoer(i, wordStateNormalizer, sbtsm, tm, eqW);
		}
		
		int rank = MPI.COMM_WORLD.getRank();
		int size = MPI.COMM_WORLD.getSize();
		
		for(int i=rank; i<tds.length;i+=size) {
			if(sbtsm[0]._MAXTESTDOCSIZE < 0 || sbtsm[0]._c._docs[i].size() > sbtsm[0]._MAXTESTDOCSIZE) {
				continue;
			}
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
		for(int i=rank; i<numDocs;i+=size) {
			ps[i] = tds[i]._probs;
		}
		
		int root;
		for(int i=0;i<numDocs;i++)
		{
			for(int j=0;j<ps[i].length;j++)
			{
				root = i%size;
				MPI.COMM_WORLD.bcast(ps[i][j], ps[i][j].length, MPI.DOUBLE, root);
			}
		}
		
		
		//multi thread 
		/*
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
		*/
		
		double [] outW = optimizeEnsembleWeights(ps);
		return outW;		
	}
	
	//computes log likelihood of model using left-to-right method
	//returns {ppl, number of tested words}
	//numDocs is number in test file
	//maxDocs is number actually tested (starting from the beginning of the file)
	public static double [] testEnsemble(File [] modelFile, String testFile, int numDocs, int maxDocs, String configFile,
			String validationFile, int validationDocs, String probsOutfile) throws Exception {
		SBTSequenceModel_mpi [] sbtsm = new SBTSequenceModel_mpi[modelFile.length];
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
		
		
		ExecutorService e = Executors.newFixedThreadPool(1);
		//ExecutorService e = Executors.newFixedThreadPool(1);
		int rank = MPI.COMM_WORLD.getRank();
		int size = MPI.COMM_WORLD.getSize();
		
		for(int i=rank; i<tds.length;i+=size) {
			if(sbtsm[0]._MAXTESTDOCSIZE < 0 || sbtsm[0]._c._docs[i].size() > sbtsm[0]._MAXTESTDOCSIZE) {
				System.out.println("skipping " + i + " size " + sbtsm[0]._c._docs[i].size());
				continue;
			}
			tds[i] = new TestEnsembleExactDoer(i, wordStateNormalizer, sbtsm, tm, ws);
		}
		for(int i=rank; i<tds.length;i+=size) {
			if(sbtsm[0]._MAXTESTDOCSIZE < 0 || sbtsm[0]._c._docs[i].size() > sbtsm[0]._MAXTESTDOCSIZE) {
				continue;
			}
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
		
		double t_LL[] = {0};
		double g_LL[] = {0};
		double t_numWords[] = {0};
		double g_numWords[] = {0};

		for(int i=rank; i<tds.length; i+=size) {
			if(tds[i]==null)
				continue;
			t_LL[0] += tds[i]._res[0];
			t_numWords[0] += tds[i]._res[1];
			psOut[i] = tds[i]._ps;
		}
		
		MPI.COMM_WORLD.allReduce(t_LL, g_LL, 1, MPI.DOUBLE, MPI.SUM);
		MPI.COMM_WORLD.allReduce(t_numWords, g_numWords, 1, MPI.DOUBLE, MPI.SUM);
		
		LL = g_LL[0];
		numWords = g_numWords[0];
		int root;
		
		for(int i=0;i<tds.length;i++)
		{
			root = i%size;
			MPI.COMM_WORLD.bcast(psOut[i], psOut[i].length, MPI.DOUBLE, root);
		}
		
		System.out.println("Test LL: " + LL);
		System.out.println("Words: " + numWords);
		System.out.println("ppl: " + Math.exp(-LL/numWords));
		if(probsOutfile != null && rank==0) {
			writeBaseProbs(probsOutfile, psOut);
		}
		
		
		//multi thread
		/*
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
		
		*/
		return new double [] {LL, numWords};
	}
	
	//computes log likelihood of model using left-to-right method
	//returns {ppl, number of tested words}
	//numDocs is number in test file
	//maxDocs is number actually tested (starting from the beginning of the file)
	public static double [] testModel(String modelFile, String testFile, int numDocs, int maxDocs, String configFile) throws Exception {
		boolean CHECKWORDDISTRIBUTION = false;
		
		SBTSequenceModel_mpi sbtsm = readModel(modelFile);
		System.out.println("forward");
		for(int i=0;i<20;i++){
				System.out.println(Arrays.toString(sbtsm._forward[i].toDoubleArray()));
		}
		System.out.println("backward");
		for(int i=0;i<20;i++){
				System.out.println(Arrays.toString(sbtsm._backward[i].toDoubleArray()));
		}
		
		return null;
		
		/*
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
		sbtsm.setThreadBreaks(1);
		double LL = 0.0;
		double numWords = 0.0;
		//TestDoer [] tds = new TestDoer[maxDocs];
		TestExactDoer [] tds = new TestExactDoer[maxDocs];
		
		//mpi
		
		ExecutorService e = Executors.newFixedThreadPool(1);
		//ExecutorService e = Executors.newFixedThreadPool(1);
		
		System.out.println("tds . length = " + tds.length + " total doc =  " + sbtsm._c._docs.length);
		
		for(int i=0; i<tds.length;i++) {
			if(sbtsm._MAXTESTDOCSIZE < 0 || sbtsm._c._docs[i].size() > sbtsm._MAXTESTDOCSIZE) {
				System.out.println("skipping " + i + " size " + sbtsm._c._docs[i].size());
				continue;
			}
			tds[i] = sbtsm.new TestExactDoer(i, wordStateNormalizer, tm);
		}
		int rank = MPI.COMM_WORLD.getRank();
		int size = MPI.COMM_WORLD.getSize();
		
		for(int i=rank; i<tds.length;i+=size) {
			if(sbtsm._MAXTESTDOCSIZE < 0 || sbtsm._c._docs[i].size() > sbtsm._MAXTESTDOCSIZE) {
				continue;
			}
			
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
		
		double t_LL[] = {0};
		double g_LL[] = {0};
		double t_numWords[] = {0};
		double g_numWords[] = {0};
		
		for(int i=rank; i<tds.length; i+=size) {
			if(tds[i]==null)
				continue;
			t_LL[0] += tds[i]._res[0];
			t_numWords[0] += tds[i]._res[1];
		}
		
		MPI.COMM_WORLD.allReduce(t_LL, g_LL, 1, MPI.DOUBLE, MPI.SUM);
		MPI.COMM_WORLD.allReduce(t_numWords, g_numWords, 1, MPI.DOUBLE, MPI.SUM);

		LL = g_LL[0];
		numWords = g_numWords[0];
		
		
		
		
		return new double [] {LL, numWords};
		*/
	}
	
	// multi thread
		/*
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
		*/
	
	public static void test(String modelFile, String inputFile, int numDocsInFile, int numDocsToTest, String configFile) throws Exception {
		double [] ll = testModel(modelFile, inputFile, numDocsInFile, numDocsToTest, configFile);
		System.out.println("Test LL: " + Arrays.toString(ll));
		System.out.println("ppl: " + Math.exp(-ll[0]/ll[1]));
	}
	
	/**
	 * Outputs the sparse word-to-topic vector proportional to P(t | w)/P(t) for words and topics with positive counts
	 * @param modelFile	Contains the model
	 * @param outputFile	Where to write the output
	 * @throws Exception
	 */
	public static void outputWordToTopicFile(String modelFile, String outputFile) throws Exception {
		BufferedWriter bwOut = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), "UTF8"));
		SBTSequenceModel_mpi sbtsm = readModel(modelFile);
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
	
	public void testAggregateCounts() throws MPIException{
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
		TIntDoubleHashMap [] hm = aggregateCounts(-1, oneZs, oneWs, 1);
		System.out.println("hm is " + hm.length);
		hm = aggregateCounts(1, oneZs, oneWs, 1);
		System.out.println("hm is " + hm.length);
		hm = aggregateCounts(0, oneZs, oneWs, 1);
		System.out.println("hm is " + hm.length);
	}
	
	public static void main(String[] args) throws Exception {
		MPI.Init(args);
		if(args.length > 0 && args[0].equalsIgnoreCase("train")) {
			if(args.length != 4) {
				int rank = MPI.COMM_WORLD.getRank();
				if(rank==0)
					System.err.println("Usage: train <input_file> <model_output_file> <configuration_file>");
				MPI.Finalize();
				return;
			}
			
			InetAddress ip;
	        String hostname;

			ip = InetAddress.getLocalHost();
            hostname = ip.getHostName();
			
			train(args[1], args[2], args[3]);
		}
		else if(args.length > 0 && args[0].equalsIgnoreCase("trainInc")) { //restart from previous module
			if(args.length != 5) {
				int rank = MPI.COMM_WORLD.getRank();
				if(rank==0)
					System.err.println("Usage: trainInc <input_file> <model_output_file> <configuration_file> <restart_model_file>");
				MPI.Finalize();
				return;
			}
			train(args[1], args[2], args[3], null,args[4]);
		}
		else if(args.length > 0 && args[0].equalsIgnoreCase("trainAug")) { //augment existing probs
			if(args.length != 5) {
				System.err.println("Usage: train <input_file> <model_output_file> <configuration_file> <probs_file>");
				MPI.Finalize();
				return;
			}
			train(args[1], args[2], args[3], args[4],null);
		}
		else if(args.length > 0 && args[0].equalsIgnoreCase("test")) {
			int rank = MPI.COMM_WORLD.getRank();
			int size = MPI.COMM_WORLD.getSize();
			
			if(args.length != 6) {
				if(rank==0)
					System.err.println("Usage: test <model_file> <test_file> <configuration_file> <num_docs_in_test_file> <num_docs_to_test>");
				MPI.Finalize();
				return;
			}
			
			if(size>1){
				if(rank==0)
					System.err.println("test phase only supports 1 MPI process currently");
				MPI.Finalize();
				return;
			}
			
			test(args[1], args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3]);
		}
		else if(args.length > 0 && args[0].equalsIgnoreCase("wordreps")) {
			if(args.length != 3) {
				System.err.println("Usage: wordreps <model_file> <output_file>");
				MPI.Finalize();
				return;
			}
			outputWordToTopicFile(args[1], args[2]);
		}
		else if(args.length > 0 && args[0].equalsIgnoreCase("testEns")) {
			if(args.length != 6 && args.length != 8) {
				System.err.println("Usage: test <model_dir> <test_file> <configuration_file> <num_docs_in_test_file> <num_docs_to_test>"
						+ " [<validation_file> <num_validation_docs>]");
				MPI.Finalize();
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
				MPI.Finalize();
				return;
			}
			File f = new File(args[1]);
			if(args.length==7)
				testEnsemble(f.listFiles(), args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3], null, 0, args[6]);
			else
				testEnsemble(f.listFiles(), args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3], args[7], 
						Integer.parseInt(args[8]), args[6]);			
		}
		else if(args.length > 0 && args[0].equalsIgnoreCase("Debug")){
			if(args.length==2)
			debug(Integer.parseInt(args[1]));
		}
		else {
			System.err.println("Usage: <train|trainAug|test|wordreps|testEns|outputEnsProbs> <args...>");
			MPI.Finalize();
			return;
		}
		
		MPI.Finalize();
	}
}
