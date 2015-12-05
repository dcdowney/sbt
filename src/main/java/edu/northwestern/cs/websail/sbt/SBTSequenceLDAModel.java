package edu.northwestern.cs.websail.sbt;

import java.util.Arrays;

import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntDoubleHashMap;

/**
 * Adds a single additional SBT at the document/sentence level, to SBTSequenceModel.
 * @author dcdowney
 *
 */
public class SBTSequenceLDAModel extends SBTSequenceModel {

	private static final long serialVersionUID = 1L;
	
	double [] _docBagDelta;
	
	double [] _forwardStateMarginal; //we will divide by this now
	
	
	SBTSequenceLDAModel(String configFile) throws Exception {
		super(configFile);
	}
	
	@Override
	public void initDeltas() {
		super.initDeltas();
		_docBagDelta =new double[_branchingFactors.length];
		Arrays.fill(_docBagDelta, this._forwardDelta[0]); 
	}
	
	public int sampleZ(int doc, int pos, boolean testing, SparseBackoffTree sbtDoc) {
		int w = _c._docs[doc].get(pos);
		int curZ = _c._z[doc].get(pos);
		SparseBackoffTree sbtForward = null;
		if(pos==0) {
			sbtForward = this._forward[0];
		}
		else {
			sbtForward = this._forward[_c._z[doc].get(pos-1)];
		}
		SparseBackoffTree sbtBackward = null;
		boolean forward = (testing && (pos==_c._docs[doc].size() - 1 || _c._z[doc].get(pos+1)==-1)); 
		if(!forward) { //don't use backward distribution if doing forward prediction
			if(pos==_c._docs[doc].size() - 1) {  
				sbtBackward = this._backward[0];
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
		if(sbtWord._totalMass == 0.0) {
			unseen = true;
		}
		double bMarginalInv = 1.0f;
		double wMarginalInv = 1.0f;
		double fMarginalInv = 1.0f;
		if(curZ >= 0) {
			bMarginalInv = 1.0 / this._backwardStateMarginal[curZ];
			wMarginalInv = 1.0 / this._wordStateMarginal[curZ];
			fMarginalInv = 1.0 / this._forwardStateMarginal[curZ];
		}
		if(testing) { //subtract from only docs at test time:
			bMarginalInv = 0.0;
			wMarginalInv = 0.0;
			fMarginalInv = 0.0;
		}
		
		SparseBackoffTree [] sbts;
		double [] subs;
		if(!forward) {
			if(unseen) {
				sbts = new SparseBackoffTree [] {sbtForward, sbtBackward, sbtDoc};
				subs = new double [] {fMarginalInv, wMarginalInv, 1.0};
			}
			else {
				sbts = new SparseBackoffTree [] {sbtForward, sbtWord, sbtBackward, sbtDoc};
				subs = new double [] {fMarginalInv, wMarginalInv, bMarginalInv, 1.0};
			}
		}
		else {
			if(unseen) {
				sbts = new SparseBackoffTree [] {sbtForward, sbtDoc};
				subs = new double [] {fMarginalInv, 1.0};
			}
			else {
				sbts = new SparseBackoffTree [] {sbtForward, sbtWord, sbtDoc};
				subs = new double [] {fMarginalInv, wMarginalInv, 1.0};
			}
		}

		SparseBackoffTreeIntersection sbti = new SparseBackoffTreeIntersection(sbts, curZ, subs, true);
		int sample = sbti.sample(_r);
		return sample;
	}

	//scans doc, resampling each variable.  Not used at test time.
	//returns number of changes
	@Override
	protected int sampleDoc(int docId) {
		int changes = 0;
		TIntArrayList zs = _c._z[docId];
		TIntArrayList doc = _c._docs[docId];
		TIntDoubleHashMap docCts = aggregateCounts(zs);
		SparseBackoffTree sbtDoc = new SparseBackoffTree(_struct);
		//System.out.println("adding " + docCts.toString());
		sbtDoc.addAllMass(docCts, this._docBagDelta);
		for(int i=0; i<doc.size(); i++) {
			int newZ = sampleZ(docId, i, false, sbtDoc);
			if(newZ != zs.get(i))
				changes++;
			zs.set(i, newZ);
		}
		return changes;
	}
	
	@Override
	/**
	 * Follows superclass, except if from=2, returns doc counts
	 */
	public TIntDoubleHashMap [] aggregateCounts(int from, TIntArrayList [] zs, TIntArrayList [] ws) {
		if(from >=-1 && from <=1)
			return super.aggregateCounts(from, zs, ws);
		else if(from==2) {
			TIntDoubleHashMap [] hm = new TIntDoubleHashMap[zs.length];
			for(int i=0; i<zs.length; i++) {
				hm[i] = aggregateCounts(zs[i]);
			}
			return hm;
		}
		else {
			System.err.println("unsupported from type.");
			return null;
		}
	}
	
	@Override
    /**
     * Optimizes parameters using gradient ascent in log likelihood
     */
    public void optimizeParameters(TIntArrayList [] zs) {
		
    	//TODO: multi-thread
    	double STEPSIZE = 0.02; //start stepping this far in L1 norm
    	double STEPDEC = 0.8; //decrease step size this much each step
    	int NUMSTEPS = 10; //number of steps to take

    	TIntDoubleHashMap [] hm = aggregateCounts(2, zs, _c._docs);
    	System.out.println("docs:");
    	double step = STEPSIZE;
    	for(int i=0; i<NUMSTEPS; i++) {
        	gradientStep(_docBagDelta, hm, step, false);
        	step *= STEPDEC;
    	}
    	long totalparams = 0L;
    	for(TIntDoubleHashMap i : hm) {
    		if(i != null)
    			totalparams += i.size();
    	}
    	super.optimizeParameters(zs);
    	System.out.println("full doc params: " + Arrays.toString(_docBagDelta));
    	System.out.println("total nonzero doc-state params: " + totalparams);
    }

	/**
	 * Updates the model given the topic assignments (_z) and divides by marginal for next sampling pass
	 */
    public void updateModel(TIntArrayList [] zs) {
    	super.updateModel(zs);
    	
		this._forwardStateMarginal = getNormalizers(_forward, _struct);
		for(int i=0; i<_wordToState.length; i++) {
			_forward[i].divideCountsBy(_forwardStateMarginal);	
		}
		System.out.println("\tdiscounts forward: " + Arrays.toString(_forwardDelta) + 
				"\tbackward " + Arrays.toString(_backwardDelta) +
						"\tword " + Arrays.toString(_wordDelta));
    }

	public static void train(String inputFile, String outputFile, String configFile) throws Exception {
		SBTSequenceLDAModel sbtsm = new SBTSequenceLDAModel(configFile);
		sbtsm.initializeForCorpus(inputFile, sbtsm._struct.numLeaves());
		sbtsm.trainModel(sbtsm._NUMITERATIONS, sbtsm._OPTIMIZEINTERVAL, outputFile);
		sbtsm.writeModel(outputFile);
	}
	
	public static void test(String modelFile, String inputFile, int numDocsInFile, int numDocsToTest, String configFile) throws Exception {
		double [] ll = testModel(modelFile, inputFile, numDocsInFile, numDocsToTest, configFile);
		System.out.println("Test LL: " + Arrays.toString(ll));
		System.out.println("ppl: " + Math.exp(-ll[0]/ll[1]));
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
				
				if(_c._pWord[w]>0.0) {
					//init topicGivenDoc with all samples:
					SparseBackoffTree topicGivenDoc = new SparseBackoffTree(_struct);					
					for(int k =0; k<_c._z[doc].size(); k++) {
						int r = _c._z[doc].get(k);
						if(r != -1)
							topicGivenDoc.addAndSmoothIfZero(r, 1.0, this._docBagDelta);
					}
					//resample everything before:
					for(int k=0; k<i; k++) {
						int r = sampleZ(doc, k, true, topicGivenDoc);
						_c._z[doc].set(k, r);
					}

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

	
	//TODO: make a runner class so I don't have this copy/paste of main (and train, test) across model classes
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
			}
			test(args[1], args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[3]);
		}
		else {
			System.err.println("Usage: <train|test> <args...>");
		}
	}
	
}
