

package org.jamesframework.gs.simulation.obj.eval;

import java.util.Arrays;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;

/**
 * Evaluation based on allele frequencies in the selection.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public abstract class FrequencyBasedEvaluation implements Evaluation {

    // subset size
    private final int subsetSize;
    // allele frequencies (ratio of 1-alleles)
    private final double[] alleleFreqs;
    // indicates which allele is favourable
    private final int[] favAlleles;

    public FrequencyBasedEvaluation(int subsetSize, double[] alleleFreqs, int[] favAlleles) {
        this.subsetSize = subsetSize;
        this.alleleFreqs = alleleFreqs;
        this.favAlleles = favAlleles;
    }
    
    public FrequencyBasedEvaluation(FrequencyBasedEvaluation toCopy){
        this.subsetSize = toCopy.subsetSize;
        this.alleleFreqs = Arrays.copyOf(toCopy.alleleFreqs, toCopy.alleleFreqs.length);
        this.favAlleles = toCopy.favAlleles; // no deep copy (doesn't change across selections!)
    }
    
    // update allel frequencies when swapping an item in the selection
    public void swap(int[] delGenome, int[] addGenome){
        int numMarkers = alleleFreqs.length;
        for(int m=0; m<numMarkers; m++){
            alleleFreqs[m] *= subsetSize;
            alleleFreqs[m] -= delGenome[m];
            alleleFreqs[m] += addGenome[m];
            alleleFreqs[m] /= subsetSize;
        }
    }

    public int getSubsetSize() {
        return subsetSize;
    }
    
    public int getNumMarkers(){
        return alleleFreqs.length;
    }

    public double getAlleleFreq(int m) {
        return alleleFreqs[m];
    }
    
    public double getMinorAllelFreq(int m){
        double p = getAlleleFreq(m);
        double maf = Math.min(p, 1.0 - p);
        return maf;
    }
    
    public double getFavourableAlleleFreq(int m){
        // get frequency of 1-allele at marker m
        double freq = getAlleleFreq(m);
        // flip if 0-allele is favourable
        if(favAlleles[m] == 0){
            freq = 1.0 - freq;
        }
        return freq;
    }
    
    public double getUnfavourableAlleleFreq(int m){
        return 1.0 - getFavourableAlleleFreq(m);
    }
    
}
