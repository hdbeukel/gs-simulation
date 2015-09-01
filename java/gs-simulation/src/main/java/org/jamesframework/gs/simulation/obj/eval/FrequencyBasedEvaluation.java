

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
    private final int n;
    // allele frequencies (ratio of 1-alleles)
    private final double[] alleleFreqs;

    public FrequencyBasedEvaluation(int n, double[] alleleFreqs) {
        this.n = n;
        this.alleleFreqs = alleleFreqs;
    }
    
    public FrequencyBasedEvaluation(FrequencyBasedEvaluation toCopy){
        this.n = toCopy.n;
        this.alleleFreqs = Arrays.copyOf(toCopy.alleleFreqs, toCopy.alleleFreqs.length);
    }
    
    // update allel frequencies when swapping an item in the selection
    public void swap(int[] delGenome, int[] addGenome){
        int numMarkers = alleleFreqs.length;
        for(int m=0; m<numMarkers; m++){
            alleleFreqs[m] *= n;
            alleleFreqs[m] -= delGenome[m];
            alleleFreqs[m] += addGenome[m];
            alleleFreqs[m] /= n;
        }
    }

    public int getSubsetSize() {
        return n;
    }

    public double[] getAlleleFreqs() {
        return alleleFreqs;
    }
    
}
