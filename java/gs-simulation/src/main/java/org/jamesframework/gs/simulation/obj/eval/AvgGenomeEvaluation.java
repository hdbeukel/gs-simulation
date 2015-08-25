

package org.jamesframework.gs.simulation.obj.eval;

import java.util.Arrays;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;

public abstract class AvgGenomeEvaluation implements Evaluation {

    // subset size
    private final int n;
    // averaged allele frequencies (ratio of 1-alleles)
    private final double[] alleleFreqs;

    public AvgGenomeEvaluation(int n, double[] alleleFreqs) {
        this.n = n;
        this.alleleFreqs = alleleFreqs;
    }
    
    public AvgGenomeEvaluation(AvgGenomeEvaluation toCopy){
        this.n = toCopy.n;
        this.alleleFreqs = Arrays.copyOf(toCopy.alleleFreqs, toCopy.alleleFreqs.length);
    }
    
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
