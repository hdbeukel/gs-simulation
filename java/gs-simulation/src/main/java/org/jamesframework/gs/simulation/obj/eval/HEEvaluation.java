

package org.jamesframework.gs.simulation.obj.eval;

import java.util.Arrays;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;


public class HEEvaluation implements Evaluation {

    // subset size
    private final int n;
    // averaged allele frequencies
    private final double[] alleleFreqs;

    public HEEvaluation(int n, double[] alleleFreqs) {
        this.n = n;
        this.alleleFreqs = alleleFreqs;
    }
    
    public HEEvaluation(HEEvaluation toCopy){
        this.n = toCopy.n;
        this.alleleFreqs = Arrays.copyOf(toCopy.alleleFreqs, n);
    }
    
    public void swap(int[] delGenome, int[] addGenome){
        for(int m=0; m<n; m++){
            alleleFreqs[m] -= ((double) delGenome[m])/n;
            alleleFreqs[m] += ((double) addGenome[m])/n;
        }
    }

    @Override
    public double getValue() {
        double he = 0.0;
        for(int m=0; m<n; m++){
            he += alleleFreqs[m] * (1.0 - alleleFreqs[m]);
        }
        he = 2.0/n * he;
        return he;
    }
    
}
