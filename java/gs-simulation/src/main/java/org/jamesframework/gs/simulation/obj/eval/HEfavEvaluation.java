

package org.jamesframework.gs.simulation.obj.eval;


public class HEfavEvaluation extends FrequencyBasedEvaluation {
    
    public HEfavEvaluation(int n, double[] alleleFreqs, int[] favAlleles) {
        super(n, alleleFreqs, favAlleles);
    }
    
    public HEfavEvaluation(HEfavEvaluation toCopy){
        super(toCopy);
    }
    
    // penalize fixation of unfavourable alleles only
    @Override
    public double getValue() {
        double val = 0.0;
        int numMarkers = getNumMarkers();
        for(int m = 0; m < numMarkers; m++){
            double p = getFavourableAlleleFreq(m);
            val += p * (2.0 - p);
        }
        val = val / numMarkers;
        return val;
    }
    
}
