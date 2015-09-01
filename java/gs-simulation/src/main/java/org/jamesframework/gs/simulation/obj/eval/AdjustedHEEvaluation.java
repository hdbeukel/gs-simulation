

package org.jamesframework.gs.simulation.obj.eval;


public class AdjustedHEEvaluation extends FavourableAlleleBasedEvaluation {
    
    public AdjustedHEEvaluation(int n, double[] alleleFreqs, int[] favAlleles) {
        super(n, alleleFreqs, favAlleles);
    }
    
    public AdjustedHEEvaluation(AdjustedHEEvaluation toCopy){
        super(toCopy);
    }
    
    // penalize fixation of unfavourable alleles only
    @Override
    public double getValue() {
        double[] alleleFreqs = getAlleleFreqs();
        double val = 0.0;
        int numMarkers = alleleFreqs.length;
        for(int m=0; m<numMarkers; m++){
            double f = getUnfavourableAlleleFreq(m);
            val += 1.0 - f*f;
        }
        val = 1.0/numMarkers * val;
        return val;
    }
    
}
