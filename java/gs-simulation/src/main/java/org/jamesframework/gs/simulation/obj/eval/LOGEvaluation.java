

package org.jamesframework.gs.simulation.obj.eval;


public class LOGEvaluation extends FavourableAlleleBasedEvaluation {
    
    private final double truncate;
    
    public LOGEvaluation(int subsetSize, double[] alleleFreqs, int[] favAlleles) {
        super(subsetSize, alleleFreqs, favAlleles);
        // compute truncation level: set difference
        // between 0 and 1 occurrence  = 1 > diff between
        // 1 and 2 occurrences = ln(2) = 0.693
        double valSingleOcc = -Math.log(subsetSize);
        truncate = valSingleOcc - 1.0;
    }
    
    public LOGEvaluation(LOGEvaluation toCopy){
        super(toCopy);
        this.truncate = toCopy.truncate;
    }
    
    // impose exponentially increasing emphasis on low frequency favourable alleles
    @Override
    public double getValue() {
        double[] alleleFreqs = getAlleleFreqs();
        double val = 0.0;
        int numMarkers = alleleFreqs.length;
        for(int m=0; m<numMarkers; m++){
            double p = getFavourableAlleleFreq(m);
            // truncate at very low value to avoid -infinity from
            // which a local search might not be able to escape
            val += Math.max(Math.log(p), truncate);
        }
        val = 1.0/numMarkers * val;
        return val;
    }
    
}
