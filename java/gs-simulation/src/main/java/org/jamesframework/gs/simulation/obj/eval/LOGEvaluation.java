

package org.jamesframework.gs.simulation.obj.eval;


public class LOGEvaluation extends FavourableAlleleBasedEvaluation {
    
    private static final double TRUNCATE = -10.0;
    
    public LOGEvaluation(int n, double[] alleleFreqs, int[] favAlleles) {
        super(n, alleleFreqs, favAlleles);
    }
    
    public LOGEvaluation(LOGEvaluation toCopy){
        super(toCopy);
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
            val += Math.max(Math.log(p), TRUNCATE);
        }
        val = 1.0/numMarkers * val;
        return val;
    }
    
}
