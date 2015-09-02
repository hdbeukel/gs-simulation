

package org.jamesframework.gs.simulation.obj.eval;


public class LOGEvaluation extends FavourableAlleleBasedEvaluation {
    
    private final double valueAtZero;
    
    public LOGEvaluation(int subsetSize, double[] alleleFreqs, int[] favAlleles, double valueAtZero) {
        super(subsetSize, alleleFreqs, favAlleles);
        this.valueAtZero = valueAtZero;
    }
    
    public LOGEvaluation(LOGEvaluation toCopy){
        super(toCopy);
        this.valueAtZero = toCopy.valueAtZero;
    }
    
    // impose exponentially increasing emphasis on low frequency favourable alleles
    @Override
    public double getValue() {
        double[] alleleFreqs = getAlleleFreqs();
        double val = 0.0;
        int numMarkers = alleleFreqs.length;
        for(int m=0; m<numMarkers; m++){
            double p = getFavourableAlleleFreq(m);
            if(p > 0.0){
                val += Math.log(p);
            } else {
                val += valueAtZero;
            }
        }
        val = 1.0/numMarkers * val;
        return val;
    }
    
}
