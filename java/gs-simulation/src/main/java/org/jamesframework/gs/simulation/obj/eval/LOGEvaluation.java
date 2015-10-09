

package org.jamesframework.gs.simulation.obj.eval;


public abstract class LOGEvaluation extends FrequencyBasedEvaluation {
    
    private final double valueAtZero;
    
    public LOGEvaluation(int subsetSize, double[] alleleFreqs, int[] favAlleles) {
        super(subsetSize, alleleFreqs, favAlleles);
        this.valueAtZero = -(Math.log(subsetSize) + 1.0);
    }
    
    public LOGEvaluation(LOGEvaluation toCopy){
        super(toCopy);
        this.valueAtZero = toCopy.valueAtZero;
    }
    
    // impose exponentially increasing emphasis on rare (favourable) alleles
    @Override
    public double getValue() {
        double val = 0.0;
        int numMarkers = getNumMarkers();
        for(int m = 0; m < numMarkers; m++){
            double p = getFrequency(m);
            if(p > 0.0){
                val += Math.log(p);
            } else {
                val += valueAtZero;
            }
        }
        val = val / numMarkers;
        return val;
    }
    
    // returns frequency on which evaluation is based
    // (MAF, favourable allele frequency, ...); depending
    // on implementation in sub-class
    abstract protected double getFrequency(int m);
    
}
