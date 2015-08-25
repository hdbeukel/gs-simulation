

package org.jamesframework.gs.simulation.obj.eval;


public class AdjustedHEEvaluation extends AvgGenomeEvaluation {

    // indicates which allele is favourable
    private final int[] favAlleles;
    
    public AdjustedHEEvaluation(int n, double[] alleleFreqs, int[] favAlleles) {
        super(n, alleleFreqs);
        this.favAlleles = favAlleles;
    }
    
    public AdjustedHEEvaluation(AdjustedHEEvaluation toCopy){
        super(toCopy);
        this.favAlleles = toCopy.favAlleles; // no deep copy (doesn't change)
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
    
    private double getUnfavourableAlleleFreq(int m){
        // get frequency of 1-allele at marker m
        double freq = getAlleleFreqs()[m];
        // flip if 1-allele is favourable
        if(favAlleles[m] == 1){
            freq = 1.0 - freq;
        }
        return freq;
    }
    
}
