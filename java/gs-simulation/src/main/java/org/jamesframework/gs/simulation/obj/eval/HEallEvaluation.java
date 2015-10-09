

package org.jamesframework.gs.simulation.obj.eval;

public class HEallEvaluation extends FrequencyBasedEvaluation {

    public HEallEvaluation(int n, double[] alleleFreqs) {
        // doesn't require to know favourable alleles
        super(n, alleleFreqs, null);
    }
    
    public HEallEvaluation(HEallEvaluation toCopy){
        super(toCopy);
    }
    
    @Override
    public double getValue() {
        double he = 0.0;
        int numMarkers = getNumMarkers();
        for(int m = 0; m < numMarkers; m++){
            double p = getAlleleFreq(m);
            he += p * (1.0 - p);
        }
        he = 2.0 * he / numMarkers;
        return he;
    }
    
}
