

package org.jamesframework.gs.simulation.obj.eval;

public class HEEvaluation extends FrequencyBasedEvaluation {

    public HEEvaluation(int n, double[] alleleFreqs) {
        super(n, alleleFreqs);
    }
    
    public HEEvaluation(HEEvaluation toCopy){
        super(toCopy);
    }
    
    @Override
    public double getValue() {
        double[] alleleFreqs = getAlleleFreqs();
        double he = 0.0;
        int numMarkers = alleleFreqs.length;
        for(int m=0; m<numMarkers; m++){
            he += alleleFreqs[m] * (1.0 - alleleFreqs[m]);
        }
        he = 2.0/numMarkers * he;
        return he;
    }
    
}
