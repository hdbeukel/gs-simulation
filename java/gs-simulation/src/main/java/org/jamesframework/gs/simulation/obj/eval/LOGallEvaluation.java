

package org.jamesframework.gs.simulation.obj.eval;

// LOG evaluation based on minor allele frequencies (general diversity)
public class LOGallEvaluation extends LOGEvaluation {
        
    public LOGallEvaluation(int subsetSize, double[] alleleFreqs) {
        // doesn't require to know favourable alleles
        super(subsetSize, alleleFreqs, null);
    }
    
    public LOGallEvaluation(LOGallEvaluation toCopy){
        super(toCopy);
    }
    
    @Override
    protected double getFrequency(int m){
        return getMinorAllelFreq(m);
    }
    
}
