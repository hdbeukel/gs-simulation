

package org.jamesframework.gs.simulation.obj.eval;

// LOG evaluation based on favourable allele frequencies
public class LOGfavEvaluation extends LOGEvaluation {
        
    public LOGfavEvaluation(int subsetSize, double[] alleleFreqs, int[] favAlleles) {
        super(subsetSize, alleleFreqs, favAlleles);
    }
    
    public LOGfavEvaluation(LOGfavEvaluation toCopy){
        super(toCopy);
    }
    
    @Override
    protected double getFrequency(int m){
        return getFavourableAlleleFreq(m);
    }
    
}
