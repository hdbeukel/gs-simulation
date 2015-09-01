

package org.jamesframework.gs.simulation.obj.eval;

/**
 * Evaluation based on favourable allele frequencies in the selection.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public abstract class FavourableAlleleBasedEvaluation extends FrequencyBasedEvaluation {

    // indicates which allele is favourable
    private final int[] favAlleles;
    
    public FavourableAlleleBasedEvaluation(int n, double[] alleleFreqs, int[] favAlleles) {
        super(n, alleleFreqs);
        this.favAlleles = favAlleles;
    }
    
    public FavourableAlleleBasedEvaluation(FavourableAlleleBasedEvaluation toCopy){
        super(toCopy);
        this.favAlleles = toCopy.favAlleles; // no deep copy (doesn't change across selections)
    }
    
    public double getFavourableAlleleFreq(int m){
        // get frequency of 1-allele at marker m
        double freq = getAlleleFreqs()[m];
        // flip if 0-allele is favourable
        if(favAlleles[m] == 0){
            freq = 1.0 - freq;
        }
        return freq;
    }
    
    public double getUnfavourableAlleleFreq(int m){
        return 1.0 - getFavourableAlleleFreq(m);
    }
    
}
