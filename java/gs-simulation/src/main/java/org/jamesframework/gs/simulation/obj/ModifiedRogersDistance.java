

package org.jamesframework.gs.simulation.obj;

/**
 * Modified Rogers' genetic distance.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class ModifiedRogersDistance implements GeneticDistanceFunction {

    @Override
    public double computeDistance(int[] markers1, int[] markers2) {
        double sum = 0.0;
        int n = markers1.length;
        for(int m=0; m < n; m++){
            // abs instead of square (same results for 0/1 markers, but faster)
            sum += Math.abs(markers1[m] - markers2[m]);
        }
        return Math.sqrt(sum/n);
    }

}
