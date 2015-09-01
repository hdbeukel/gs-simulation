

package org.jamesframework.gs.simulation.obj;

import java.util.Set;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.eval.FrequencyBasedEvaluation;


public abstract class FrequencyBasedObjective implements Objective<SubsetSolution, PopulationData>{

    public double[] computeAlleleFrequencies(SubsetSolution sol, PopulationData data){
        
        // retrieve selection
        Set<Integer> selection = sol.getSelectedIDs();
        int n = selection.size();
        
        // compute allele frequencies
        int numMarkers = data.numMarkers();
        double[] avgMarkers = new double[numMarkers];
        selection.forEach(sel -> {
            // retrieve markers of selected accession
            int[] markers = data.getMarkers(sel);
            // aggregate with others
            for(int m=0; m<numMarkers; m++){
                avgMarkers[m] += markers[m];
            }
        });
        for(int m=0; m<numMarkers; m++){
            avgMarkers[m] = avgMarkers[m]/n;
        }
        
        return avgMarkers;
        
    }
    
    public void updateEvaluation(FrequencyBasedEvaluation eval, SwapMove move, PopulationData data){
        
        // extract genome of added/deleted individual
        int[] delGenome = data.getMarkers(move.getDeletedID());
        int[] addGenome = data.getMarkers(move.getAddedID());
        
        // update evaluation
        eval.swap(delGenome, addGenome);
        
    }
    
}
