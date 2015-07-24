

package org.jamesframework.gs.simulation.obj;

import java.util.Set;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.problems.objectives.evaluations.SimpleEvaluation;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.gs.simulation.data.PopulationData;

/**
 * Evaluates selection by computing the average distance from each entry to the closest
 * other entry. Requires that a distance matrix has been precomputed in the data class.
 * The computed value is to be maximized to retain diversity in the selection.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class EntryToNearestEntryDistance implements Objective<SubsetSolution, PopulationData>{

    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection and distance matrix
        Set<Integer> selection = solution.getSelectedIDs();
        double[][] dist = data.getDistanceMatrix();
        
        // sum distance from each item to closest other item
        double sum = selection.stream().mapToDouble(
            sel -> 
            // compute distance to closest other item
            selection.stream().filter(other -> (int) other != sel)
                              .mapToDouble(other -> dist[sel][other])
                              .min()
                              .getAsDouble()
        ).sum();
        
        // average and wrap in simple evaluation 
        return SimpleEvaluation.WITH_VALUE(sum/selection.size());
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }
    
}
