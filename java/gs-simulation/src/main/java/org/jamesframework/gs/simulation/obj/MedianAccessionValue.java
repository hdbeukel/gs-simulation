
package org.jamesframework.gs.simulation.obj;

import java.util.Set;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.problems.objectives.evaluations.SimpleEvaluation;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.gs.simulation.data.PopulationData;

/**
 * Evaluates the selection by computing the median value of the selected accessions.
 * This value is to be maximized to select good quality parents to be crossed.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class MedianAccessionValue implements Objective<SubsetSolution, PopulationData>{

    private static final Median MEDIAN = new Median();
    
    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection
        Set<Integer> selection = solution.getSelectedIDs();
        
        // retrieve values of selected accessions
        double[] values = selection.stream().mapToDouble(data::getValue).toArray();
        
        // compute median value
        double median = MEDIAN.evaluate(values);
        
        // wrap in simple evaluation
        return SimpleEvaluation.WITH_VALUE(median);
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }

}
