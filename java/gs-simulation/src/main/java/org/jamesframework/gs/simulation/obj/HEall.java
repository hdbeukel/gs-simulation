

package org.jamesframework.gs.simulation.obj;

import org.jamesframework.core.exceptions.IncompatibleDeltaEvaluationException;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.search.neigh.Move;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.eval.HEallEvaluation;

/**
 * Evaluates selection by computing the expected proportion of heterozygotes (HEall).
 * The computed value is to be maximized to retain diversity in the selection.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class HEall extends FrequencyBasedObjective{
    
    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection size
        int n = solution.getNumSelectedIDs();
        // compute average genome
        double[] avgMarkers = computeAlleleFrequencies(solution, data);
        
        // wrap in evaluation
        return new HEallEvaluation(n, avgMarkers);
        
    }
    
    @Override
    public Evaluation evaluate(Move move, SubsetSolution curSolution, Evaluation curEvaluation, PopulationData data){
        
        // check move type
        if(!(move instanceof SwapMove)){
            throw new IncompatibleDeltaEvaluationException("HEall objective should be used in combination"
                                                         + "with neighbourhoods that generate swap moves.");
        }
        // cast move
        SwapMove swapMove = (SwapMove) move;

        // initialize new evaluation
        HEallEvaluation newEval = new HEallEvaluation((HEallEvaluation) curEvaluation);
        // update
        updateEvaluation(newEval, swapMove, data);
        
        return newEval;
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }

}
