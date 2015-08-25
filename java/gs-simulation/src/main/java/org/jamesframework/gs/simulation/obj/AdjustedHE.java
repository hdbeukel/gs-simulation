

package org.jamesframework.gs.simulation.obj;

import org.jamesframework.core.exceptions.IncompatibleDeltaEvaluationException;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.search.neigh.Move;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.eval.AdjustedHEEvaluation;

/**
 * Adjusted HE that only penalizes fixation of unfavourable alleles.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class AdjustedHE extends AvgGenomeObjective{

    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection size
        int n = solution.getNumSelectedIDs();
        // compute average genome
        double[] avgMarkers = computeAvgGenome(solution, data);
        
        // wrap in adjusted HE evaluation
        return new AdjustedHEEvaluation(n, avgMarkers, data.getFavourableAlleles());
        
    }
    
    @Override
    public Evaluation evaluate(Move move, SubsetSolution curSolution, Evaluation curEvaluation, PopulationData data){
        
        // check move type
        if(!(move instanceof SwapMove)){
            throw new IncompatibleDeltaEvaluationException("Expected proportion of heterozygous loci objective should be used in "
                                                         + "combination with neighbourhoods that generate swap moves.");
        }
        // cast move
        SwapMove swapMove = (SwapMove) move;

        // initialize new evaluation
        AdjustedHEEvaluation newEval = new AdjustedHEEvaluation((AdjustedHEEvaluation) curEvaluation);
        // update
        updateEvaluation(newEval, swapMove, data);
        
        return newEval;
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }

}
