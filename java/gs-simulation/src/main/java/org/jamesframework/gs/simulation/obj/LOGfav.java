

package org.jamesframework.gs.simulation.obj;

import org.jamesframework.core.exceptions.IncompatibleDeltaEvaluationException;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.search.neigh.Move;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.eval.LOGfavEvaluation;

/**
 * Computes averaged logarithm of FAVOURABLE allele frequencies to put exponentially
 * increasing emphasis on amplifying rare FAVOURABLE alleles. See Li et al., Selection
 * on multiple QTL with control of gene diversity and inbreeding for long-term benefit,
 * 2008, J. Anim. Breed. Genet.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class LOGfav extends FrequencyBasedObjective{
    
    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection size
        int n = solution.getNumSelectedIDs();
        // compute allele frequencies
        double[] avgMarkers = computeAlleleFrequencies(solution, data);
        
        // wrap in evaluation
        return new LOGfavEvaluation(n, avgMarkers, data.getFavourableAlleles());
        
    }
    
    @Override
    public Evaluation evaluate(Move move, SubsetSolution curSolution, Evaluation curEvaluation, PopulationData data){
        
        // check move type
        if(!(move instanceof SwapMove)){
            throw new IncompatibleDeltaEvaluationException("LOGfav objective should be used in "
                                                         + "combination with neighbourhoods that generate swap moves.");
        }
        // cast move
        SwapMove swapMove = (SwapMove) move;

        // initialize new evaluation
        LOGfavEvaluation newEval = new LOGfavEvaluation((LOGfavEvaluation) curEvaluation);
        // update
        updateEvaluation(newEval, swapMove, data);
        
        return newEval;
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }

}
