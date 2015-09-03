

package org.jamesframework.gs.simulation.obj;

import java.util.function.BiFunction;
import org.jamesframework.core.exceptions.IncompatibleDeltaEvaluationException;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.search.neigh.Move;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.eval.LOGEvaluation;

/**
 * Computes averaged logarithm of favourable allele frequencies to put exponentially
 * increasing emphasis on low frequency alleles. See Li et al., Selection on multiple
 * QTL with control of gene diversity and inbreeding for long-term benefit, 2008, J.
 * Anim. Breed. Genet.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class LOGFrequency extends FrequencyBasedObjective{

    // function used to produce the value at zero:
    //  - two integer arguments: (1) subset size, (2) number of markers
    //  - double results (value at zero)
    private final BiFunction<Integer, Integer, Double> valueAtZeroFunction;

    public LOGFrequency(BiFunction<Integer, Integer, Double> valueAtZeroFunction) {
        this.valueAtZeroFunction = valueAtZeroFunction;
    }
    
    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection size
        int n = solution.getNumSelectedIDs();
        // compute allele frequencies
        double[] avgMarkers = computeAlleleFrequencies(solution, data);
        
        // wrap in LOG evaluation
        return new LOGEvaluation(n, avgMarkers, data.getFavourableAlleles(), valueAtZeroFunction.apply(n, data.numMarkers()));
        
    }
    
    @Override
    public Evaluation evaluate(Move move, SubsetSolution curSolution, Evaluation curEvaluation, PopulationData data){
        
        // check move type
        if(!(move instanceof SwapMove)){
            throw new IncompatibleDeltaEvaluationException("LOG frequency objective should be used in "
                                                         + "combination with neighbourhoods that generate swap moves.");
        }
        // cast move
        SwapMove swapMove = (SwapMove) move;

        // initialize new evaluation
        LOGEvaluation newEval = new LOGEvaluation((LOGEvaluation) curEvaluation);
        // update
        updateEvaluation(newEval, swapMove, data);
        
        return newEval;
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }

}