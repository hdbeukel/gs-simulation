
package org.jamesframework.gs.simulation.obj;

import java.util.Set;
import org.jamesframework.core.exceptions.IncompatibleDeltaEvaluationException;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.problems.objectives.evaluations.SimpleEvaluation;
import org.jamesframework.core.search.neigh.Move;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.gs.simulation.data.PopulationData;

/**
 * Evaluates the selection by computing the mean breeding value of the selected accessions.
 * This value is to be maximized to select good quality parents to be crossed.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class MeanBreedingValue implements Objective<SubsetSolution, PopulationData>{
    
    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection
        Set<Integer> selection = solution.getSelectedIDs();
        
        // compute average breeding value
        double avgValue = selection.stream().mapToDouble(data::getValue).average().getAsDouble();
        
        // wrap in simple evaluation
        return SimpleEvaluation.WITH_VALUE(avgValue);
        
    }
    
    @Override
    public Evaluation evaluate(Move move, SubsetSolution curSolution, Evaluation curEval, PopulationData data){
        
        // check move type
        if(!(move instanceof SwapMove)){
            throw new IncompatibleDeltaEvaluationException("Mean breeding value objective should be used in "
                                                  + "combination with neighbourhoods that generate swap moves.");
        }
        // cast move
        SwapMove swapMove = (SwapMove) move;
        
        // extract added and deleted item ID
        int add = swapMove.getAddedID();
        int del = swapMove.getDeletedID();
        
        // extract current value
        double val = curEval.getValue();
        
        // extract selection size
        int n = curSolution.getNumSelectedIDs();
        
        // update
        val += data.getValue(add)/n;
        val -= data.getValue(del)/n;
        
        // wrap updated value in simple evaluation
        return SimpleEvaluation.WITH_VALUE(val);
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }

}
