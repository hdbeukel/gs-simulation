
package org.jamesframework.gs.simulation.obj;

import java.util.Set;
import org.jamesframework.core.exceptions.IncompatibleDeltaEvaluationException;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.search.neigh.Move;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.eval.MedianEvaluation;

/**
 * Evaluates the selection by computing the median value of the selected accessions.
 * This value is to be maximized to select good quality parents to be crossed.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class MedianAccessionValue implements Objective<SubsetSolution, PopulationData>{
    
    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection
        Set<Integer> selection = solution.getSelectedIDs();
        
        // initialize evaluation
        MedianEvaluation eval = new MedianEvaluation();
        
        // register values of selected accessions
        selection.stream().mapToDouble(data::getValue).forEach(eval::add);
        
        // return evaluation
        return eval;
        
    }
    
    @Override
    public Evaluation evaluate(Move move, SubsetSolution curSolution, Evaluation curEval, PopulationData data){
        
        // check move type
        if(!(move instanceof SwapMove)){
            throw new IncompatibleDeltaEvaluationException("Entry to nearest entry distance should be used in "
                                                  + "combination with neighbourhoods that generate swap moves.");
        }
        // cast move
        SwapMove swapMove = (SwapMove) move;
        
        // extract added and deleted item ID
        int add = swapMove.getAddedID();
        int del = swapMove.getDeletedID();
        
        // cast and copy evaluation object
        MedianEvaluation newEval = new MedianEvaluation((MedianEvaluation) curEval);
        
        // update
        newEval.remove(data.getValue(del));
        newEval.add(data.getValue(add));
        
        // return updated evaluation
        return newEval;
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }

}
