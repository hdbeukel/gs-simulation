
package org.jamesframework.gs.simulation.obj;

import java.util.Set;
import org.jamesframework.core.exceptions.IncompatibleDeltaEvaluationException;
import org.jamesframework.core.exceptions.SearchException;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.problems.objectives.evaluations.SimpleEvaluation;
import org.jamesframework.core.search.neigh.Move;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.gs.simulation.data.PopulationData;

/**
 * Evaluates the selection by computing the average coancestry of the selected individuals.
 * This value is to be minimized to limit inbreeding (obtained by maximizing -cGc/2).
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class AverageCoancestry implements Objective<SubsetSolution, PopulationData> {
    
    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection and size
        Set<Integer> selection = solution.getSelectedIDs();
        int n = solution.getNumSelectedIDs();
        // retrieve genomic relationship matrix
        double[][] G = data.getG();
        if(G == null){
            throw new SearchException("G matrix needed for average coancestry diversity measure!");
        }
        
        // compute average coancestry (including coancestry of plant with itself)
        double cGc = 0.0;
        for(int sel1 : selection){
            for(int sel2 : selection){
                cGc += G[sel1][sel2];
            }
        }
        cGc /= 2*n*n;
        
        // wrap in simple evaluation
        return SimpleEvaluation.WITH_VALUE(-cGc);
        
    }
    
    @Override
    public Evaluation evaluate(Move move, SubsetSolution curSolution, Evaluation curEval, PopulationData data){
        
        // check move type
        if(!(move instanceof SwapMove)){
            throw new IncompatibleDeltaEvaluationException("Average coancestry diversity measure should be used in "
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
        
        // undo average and minus
        val *= -2*n*n;
        
        // update
        double[][] G = data.getG();
        for(int sel : curSolution.getSelectedIDs()){
            if(sel != del){
                val -= 2 * G[sel][del];
                val += 2 * G[sel][add];
            } else {
                val -= G[sel][sel];
            }
        }
        val += G[add][add];
        
        // redo average and minus
        val /= -2*n*n;
        
        // wrap updated value in simple evaluation
        return SimpleEvaluation.WITH_VALUE(val);
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }

}
