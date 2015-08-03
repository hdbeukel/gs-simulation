

package org.jamesframework.gs.simulation.obj;

import java.util.Set;
import org.jamesframework.core.exceptions.IncompatibleDeltaEvaluationException;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.search.neigh.Move;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.eval.HEEvaluation;

/**
 * Evaluates selection by computing the expected proportion of heterozygotes.
 * This measures is also known as Nei's index of variation. The computed value
 * is to be maximized to retain diversity in the selection.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class ExpectedProportionOfHeterozygousLoci implements Objective<SubsetSolution, PopulationData>{

    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection
        Set<Integer> selection = solution.getSelectedIDs();
        int n = selection.size();
        
        // compute average genome
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
        
        // wrap in HE evaluation
        return new HEEvaluation(n, avgMarkers);
        
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
        
        // cast current evaluation
        HEEvaluation eval = (HEEvaluation) curEvaluation;
        // initialize new evaluation
        HEEvaluation newEval = new HEEvaluation(eval);
        
        // extract genome of added/deleted individual
        int[] delGenome = data.getMarkers(swapMove.getDeletedID());
        int[] addGenome = data.getMarkers(swapMove.getAddedID());
        
        // update evaluation
        newEval.swap(delGenome, addGenome);
        
        return newEval;
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }

}
