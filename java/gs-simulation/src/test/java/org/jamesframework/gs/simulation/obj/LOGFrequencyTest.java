
package org.jamesframework.gs.simulation.obj;

import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.core.util.SetUtilities;
import org.junit.Test;
import static org.junit.Assert.*;

public class LOGFrequencyTest extends ObjectiveTest {
    
    /**
     * Test delta evaluation.
     */
    @Test
    public void testDeltaEvaluation() {
        System.out.println(" - test delta evaluation");
       
        // create objective
        LOGFrequency obj = new LOGFrequency((n,m) -> - (m * Math.log(n) + 1.0));
        
        // sample random subsets
        for (int i=0; i<1000; i++){
            SubsetSolution sol = new SubsetSolution(DATA.getIDs(),
                                                    SetUtilities.getRandomSubset(DATA.getIDs(), 20, RG));
            // full evaluation
            Evaluation eval = obj.evaluate(sol, DATA);
            // apply series of random swap moves
            for(int m=0; m<100; m++){
                int add = SetUtilities.getRandomElement(sol.getUnselectedIDs(), RG);
                int del = SetUtilities.getRandomElement(sol.getSelectedIDs(), RG);
                SwapMove move = new SwapMove(add, del);
                // delta eval
                Evaluation deltaEval = obj.evaluate(move, sol, eval, DATA);
                // apply-undo full eval
                move.apply(sol);
                Evaluation fullEval = obj.evaluate(sol, DATA);
                move.undo(sol);
                // verify
                assertEquals(fullEval.getValue(), deltaEval.getValue(), 1e-8);
            }
        }
        
    }
    
}
