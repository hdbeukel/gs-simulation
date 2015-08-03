
package org.jamesframework.gs.simulation.obj;

import java.util.Random;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.core.util.SetUtilities;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class ExpectedProportionOfHeterozygousLociTest {

    private static final Random RG = new Random();
    
    private static PopulationData data;
    
    @BeforeClass
    public static void beforeClass() {
        
        // generate random data
        
        int individuals = 200;
        int markers = 800;
        
        String[] names = new String[individuals];
        double[] values = new double[individuals];
        int[][] markerdata = new int[individuals][markers];
        
        for(int i=0; i<individuals; i++){
            names[i] = "ind-" + i;
            values[i] = RG.nextDouble();
            for(int m=0; m<markers; m++){
                markerdata[i][m] = RG.nextBoolean() ? 1 : 0;
            }
        }
        
        data = new PopulationData(names, values, markerdata);
        data.precomputeDistanceMatrix(new ModifiedRogersDistance());
        
    }
    
    /**
     * Test delta evaluation.
     */
    @Test
    public void testDeltaEvaluation() {
        System.out.println(" - test delta evaluation");
       
        // create objective
        ExpectedProportionOfHeterozygousLoci obj = new ExpectedProportionOfHeterozygousLoci();
        
        // sample random subsets
        for (int i=0; i<1000; i++){
            SubsetSolution sol = new SubsetSolution(data.getIDs(),
                                                    SetUtilities.getRandomSubset(data.getIDs(), 20, RG));
            // full evaluation
            Evaluation eval = obj.evaluate(sol, data);
            // apply series of random swap moves
            for(int m=0; m<100; m++){
                int add = SetUtilities.getRandomElement(sol.getUnselectedIDs(), RG);
                int del = SetUtilities.getRandomElement(sol.getSelectedIDs(), RG);
                SwapMove move = new SwapMove(add, del);
                // delta eval
                Evaluation deltaEval = obj.evaluate(move, sol, eval, data);
                // apply-undo full eval
                move.apply(sol);
                Evaluation fullEval = obj.evaluate(sol, data);
                move.undo(sol);
                // verify
                assertEquals(fullEval.getValue(), deltaEval.getValue(), 1e-8);
            }
        }
        
    }
    
}
