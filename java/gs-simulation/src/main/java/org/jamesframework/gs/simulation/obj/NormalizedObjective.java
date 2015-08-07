

package org.jamesframework.gs.simulation.obj;

import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.problems.sol.Solution;
import org.jamesframework.core.search.neigh.Move;
import org.jamesframework.gs.simulation.obj.eval.NormalizedEvaluation;


public class NormalizedObjective<SolutionType extends Solution, DataType> implements Objective<SolutionType, DataType> {

    // wrapped objective
    private final Objective<? super SolutionType, ? super DataType> obj;
    // bounds of interval [min,max] which is normalized to [0,1] (linear)
    private final double min, max;

    public NormalizedObjective(Objective<? super SolutionType, ? super DataType> obj, double min, double max) {
        this.obj = obj;
        this.min = min;
        this.max = max;
    }

    public Objective<? super SolutionType, ? super DataType> getObjective() {
        return obj;
    }

    public double getMin() {
        return min;
    }

    public double getMax() {
        return max;
    }

    @Override
    public NormalizedEvaluation evaluate(SolutionType solution, DataType data) {
        // evaluate objective
        Evaluation eval = obj.evaluate(solution, data);
        // normalize
        NormalizedEvaluation normalizedEVal = new NormalizedEvaluation(eval, min, max);
        return normalizedEVal;
    }

    @Override
    public <ActualSolutionType extends SolutionType> NormalizedEvaluation evaluate(Move<? super ActualSolutionType> move,
                                                                                   ActualSolutionType curSolution,
                                                                                   Evaluation curEvaluation,
                                                                                   DataType data) {
        // cast evaluation
        NormalizedEvaluation normalizedEval = (NormalizedEvaluation) curEvaluation;
        // retrieve non-normalized evaluation
        Evaluation eval = normalizedEval.getEvaluation();
        // delta evaluation for non-normalized objective
        Evaluation newEval = obj.evaluate(move, curSolution, eval, data);
        // wrap in normalized evaluation
        NormalizedEvaluation newNormalizedEval = new NormalizedEvaluation(newEval, min, max);
        
        return newNormalizedEval;
        
    }
    
    @Override
    public boolean isMinimizing() {
        return obj.isMinimizing();
    }
    
}
