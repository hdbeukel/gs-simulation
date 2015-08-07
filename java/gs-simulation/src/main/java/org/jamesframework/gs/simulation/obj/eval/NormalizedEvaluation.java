
package org.jamesframework.gs.simulation.obj.eval;

import org.jamesframework.core.problems.objectives.evaluations.Evaluation;

public class NormalizedEvaluation implements Evaluation {

    // wrapped evaluation to normalize
    private final Evaluation eval;
    // normalized value
    private final double normalizedValue;

    public NormalizedEvaluation(Evaluation eval, double min, double max) {
        // store evaluation
        this.eval = eval;
        // normalize value from [min,max] to [0,1] (linear)
        normalizedValue = (eval.getValue() - min) / (max - min);
    }

    @Override
    public double getValue() {
        return normalizedValue;
    }
    
    public Evaluation getEvaluation(){
        return eval;
    }
    
}
