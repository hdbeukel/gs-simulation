

package org.jamesframework.gs.simulation.obj.eval;

import java.util.ArrayList;
import java.util.List;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;


public class MedianEvaluation implements Evaluation {

    private final List<Double> values;

    public MedianEvaluation() {
        values = new ArrayList<>();
    }
    
    public MedianEvaluation(MedianEvaluation other){
        values = new ArrayList<>(other.values);
    }
    
    // insert into list so that order is maintained (ascending)
    public void add(double value){
        int i = 0;
        while(i < values.size() && values.get(i) < value){
            i++;
        }
        values.add(i, value);
    }
    
    public void remove(double value){
        values.remove(value);
    }
    
    public double getMedian(){
        int n = values.size();
        if(n % 2 == 1){
            // odd number of values: get middle one
            return values.get(n/2);
        } else {
            // even number of values: get average of the 2 central values
            double val1 = values.get(n/2 - 1);
            double val2 = values.get(n/2);
            return (val1 + val2) / 2;
        }
    }

    @Override
    public double getValue() {
        return getMedian();
    }
    
}
