

package org.jamesframework.gs.simulation;

public class CLI {

    public static void main(String[] args) {
        System.out.println("Run one of the following commands:");
        System.out.println();
        System.out.println("java -cp gs-simulation.jar org.jamesframework.gs.simulation.ParetoFrontApproximation "
                         + "<value-file> <marker-file> <diversity-measure> <weight-delta> "
                         + "<subset-size> <repeats> <algo-time-limit>");
        System.out.println("java -cp gs-simulation.jar org.jamesframework.gs.simulation.WeightedOptimization "
                         + "<value-file> <marker-file> <diversity-measure> <diversity-weight> "
                         + "<subset-size> <time-limit>");
        System.out.println();
    }
    
}
