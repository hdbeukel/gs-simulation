

package org.jamesframework.gs.simulation;

public class CLI {

    public static void main(String[] args) {
        System.out.println("Run one of the following commands:");
        System.out.println();
        System.out.println("java -cp gs-simulation.jar org.jamesframework.gs.simulation.ParetoFrontApproximation "
                         + "<value-file> <marker-file> <fav-allele-file> <G-file> <diversity-measure> <weight-delta> "
                         + "<subset-size> <repeats> <time-limit-without-impr> <normalized>");
        System.out.println("java -cp gs-simulation.jar org.jamesframework.gs.simulation.WeightedOptimization "
                         + "<value-file> <marker-file> <fav-allele-file> <G-file> <diversity-measure> "
                         + "<diversity-weight> <subset-size> <time-without-impr>");
        System.out.println();
    }
    
}
