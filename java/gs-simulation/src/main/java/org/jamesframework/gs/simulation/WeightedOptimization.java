/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.jamesframework.gs.simulation;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.search.Search;
import org.jamesframework.core.search.algo.ParallelTempering;
import org.jamesframework.core.search.neigh.Neighbourhood;
import org.jamesframework.core.search.stopcriteria.MaxRuntime;
import org.jamesframework.core.subset.SubsetProblem;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.SingleSwapNeighbourhood;
import org.jamesframework.examples.util.ProgressSearchListener;
import org.jamesframework.ext.problems.objectives.WeightedIndex;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.data.PopulationReader;
import org.jamesframework.gs.simulation.obj.EntryToNearestEntryDistance;
import org.jamesframework.gs.simulation.obj.ExpectedProportionOfHeterozygotes;
import org.jamesframework.gs.simulation.obj.MedianAccessionValue;
import org.jamesframework.gs.simulation.obj.ModifiedRogersDistance;

public class WeightedOptimization {

    public static void main(String[] args) throws IOException {
        
        // get arguments
        String valueFile = args[0];
        String markerFile = args[1];
        String divMeasure = args[2].toUpperCase();
        double divWeight = Double.parseDouble(args[3]);
        double valueWeight = 1.0 - divWeight;
        int subsetSize = Integer.parseInt(args[4]);
        int timeLimit = Integer.parseInt(args[5]);
        // read data
        PopulationReader reader = new PopulationReader();
        PopulationData data = reader.read(Paths.get(valueFile), Paths.get(markerFile));
        // precompute MR distance matrix
        data.precomputeDistanceMatrix(new ModifiedRogersDistance());

        // create diversity objective
        System.out.println("Diversity measure: " + divMeasure);
        Objective<SubsetSolution, PopulationData> divObj;
        switch(divMeasure){
            case "HE": divObj = new ExpectedProportionOfHeterozygotes();
                break;
            case "ENE": divObj = new EntryToNearestEntryDistance();
                break;
            default: throw new IllegalArgumentException("Unknown diversity measure: " + divMeasure + ". "
                                                      + "Please specify one of HE or ENE.");
        }

        // create value objective
        Objective<SubsetSolution, PopulationData> valueObj = new MedianAccessionValue();

        // compose weighted index
        System.out.println("Diversity weight: " + divWeight);
        System.out.println("Quality weight: " + valueWeight);
        WeightedIndex<SubsetSolution, PopulationData> index = new WeightedIndex<>();
        if(divWeight > 1e-8){
            index.addObjective(divObj, divWeight);
        }
        if(valueWeight > 1e-8){
            index.addObjective(valueObj, valueWeight);
        }
        
        // create problem
        SubsetProblem<PopulationData> problem = new SubsetProblem<>(data, index, subsetSize);
        
        // create single swap neighbourhod
        Neighbourhood<SubsetSolution> neigh = new SingleSwapNeighbourhood();
        // create parallel tempering search
        double minTemp = 1e-8;
        double maxTemp = 1e-3;
        int numReplicas = 10;
        Search<SubsetSolution> search = new ParallelTempering<>(problem, neigh,
                                                                numReplicas,
                                                                minTemp, maxTemp);
        // set maximum runtime
        search.addStopCriterion(new MaxRuntime(timeLimit, TimeUnit.SECONDS));
        // track progress
        search.addSearchListener(new ProgressSearchListener());
        
        // run search
        search.start();
        // output results
        SubsetSolution bestSol = search.getBestSolution();
        Integer[] selection = bestSol.getSelectedIDs().toArray(new Integer[0]);
        Arrays.sort(selection);
        Evaluation bestEval = search.getBestSolutionEvaluation();
        Evaluation divEval = divObj.evaluate(bestSol, data);
        Evaluation valueEval = valueObj.evaluate(bestSol, data);
        System.out.println("Final selection: " + Arrays.toString(selection));
        System.out.println("Best solution value: " + bestEval.getValue());
        System.out.println("Diversity score: " + divEval.getValue());
        System.out.println("Median accession value: " + valueEval.getValue());
        
        // dispose search
        search.dispose();

    }
    
}
