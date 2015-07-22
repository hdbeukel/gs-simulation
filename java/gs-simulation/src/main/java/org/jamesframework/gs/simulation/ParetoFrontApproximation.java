/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.jamesframework.gs.simulation;

import java.io.IOException;
import java.nio.file.Paths;
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
import org.jamesframework.ext.problems.objectives.WeightedIndex;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.data.PopulationReader;
import org.jamesframework.gs.simulation.obj.EntryToNearestEntryDistance;
import org.jamesframework.gs.simulation.obj.ExpectedProportionOfHeterozygotes;
import org.jamesframework.gs.simulation.obj.MedianAccessionValue;
import org.jamesframework.gs.simulation.obj.ModifiedRogersDistance;

public class ParetoFrontApproximation {

    public static void main(String[] args) throws IOException {
        
        // get arguments
        String valueFile = args[0];
        String markerFile = args[1];
        String divMeasure = args[2].toUpperCase();
        double weightDelta = Double.parseDouble(args[3]);
        int subsetSize = Integer.parseInt(args[4]);
        int repeats = Integer.parseInt(args[5]);
        int timeLimit = Integer.parseInt(args[6]);
        // read data
        PopulationReader reader = new PopulationReader();
        PopulationData data = reader.read(Paths.get(valueFile), Paths.get(markerFile));
        // precompute MR distance matrix
        data.precomputeDistanceMatrix(new ModifiedRogersDistance());

        // create diversity objective
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

        // run optimizations with different weighted objectives
        System.out.println("ID, repeat, divWeight, valueWeight, div, value, weighted");

        double divWeight = 0.0, valueWeight;
        int experimentID = 0;
        while(divWeight <= 1.0){

            // set value weight
            valueWeight = 1.0 - divWeight;

            // compose weighted index
            WeightedIndex<SubsetSolution, PopulationData> index = new WeightedIndex<>();
            if(divWeight > 1e-8){
                index.addObjective(divObj, divWeight);
            }
            if(valueWeight > 1e-8){
                index.addObjective(valueObj, valueWeight);
            }

            // create problem
            SubsetProblem<PopulationData> problem = new SubsetProblem<>(data, index, subsetSize);

            // repeatedly run optimization
            for(int r=0; r<repeats; r++){
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
                // run search
                search.start();
                // output results
                SubsetSolution bestSol = search.getBestSolution();
                Evaluation weightedEval = search.getBestSolutionEvaluation();
                Evaluation divEval = divObj.evaluate(bestSol, data);
                Evaluation valueEval = valueObj.evaluate(bestSol, data);
                System.out.format("%d, %d, %f, %f, %f, %f, %f\n",
                                  experimentID,
                                  r,
                                  divWeight,
                                  valueWeight,
                                  divEval.getValue(),
                                  valueEval.getValue(),
                                  weightedEval.getValue());
                // dispose search
                search.dispose();
            }

            // update diversity weight and ID for next experiment
            divWeight += weightDelta;
            experimentID++;

        }

    }
    
}
