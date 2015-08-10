

package org.jamesframework.gs.simulation;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.search.Search;
import org.jamesframework.core.search.stopcriteria.MaxTimeWithoutImprovement;
import org.jamesframework.core.subset.SubsetProblem;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.ext.problems.objectives.WeightedIndex;
import org.jamesframework.gs.simulation.api.API;
import org.jamesframework.gs.simulation.api.NormalizedObjectives;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.data.PopulationReader;
import org.jamesframework.gs.simulation.obj.EntryToNearestEntryDistance;
import org.jamesframework.gs.simulation.obj.ExpectedProportionOfHeterozygousLoci;
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
        int timeWithoutImpr = Integer.parseInt(args[6]);
        // read data
        PopulationReader reader = new PopulationReader();
        PopulationData data = reader.read(Paths.get(valueFile), Paths.get(markerFile));

        // create diversity objective
        Objective<SubsetSolution, PopulationData> divObj;
        switch(divMeasure){
            case "HE": divObj = new ExpectedProportionOfHeterozygousLoci();
                break;
            case "MR": divObj = new EntryToNearestEntryDistance();
                       // precompute MR distance matrix
                       data.precomputeDistanceMatrix(new ModifiedRogersDistance());
                break;
            default: throw new IllegalArgumentException("Unknown diversity measure: " + divMeasure + ". "
                                                      + "Please specify one of HE or MR.");
        }
        
        // normalize objectives
        NormalizedObjectives normObjs = API.get().getNormalizedObjectives(divObj, data, subsetSize);

        // run optimizations with different weighted objectives
        System.out.println("ID, repeat, divWeight, valueWeight, div, value, div (normalized), value (normalized), weighted (normalized)");

        double divWeight = 0.0, valueWeight;
        int experimentID = 0;
        while(divWeight <= 1.0){

            // set value weight
            valueWeight = 1.0 - divWeight;

            // create weighted index
            List<Objective<SubsetSolution, PopulationData>> objs = Arrays.asList(
                    normObjs.getDivObj(),
                    normObjs.getValueObj()
            );
            List<Double> weights = Arrays.asList(divWeight, valueWeight);
            WeightedIndex<SubsetSolution, PopulationData> index = API.get().getWeightedIndex(objs, weights);

            // create problem
            SubsetProblem<PopulationData> problem = new SubsetProblem<>(data, index, subsetSize);

            // repeatedly run optimization
            for(int r=0; r<repeats; r++){
                // create search
                Search<SubsetSolution> search = API.get().createParallelTempering(problem);
                // set maximum runtime
                search.addStopCriterion(new MaxTimeWithoutImprovement(timeWithoutImpr, TimeUnit.SECONDS));
                // run search
                search.start();
                // output results (!! non-normalized values)
                SubsetSolution bestSol = search.getBestSolution();
                Evaluation weightedEval = search.getBestSolutionEvaluation();
                Evaluation divEval = normObjs.getDivObj().getObjective().evaluate(bestSol, data);
                Evaluation valueEval = normObjs.getValueObj().getObjective().evaluate(bestSol, data);
                Evaluation divNormEval = normObjs.getDivObj().evaluate(bestSol, data);
                Evaluation valueNormEval = normObjs.getValueObj().evaluate(bestSol, data);
                System.out.format("%d, %d, %f, %f, %f, %f, %f, %f, %f\n",
                                  experimentID,
                                  r,
                                  divWeight,
                                  valueWeight,
                                  divEval.getValue(),
                                  valueEval.getValue(),
                                  divNormEval.getValue(),
                                  valueNormEval.getValue(),
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
