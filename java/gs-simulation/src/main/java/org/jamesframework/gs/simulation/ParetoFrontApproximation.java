

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
import org.jamesframework.gs.simulation.obj.AdjustedHE;
import org.jamesframework.gs.simulation.obj.EntryToNearestEntryDistance;
import org.jamesframework.gs.simulation.obj.ExpectedProportionOfHeterozygousLoci;
import org.jamesframework.gs.simulation.obj.ModifiedRogersDistance;

public class ParetoFrontApproximation {

    public static void main(String[] args) throws IOException {
        
        API api = new API();
        
        // get arguments
        String valueFile = args[0];
        String markerFile = args[1];
        String favAlleleFile = args[2];
        String divMeasure = args[3].toUpperCase();
        double weightDelta = Double.parseDouble(args[4]);
        int subsetSize = Integer.parseInt(args[5]);
        int repeats = Integer.parseInt(args[6]);
        int timeWithoutImpr = Integer.parseInt(args[7]);
        // read data
        PopulationReader reader = new PopulationReader();
        PopulationData data = reader.read(Paths.get(valueFile), Paths.get(markerFile), Paths.get(favAlleleFile));

        // create diversity objective
        Objective<SubsetSolution, PopulationData> divObj;
        switch(divMeasure){
            case "LOG":
                divObj = api.createLOGobjective();
                break;
            case "LOG2":
                divObj = api.createLOG2objective();
                break;
            case "HEADJ": divObj = new AdjustedHE();
                break;
            case "HE": divObj = new ExpectedProportionOfHeterozygousLoci();
                break;
            case "MR": divObj = new EntryToNearestEntryDistance();
                       // precompute MR distance matrix
                       data.precomputeDistanceMatrix(new ModifiedRogersDistance());
                break;
            default: throw new IllegalArgumentException("Unknown diversity measure: " + divMeasure + ". "
                                                      + "Please specify one of HE, MR, HEadj, LOG or LOG2.");
        }
        
        // normalize objectives
        NormalizedObjectives normObjs = api.getNormalizedObjectives(divObj, data, subsetSize);

        // run optimizations with different weighted objectives
        System.out.println("ID, repeat, divWeight, valueWeight, div, value, div (normalized), value (normalized), weighted (normalized)");

        double divWeight = 0.0, valueWeight;
        int numExperiments = (int)(Math.round(1/weightDelta)) + 1;
        for(int experimentID = 0; experimentID < numExperiments; experimentID++){

            // set value weight
            valueWeight = 1.0 - divWeight;

            // create weighted index
            List<Objective<SubsetSolution, PopulationData>> objs = Arrays.asList(
                    normObjs.getDivObj(),
                    normObjs.getValueObj()
            );
            List<Double> weights = Arrays.asList(divWeight, valueWeight);
            WeightedIndex<SubsetSolution, PopulationData> index = api.getWeightedIndex(objs, weights);

            // create problem
            SubsetProblem<PopulationData> problem = new SubsetProblem<>(data, index, subsetSize);

            // repeatedly run optimization
            for(int r=0; r<repeats; r++){
                // create search
                Search<SubsetSolution> search = api.createParallelTempering(problem);
                // set maximum runtime
                search.addStopCriterion(new MaxTimeWithoutImprovement(timeWithoutImpr, TimeUnit.SECONDS));
                // run search
                search.start();
                // output results
                SubsetSolution bestSol = search.getBestSolution();
                Evaluation weightedEval = search.getBestSolutionEvaluation();
                Evaluation divEval = normObjs.getDivObj().getUnnormalizedObjective().evaluate(bestSol, data);
                Evaluation valueEval = normObjs.getValueObj().getUnnormalizedObjective().evaluate(bestSol, data);
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

        }

    }
    
}
