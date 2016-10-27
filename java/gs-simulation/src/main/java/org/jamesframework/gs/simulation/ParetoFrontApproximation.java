

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
import org.jamesframework.gs.simulation.obj.AverageCoancestry;
import org.jamesframework.gs.simulation.obj.HEfav;
import org.jamesframework.gs.simulation.obj.HEall;
import org.jamesframework.gs.simulation.obj.LOGall;
import org.jamesframework.gs.simulation.obj.LOGfav;
import org.jamesframework.gs.simulation.obj.MeanBreedingValue;

public class ParetoFrontApproximation {

    public static void main(String[] args) throws IOException {
        
        API api = new API();
        
        // get arguments
        String valueFile = args[0];
        String markerFile = args[1];
        String favAlleleFile = args[2];
        String GFile = args[3];
        String divMeasure = args[4].toUpperCase();
        double weightDelta = Double.parseDouble(args[5]);
        int subsetSize = Integer.parseInt(args[6]);
        int repeats = Integer.parseInt(args[7]);
        int timeWithoutImpr = Integer.parseInt(args[8]);
        boolean normalized = Boolean.parseBoolean(args[9]);
        // read data
        PopulationReader reader = new PopulationReader();
        PopulationData data = reader.read(
                Paths.get(valueFile), 
                Paths.get(markerFile), 
                Paths.get(favAlleleFile),
                Paths.get(GFile)
        );

        // create diversity objective
        Objective<SubsetSolution, PopulationData> divObj;
        switch(divMeasure){
            case "HEALL": divObj = new HEall();
                break;
            case "HEFAV": divObj = new HEfav();
                break;
            case "LOGALL": divObj = new LOGall();
                break;
            case "LOGFAV": divObj = new LOGfav();
                break;
            case "OC": divObj = new AverageCoancestry();
                break;
            default: throw new IllegalArgumentException("Unknown diversity measure: " + divMeasure + ". "
                                                      + "Please specify one of HEall, HEfav, LOGall, LOGfav or OC.");
        }
        
        Objective<SubsetSolution, PopulationData> valObj = new MeanBreedingValue();
        NormalizedObjectives normObjs = null;
        if(normalized){
            // normalize objectives
            normObjs = api.getNormalizedObjectives(divObj, data, subsetSize);
        }

        // run optimizations with different weighted objectives
        if(normalized){
            System.out.println("ID, repeat, divWeight, valueWeight, div, value, div (normalized), value (normalized), weighted (normalized)");
        } else {
            System.out.println("ID, repeat, divWeight, valueWeight, div, value, weighted");
        }

        double divWeight = 0.0, valueWeight;
        int numExperiments = (int)(Math.round(1/weightDelta)) + 1;
        for(int experimentID = 0; experimentID < numExperiments; experimentID++){

            // set value weight
            valueWeight = 1.0 - divWeight;

            // create weighted index
            List<Objective<SubsetSolution, PopulationData>> objs = Arrays.asList(
                    normObjs != null ? normObjs.getDivObj() : divObj,
                    normObjs != null ? normObjs.getValueObj() : valObj 
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
                Evaluation divEval = divObj.evaluate(bestSol, data);
                Evaluation valueEval = valObj.evaluate(bestSol, data);
                if(normObjs != null){
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
                } else {
                    System.out.format("%d, %d, %f, %f, %f, %f, %f\n",
                                      experimentID,
                                      r,
                                      divWeight,
                                      valueWeight,
                                      divEval.getValue(),
                                      valueEval.getValue(),
                                      weightedEval.getValue());
                }
                // dispose search
                search.dispose();
            }

            // update diversity weight and ID for next experiment
            divWeight += weightDelta;

        }

    }
    
}
