

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
import org.jamesframework.examples.util.ProgressSearchListener;
import org.jamesframework.ext.problems.objectives.WeightedIndex;
import org.jamesframework.gs.simulation.api.API;
import org.jamesframework.gs.simulation.api.NormalizedObjectives;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.data.PopulationReader;
import org.jamesframework.gs.simulation.obj.EntryToNearestEntryDistance;
import org.jamesframework.gs.simulation.obj.ExpectedProportionOfHeterozygousLoci;
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
        int timeWithoutImpr = Integer.parseInt(args[5]);
        // read data
        PopulationReader reader = new PopulationReader();
        PopulationData data = reader.read(Paths.get(valueFile), Paths.get(markerFile));

        // create diversity objective
        System.out.println("Diversity measure: " + divMeasure);
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
        // create index
        System.out.println("Diversity weight: " + divWeight);
        System.out.println("Quality weight: " + valueWeight);
        List<Objective<SubsetSolution, PopulationData>> objs = Arrays.asList(
                normObjs.getDivObj(),
                normObjs.getValueObj()
        );
        List<Double> weights = Arrays.asList(divWeight, valueWeight);
        WeightedIndex<SubsetSolution, PopulationData> index = API.get().getWeightedIndex(objs, weights);
        
        // create problem
        SubsetProblem<PopulationData> problem = new SubsetProblem<>(data, index, subsetSize);
        
        // create parallel tempering search
        Search<SubsetSolution> search = API.get().createParallelTempering(problem);
        // set maximum runtime
        search.addStopCriterion(new MaxTimeWithoutImprovement(timeWithoutImpr, TimeUnit.SECONDS));
        // track progress
        search.addSearchListener(new ProgressSearchListener());
        
        // run search
        search.start();
        // output results (!! non-normalized values)
        SubsetSolution bestSol = search.getBestSolution();
        Integer[] selection = bestSol.getSelectedIDs().toArray(new Integer[0]);
        Arrays.sort(selection);
        Evaluation bestEval = search.getBestSolutionEvaluation();
        Evaluation divNormEval = normObjs.getDivObj().evaluate(bestSol, data);
        Evaluation valueNormEval = normObjs.getValueObj().evaluate(bestSol, data);
        Evaluation divEval = normObjs.getDivObj().getObjective().evaluate(bestSol, data);
        Evaluation valueEval = normObjs.getValueObj().getObjective().evaluate(bestSol, data);
        System.out.println("Final selection: " + Arrays.toString(selection));
        System.out.println("Best weighted value (normalized): " + bestEval.getValue());
        System.out.println("Diversity score (normalized): " + divNormEval.getValue());
        System.out.println("Mean breeding value (normalized): " + valueNormEval.getValue());
        System.out.println("Diversity score: " + divEval.getValue());
        System.out.println("Mean breeding value: " + valueEval.getValue());
        
        // dispose search
        search.dispose();

    }
    
}
