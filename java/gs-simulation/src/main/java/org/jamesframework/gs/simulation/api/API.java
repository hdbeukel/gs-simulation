

package org.jamesframework.gs.simulation.api;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeUnit;
import org.jamesframework.core.problems.Problem;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.search.Search;
import org.jamesframework.core.search.algo.ParallelTempering;
import org.jamesframework.core.search.algo.RandomDescent;
import org.jamesframework.core.search.neigh.Neighbourhood;
import org.jamesframework.core.search.stopcriteria.MaxTimeWithoutImprovement;
import org.jamesframework.core.subset.SubsetProblem;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.SingleSwapNeighbourhood;
import org.jamesframework.examples.util.ProgressSearchListener;
import org.jamesframework.ext.problems.objectives.WeightedIndex;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.MeanBreedingValue;
import org.jamesframework.gs.simulation.obj.ModifiedRogersDistance;
import org.jamesframework.gs.simulation.obj.NormalizedObjective;

public class API {
    
    private static final double FLOATING_POINT_TOL = 1e-10;
    
    private final boolean verbose;
    
    public API(){
        this(false);
    }
    
    public API(boolean verbose){
        this.verbose = verbose;
    }
    
    private void message(String mess){
        if(verbose){
            System.out.println(mess);
        }
    }
    
    /**
     * Select based on weighted index.
     * 
     * @param subsetSize subset size
     * @param values (estimated) breeding values
     * @param markers marker matrix Z (0/1)
     * @param divWeight weight of the diversity component (real number in [0,1])
     * @param divObj diversity objective (MR-ENE or HE)
     * @param secondsWithoutImprovement number of seconds without improvement after which the search stops
     * @return array of selected individuals (names)
     */
    public String[] selectWeighted(int subsetSize, String[] names, double[] values, int[][] markers,
                                   double divWeight, Objective<SubsetSolution, PopulationData> divObj,
                                   int secondsWithoutImprovement){
        
        // wrap data
        PopulationData data = new PopulationData(names, values, markers);
        
        return selectWeighted(subsetSize, data, divWeight, divObj, secondsWithoutImprovement);
        
    }
    
    public String[] selectWeighted(int subsetSize, PopulationData data,
                                   double divWeight, Objective<SubsetSolution, PopulationData> divObj,
                                   int secondsWithoutImprovement){
        
        message("############################");
        message("# WEIGHTED INDEX SELECTION #");
        message("############################");
        
        // precompute distance matrix
        data.precomputeDistanceMatrix(new ModifiedRogersDistance());
        // set value weight
        double valueWeight = 1.0 - divWeight;
        message("Diversity weight: " + divWeight);
        message("Breeding value weight: " + valueWeight);
        
        // normalize objectives
        NormalizedObjectives normObjs = getNormalizedObjectives(divObj, data, subsetSize);
        // create index
        List<Objective<SubsetSolution, PopulationData>> objs = Arrays.asList(
                normObjs.getDivObj(),
                normObjs.getValueObj()
        );
        List<Double> weights = Arrays.asList(divWeight, valueWeight);
        WeightedIndex<SubsetSolution, PopulationData> index = getWeightedIndex(objs, weights);
        
        // create problem
        SubsetProblem<PopulationData> problem = new SubsetProblem<>(data, index, subsetSize);
        
        // create parallel tempering search
        Search<SubsetSolution> search = createParallelTempering(problem);
        // set maximum runtime
        search.addStopCriterion(new MaxTimeWithoutImprovement(secondsWithoutImprovement, TimeUnit.SECONDS));
        // track progress (if verbose only)
        if(verbose){
            search.addSearchListener(new ProgressSearchListener());
        }
        
        // run search
        search.start();
        // output results
        SubsetSolution bestSol = search.getBestSolution();
        Integer[] selection = bestSol.getSelectedIDs().toArray(new Integer[0]);
        Arrays.sort(selection);
        Evaluation bestEval = search.getBestSolutionEvaluation();
        Evaluation divNormEval = normObjs.getDivObj().evaluate(bestSol, data);
        Evaluation valueNormEval = normObjs.getValueObj().evaluate(bestSol, data);
        Evaluation divEval = normObjs.getDivObj().getObjective().evaluate(bestSol, data);
        Evaluation valueEval = normObjs.getValueObj().getObjective().evaluate(bestSol, data);
        message("Final selection: " + Arrays.toString(selection));
        message("Best weighted value (normalized): " + bestEval.getValue());
        message("Diversity score (normalized): " + divNormEval.getValue());
        message("Mean breeding value (normalized): " + valueNormEval.getValue());
        message("Diversity score: " + divEval.getValue());
        message("Mean breeding value: " + valueEval.getValue());
        
        // dispose search
        search.dispose();
        
        return search.getBestSolution().getSelectedIDs().stream().map(data::getName).toArray(n -> new String[n]);
        
    }

    public Neighbourhood<SubsetSolution> getNeighbourhood(){
        return new SingleSwapNeighbourhood();
    }
    
    public Search<SubsetSolution> createRandomDescent(Problem<SubsetSolution> problem){
        // create random descent search
        Search<SubsetSolution> search = new RandomDescent<>(problem, getNeighbourhood());
        // return created search
        return search;
    }
    
    public Search<SubsetSolution> createParallelTempering(Problem<SubsetSolution> problem){
        // create parallel tempering search
        double minTemp = 1e-8;
        double maxTemp = 1e-4;
        int numReplicas = 10;
        Search<SubsetSolution> search = new ParallelTempering<>(problem,
                                                                getNeighbourhood(),
                                                                numReplicas, minTemp, maxTemp);
        // return created search
        return search;
    }
    
    public NormalizedObjectives getNormalizedObjectives(Objective<SubsetSolution, PopulationData> divObj,
                                                        PopulationData data, int subsetSize){

        message("Normalizing objectives ...");
        
        // STEP 1: infer subset with highest possible average value (sorting)
        
        // sort individual IDs on breeding value
        List<Integer> ids = new ArrayList<>(data.getIDs());
        ids.sort(Comparator.comparing(data::getValue));
        
        // create solution with highest value IDs
        SubsetSolution highestValSol = new SubsetSolution(data.getIDs());
        highestValSol.selectAll(ids.subList(ids.size() - subsetSize, ids.size()));

        // STEP 2: approximate subset with highest possible diversity (requires preliminary optimization)
        
        // create problem: maximize diversity
        SubsetProblem<PopulationData> maxDivProblem = new SubsetProblem<>(data, divObj, subsetSize);
        // run short parallel tempering search (stop after 3 seconds without improvement)
        Search<SubsetSolution> pt = createParallelTempering(maxDivProblem);
        pt.addStopCriterion(new MaxTimeWithoutImprovement(3, TimeUnit.SECONDS));
        pt.start();
        pt.dispose();
        // retrieve solution with maximized diversity and corresponding diversity score
        SubsetSolution highestDivSol = pt.getBestSolution();
        
        // STEP 3: get bounds of Pareto front
        
        MeanBreedingValue valueObj = new MeanBreedingValue();
        
        // get maximum value
        double maxVal = valueObj.evaluate(highestValSol, data).getValue();
        // get value of highest diversity solution
        double minVal = valueObj.evaluate(highestDivSol, data).getValue();
        // get maximum diversity
        double maxDiv = divObj.evaluate(highestDivSol, data).getValue();
        // get diversity of highest value solution
        double minDiv = divObj.evaluate(highestValSol, data).getValue();
        
        message(String.format("Mean breeding value range: [%f, %f]", minVal, maxVal));
        message(String.format("Diversity score range: [%f, %f]", minDiv, maxDiv));
        
        // STEP 4: create normalized objectives
        
        NormalizedObjective<SubsetSolution, PopulationData> normValueObj = new NormalizedObjective<>(valueObj, minVal, maxVal);
        NormalizedObjective<SubsetSolution, PopulationData> normDivObj = new NormalizedObjective<>(divObj, minDiv, maxDiv);
        
        NormalizedObjectives normObjs = new NormalizedObjectives(normValueObj, normDivObj);
        
        return normObjs;
        
    }
    
    public WeightedIndex<SubsetSolution, PopulationData> getWeightedIndex(List<? extends Objective<SubsetSolution, PopulationData>> objs,
                                                                          List<Double> weights){
        
        WeightedIndex<SubsetSolution, PopulationData> index = new WeightedIndex<>();
       
        for(int i=0; i<objs.size(); i++){
            double weight = weights.get(i);
            if(weight > FLOATING_POINT_TOL){
                index.addObjective(objs.get(i), weight);
            }
        }
        
        return(index);
        
    }
    
}
