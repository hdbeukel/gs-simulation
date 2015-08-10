

package org.jamesframework.gs.simulation.api;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeUnit;
import org.jamesframework.core.problems.Problem;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.search.Search;
import org.jamesframework.core.search.algo.ParallelTempering;
import org.jamesframework.core.search.algo.RandomDescent;
import org.jamesframework.core.search.neigh.Neighbourhood;
import org.jamesframework.core.search.stopcriteria.MaxRuntime;
import org.jamesframework.core.search.stopcriteria.MaxTimeWithoutImprovement;
import org.jamesframework.core.subset.SubsetProblem;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.SingleSwapNeighbourhood;
import org.jamesframework.ext.problems.objectives.WeightedIndex;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.MeanBreedingValue;
import org.jamesframework.gs.simulation.obj.NormalizedObjective;

public class API {
    
    private static final double FLOATING_POINT_TOL = 1e-10;
    
    // singleton instance
    private static final API API = new API();
    
    private API(){}
    
    public static API get(){
        return API;
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

        System.err.println("Normalizing objectives ...");
        
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
        
        System.err.format("Mean breeding value range: [%f, %f]\n", minVal, maxVal);
        System.err.format("Diversity score range: [%f, %f]\n", minDiv, maxDiv);
        
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
