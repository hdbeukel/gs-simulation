

package org.jamesframework.gs.simulation.api;

import org.jamesframework.core.problems.Problem;
import org.jamesframework.core.search.Search;
import org.jamesframework.core.search.algo.ParallelTempering;
import org.jamesframework.core.search.neigh.Neighbourhood;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.SingleSwapNeighbourhood;

public class API {
    
    // singleton instance
    private static final API API = new API();
    
    private API(){}
    
    public static API get(){
        return API;
    }

    public Search<SubsetSolution> createSearch(Problem<SubsetSolution> problem){
        // create single swap neighbourhod
        Neighbourhood<SubsetSolution> neigh = new SingleSwapNeighbourhood();
        // create parallel tempering search
        double minTemp = 1e-8;
        double maxTemp = 1e-4;
        int numReplicas = 10;
        Search<SubsetSolution> search = new ParallelTempering<>(problem, neigh,
                                                                numReplicas,
                                                                minTemp, maxTemp);
        // return created search
        return search;
    }
    
}
