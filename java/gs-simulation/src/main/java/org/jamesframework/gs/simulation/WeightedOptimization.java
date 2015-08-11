

package org.jamesframework.gs.simulation;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.gs.simulation.api.API;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.data.PopulationReader;
import org.jamesframework.gs.simulation.obj.EntryToNearestEntryDistance;
import org.jamesframework.gs.simulation.obj.ExpectedProportionOfHeterozygousLoci;
import org.jamesframework.gs.simulation.obj.ModifiedRogersDistance;

public class WeightedOptimization {

    public static void main(String[] args) throws IOException {
        
        API api = new API(true);
        
        // get arguments
        String valueFile = args[0];
        String markerFile = args[1];
        String divMeasure = args[2].toUpperCase();
        double divWeight = Double.parseDouble(args[3]);
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

        String[] selectedNames = api.selectWeighted(subsetSize, data, divWeight, divObj, timeWithoutImpr);
        Arrays.sort(selectedNames);
        
        System.out.println("Selected names: " + Arrays.toString(selectedNames));

    }
    
}
