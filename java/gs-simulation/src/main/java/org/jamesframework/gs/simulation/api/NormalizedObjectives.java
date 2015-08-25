
package org.jamesframework.gs.simulation.api;

import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.NormalizedObjective;

public class NormalizedObjectives {

    private final NormalizedObjective<SubsetSolution, PopulationData> valueObj, divObj;

    public NormalizedObjectives(NormalizedObjective<SubsetSolution, PopulationData> valueObj,
                                NormalizedObjective<SubsetSolution, PopulationData> divObj) {
        this.valueObj = valueObj;
        this.divObj = divObj;
    }

    public NormalizedObjective<SubsetSolution, PopulationData> getValueObj() {
        return valueObj;
    }

    public NormalizedObjective<SubsetSolution, PopulationData> getDivObj() {
        return divObj;
    }
    
}
