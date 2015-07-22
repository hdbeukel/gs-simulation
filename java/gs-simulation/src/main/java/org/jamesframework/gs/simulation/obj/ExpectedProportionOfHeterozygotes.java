/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.jamesframework.gs.simulation.obj;

import java.util.Set;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.problems.objectives.evaluations.SimpleEvaluation;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.gs.simulation.data.PopulationData;

/**
 * Evaluates selection by computing the expected proportion of heterozygotes.
 * This measures is also known as Nei's index of variation. The computed value
 * is to be maximized to retain diversity in the selection.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class ExpectedProportionOfHeterozygotes implements Objective<SubsetSolution, PopulationData>{

    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        
        // retrieve selection
        Set<Integer> selection = solution.getSelectedIDs();
        int n = selection.size();
        
        // compute average genome
        int numMarkers = data.numMarkers();
        double[] avgMarkers = new double[numMarkers];
        selection.forEach(sel -> {
            // retrieve markers of selected accession
            int[] markers = data.getMarkers(sel);
            // aggregate with others
            for(int m=0; m<numMarkers; m++){
                avgMarkers[m] += markers[m];
            }
        });
        for(int m=0; m<numMarkers; m++){
            avgMarkers[m] = avgMarkers[m]/n;
        }
        
        // compute expected proportion of heterozygotes
        double he = 0.0;
        for(int m=0; m<numMarkers; m++){
            he += avgMarkers[m] * (1.0 - avgMarkers[m]);
        }
        he = 4.0/numMarkers * he;
        
        // wrap in simple evaluation
        return SimpleEvaluation.WITH_VALUE(he);
        
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }

}
