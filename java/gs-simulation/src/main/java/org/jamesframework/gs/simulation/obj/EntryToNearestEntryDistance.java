

package org.jamesframework.gs.simulation.obj;

import java.util.Collection;
import java.util.Collections;
import java.util.Set;
import org.jamesframework.core.exceptions.IncompatibleDeltaEvaluationException;
import org.jamesframework.core.problems.objectives.Objective;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;
import org.jamesframework.core.search.neigh.Move;
import org.jamesframework.core.subset.SubsetSolution;
import org.jamesframework.core.subset.neigh.moves.SwapMove;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.jamesframework.gs.simulation.obj.eval.EntryToNearestEntryEvaluation;

/**
 * Evaluates selection by computing the average distance from each entry to the closest
 * other entry. Requires that a distance matrix has been precomputed in the data class.
 * The computed value is to be maximized to retain diversity in the selection.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class EntryToNearestEntryDistance implements Objective<SubsetSolution, PopulationData>{
    
    @Override
    public Evaluation evaluate(SubsetSolution solution, PopulationData data) {
        // extract distance matrix
        double[][] distMatrix = data.getDistanceMatrix();
        // initialize evaluation object
        EntryToNearestEntryEvaluation eval = new EntryToNearestEntryEvaluation();
        // find closest neighbour of each item in the selection
        Set<Integer> selected = solution.getSelectedIDs();
        for(int sel : selected){
            // find closest other selected item
            int closestOther = findClosest(sel, selected, distMatrix);
            // register closest item in evaluation object
            eval.add(sel, closestOther, distMatrix[sel][closestOther]);
        }
        return eval;
    }
    
    private int findClosest(int item, Collection<Integer> group, double[][] distMatrix){
        return findClosest(item, group, Collections.emptySet(), distMatrix);
    }
    
    private int findClosest(int item, Collection<Integer> group, Collection<Integer> skip, double[][] distMatrix){
        double dist;
        Double minDist = null;
        Integer closestOther = null;
        for(int other : group){
            if(other != item && !skip.contains(other)){
                dist = distMatrix[item][other];
                if(minDist == null || dist < minDist){
                    minDist = dist;
                    closestOther = other;
                }
            }
        }
        return closestOther;
    }
    
    @Override
    public Evaluation evaluate(Move move, SubsetSolution curSolution, Evaluation curEvaluation, PopulationData data){
        // check move type
        if(!(move instanceof SwapMove)){
            throw new IncompatibleDeltaEvaluationException("Entry to nearest entry distance should be used in "
                                                  + "combination with neighbourhoods that generate swap moves.");
        }
        // cast move
        SwapMove swapMove = (SwapMove) move;

        // extract distance matrix
        double[][] dist = data.getDistanceMatrix();
        
        // cast evaluation (cannot fail as both evaluate methods return such evaluation object)
        EntryToNearestEntryEvaluation eval = (EntryToNearestEntryEvaluation) curEvaluation;
        // copy to initialize new evaluation
        EntryToNearestEntryEvaluation newEval = new EntryToNearestEntryEvaluation(eval);

        // get added and deleted ID from move
        int add = swapMove.getAddedID();
        int del = swapMove.getDeletedID();
        // get current selection from solution
        Set<Integer> curSelection = curSolution.getSelectedIDs();

        // update closest items (and find closest neighbour of newly added item)
        int addClosest = -1;
        double addClosestDistance = Double.MAX_VALUE;
        for(int item : curSelection){
            if(item == del){
                // discard contribution of removed item
                newEval.remove(item);
            } else {
                
                // retained item: retrieve current closest neighbour
                int curClosest = newEval.getClosest(item);
                if(dist[item][add] < dist[item][curClosest]){
                    // new item is closer: update
                    newEval.update(item, add, dist[item][add]);
                } else if (curClosest == del) {
                    // current closest is removed
                    // --> scan entire current selection
                    // --> NOTE! skip current closest item
                    int newClosest = findClosest(item, curSelection, Collections.singleton(curClosest), dist);
                    // check if added item is closer
                    if(dist[item][add] < dist[item][newClosest]){
                        newClosest = add;
                    }
                    // update
                    newEval.update(item, newClosest, dist[item][newClosest]);
                } else {
                    // current closest is retained and new item is *not* closer: do nothing
                }
                
                // update closest neighbour of newly added item
                if(dist[add][item] < addClosestDistance){
                    addClosestDistance = dist[add][item];
                    addClosest = item;
                }
                
            }
            
        }
        // register closest neighbour of newly added item
        newEval.add(add, addClosest, addClosestDistance);

        return newEval;
    }

    @Override
    public boolean isMinimizing() {
        return false;
    }
    
}
