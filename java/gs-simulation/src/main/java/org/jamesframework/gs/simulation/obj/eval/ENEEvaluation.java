
package org.jamesframework.gs.simulation.obj.eval;

import java.util.HashMap;
import java.util.Map;
import org.jamesframework.core.problems.objectives.evaluations.Evaluation;

public class ENEEvaluation implements Evaluation {

    // maps items to closest other items (IDs)
    private final Map<Integer, Integer> closestItemMap;
    // maps items to distance to respective closest item
    private final Map<Integer, Double> minDistMap;

    // sum of distances from items to respective closest items
    private double minDistSum;

    public ENEEvaluation() {
        closestItemMap = new HashMap<>();
        minDistMap = new HashMap<>();
        minDistSum = 0.0;
    }

    // deep copy constructor
    public ENEEvaluation(ENEEvaluation toCopy){
        closestItemMap = new HashMap<>(toCopy.closestItemMap);
        minDistMap = new HashMap<>(toCopy.minDistMap);
        minDistSum = toCopy.minDistSum;
    }

    // add item
    public void add(int itemID, int closestOtherItemID, double distance){
        // update minimum distance sum
        minDistSum += distance;
        // update metadata
        closestItemMap.put(itemID, closestOtherItemID);
        minDistMap.put(itemID, distance);
    }

    // remove item
    public boolean remove(int itemID){
        if(closestItemMap.containsKey(itemID)){
            // update minimum distance sum
            minDistSum -= minDistMap.get(itemID);
            // update metadata
            closestItemMap.remove(itemID);
            minDistMap.remove(itemID);
            return true;
        }
        return false;
    }

    // update closest item
    public boolean update(int itemID, int closestOtherItemID, double distance){
        if(closestItemMap.containsKey(itemID)){
            // update minimum distance sum
            minDistSum -= minDistMap.get(itemID);
            minDistSum += distance;
            // update metadata
            closestItemMap.put(itemID, closestOtherItemID);
            minDistMap.put(itemID, distance);
            return true;
        }
        return false;
    }

    // get closest item (null of no closest item registered)
    public Integer getClosest(int itemID){
        return closestItemMap.get(itemID);
    }

    // return average distance from each item to closest item; 0.0 if no distances
    @Override
    public double getValue() {
        int numDistances = minDistMap.size();
        if(numDistances > 0){
            return minDistSum/numDistances;
        } else {
            return 0.0;
        }
    }

}
