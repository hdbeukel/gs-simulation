/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.jamesframework.gs.simulation.obj;

/**
 * Interface of a distance function used to calculate the genetic distance between two accessions.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public interface GeneticDistanceFunction {
    
    public double computeDistance(int[] markers1, int[] markers2);
    
}
