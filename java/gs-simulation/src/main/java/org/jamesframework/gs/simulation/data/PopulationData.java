/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.jamesframework.gs.simulation.data;

import java.util.Arrays;
import java.util.DoubleSummaryStatistics;
import java.util.Locale;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.jamesframework.core.problems.datatypes.IntegerIdentifiedData;
import org.jamesframework.gs.simulation.obj.GeneticDistanceFunction;

/**
 * Population data consisting of (1) accession names, (2) accession values and
 * (3) 0/1 SNP marker data (DH).
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class PopulationData implements IntegerIdentifiedData {

    private final String[] names;
    private final double[] values;
    private final int[][] markers;
    
    private final Set<Integer> ids;
    
    // inter-accession distance matrix
    private double[][] dist;

    public PopulationData(String[] names, double[] values, int[][] markers) {
        // check sizes
        int n = names.length;
        if(values.length != n){
            throw new IllegalArgumentException("Number of names and values do not correspond.");
        }
        if(markers.length != n){
            throw new IllegalArgumentException("Number of names does not correspond with number of rows "
                                             + "in marker matrix.");
        }
        // store
        this.names = names;
        this.values = values;
        this.markers = markers;
        // infer IDs
        ids = IntStream.range(0, names.length)
                       .boxed()
                       .collect(Collectors.toSet());
    }

    @Override
    public Set<Integer> getIDs() {
        return ids;
    }
    
    // note: overwrites any previously computed distance matrix
    public void precomputeDistanceMatrix(GeneticDistanceFunction f){
        
        // initialize distance matrix
        dist = new double[numAccessions()][numAccessions()];
        // fill distance matrix
        for(int i=0; i<numAccessions(); i++){
            for(int j=0; j<i; j++){
                double d = f.computeDistance(markers[i], markers[j]);
                dist[i][j] = dist[j][i] = d;
            }
        }
        
    }
    
    // note: returns null if not precomputed
    public double[][] getDistanceMatrix(){
        return dist;
    }
    
    public void normalizeValues(){
        // infer current minimum and maximum value
        DoubleSummaryStatistics valueStats = Arrays.stream(values).summaryStatistics();
        double min = valueStats.getMin();
        double max = valueStats.getMax();
        for(int i=0; i<values.length; i++){
            values[i] = (values[i] - min) / (max - min);
        }
    }
    
    public String getName(int id){
        return names[id];
    }
    
    public double getValue(int id){
        return values[id];
    }
    
    public int[] getMarkers(int id){
        return markers[id];
    }
    
    public int numAccessions(){
        return ids.size();
    }
    
    public int numMarkers(){
        return markers[0].length;
    }
    
    @Override
    public String toString(){
        DoubleSummaryStatistics valueStats = Arrays.stream(values).summaryStatistics();
        return String.format(Locale.US,
                             "Population with %d accessions and %d markers",
                             numAccessions(), numMarkers());
    }
    
}
