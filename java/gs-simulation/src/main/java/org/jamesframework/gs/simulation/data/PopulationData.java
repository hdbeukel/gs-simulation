

package org.jamesframework.gs.simulation.data;

import java.util.Locale;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.jamesframework.core.problems.datatypes.IntegerIdentifiedData;

/**
 * Population data consisting of (1) accession names, (2) accession values,
 * (3) 0/1 SNP marker data (DH), and (4) favourable allele reference.
 * 
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class PopulationData implements IntegerIdentifiedData {

    private final String[] names;
    private final double[] values;
    private final int[][] markers;
    private final int[] favourable;
        
    private final Set<Integer> ids;

    public PopulationData(String[] names, double[] values, int[][] markers, int[] favourable) {
        // check sizes
        int numInd = names.length;
        int numMarkers = markers[0].length;
        if(values.length != numInd){
            throw new IllegalArgumentException("Number of names and values do not correspond.");
        }
        if(markers.length != numInd){
            throw new IllegalArgumentException("Number of names does not correspond with number of rows "
                                             + "in marker matrix.");
        }
        if(favourable != null && favourable.length != numMarkers){
            throw new IllegalArgumentException("Length of favourable allele array does not correspond to "
                                             + "number of columns in marker matrix.");
        }
        // store
        this.names = names;
        this.values = values;
        this.markers = markers;
        this.favourable = favourable;
        // infer IDs
        ids = IntStream.range(0, names.length)
                       .boxed()
                       .collect(Collectors.toSet());
    }

    @Override
    public Set<Integer> getIDs() {
        return ids;
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
    
    public int[] getFavourableAlleles(){
        return favourable;
    }
    
    public int numAccessions(){
        return ids.size();
    }
    
    public int numMarkers(){
        return markers[0].length;
    }
    
    @Override
    public String toString(){
        return String.format(Locale.US,
                             "Population with %d accessions and %d markers",
                             numAccessions(), numMarkers());
    }
    
}
