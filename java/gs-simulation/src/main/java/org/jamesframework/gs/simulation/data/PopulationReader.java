
package org.jamesframework.gs.simulation.data;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;

/**
 * Read population data.
 *  
 * @author <a href="mailto:herman.debeukelaer@ugent.be">Herman De Beukelaer</a>
 */
public class PopulationReader {

    public PopulationData read(Path valueFile, Path markerFile) throws IOException {
        
        // read and temporarily store values
        Map<String, Double> valueMap = new HashMap<>();
        Scanner sc = new Scanner(valueFile);
        sc.useLocale(Locale.US);
        while(sc.hasNext()){
            String name = stripQuotes(sc.next());
            double value = sc.nextDouble();
            valueMap.put(name, value);
        }
        // get number of markers
        sc = new Scanner(markerFile);
        sc.useLocale(Locale.US);
        int numMarkers = sc.nextLine().split(" ").length;
        
        // initialize marker matrix, name array and value array
        int numAccessions = valueMap.size();
        String[] names = new String[numAccessions];
        double[] values = new double[numAccessions];
        int[][] markers = new int[numAccessions][numMarkers];
        
        // fill data structures ordered by accession name
        int i=0;
        while(sc.hasNext()){
            // read name of next accession in marker matrix
            String name = stripQuotes(sc.next());
            // store name and corresponding value
            names[i] = name;
            values[i] = valueMap.get(name);
            // read and store marker data
            for(int m=0; m<numMarkers; m++){
                markers[i][m] = sc.nextInt();
            }
            i++;
        }
        
        // initialize data (no favourable allele references)
        PopulationData data = new PopulationData(names, values, markers, null);
        
        return data;
        
    }
    
    private String stripQuotes(String str){
        return str.replaceAll("^\"|\"$", "");
    }
    
}
