
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

    public PopulationData read(Path valueFile, Path markerFile, Path favAlleleFile, Path GFile) throws IOException {
        
        Scanner sc;
        
        // read and temporarily store individual values
        Map<String, Double> valueMap = new HashMap<>();
        sc = new Scanner(valueFile);
        sc.useLocale(Locale.US);
        while(sc.hasNext()){
            String name = stripQuotes(sc.next());
            double value = sc.nextDouble();
            valueMap.put(name, value);
        }
        // read and temporarily store marker favourable alleles
        Map<String, Integer> favAlleleMap = new HashMap<>();
        sc = new Scanner(favAlleleFile);
        sc.useLocale(Locale.US);
        while(sc.hasNext()){
            String name = stripQuotes(sc.next());
            int favAllele = sc.nextInt();
            favAlleleMap.put(name, favAllele);
        }
        
        // initialize marker matrix + name, value and fav allele arrays, G matrix
        int numAccessions = valueMap.size();
        int numMarkers = favAlleleMap.size();
        int[][] markers = new int[numAccessions][numMarkers];
        String[] names = new String[numAccessions];
        double[] values = new double[numAccessions];
        int[] favAlleles = new int[numMarkers];
        double[][] G = new double[numAccessions][numAccessions];
        
        // read marker names in order in which they occur in marker matrix
        sc = new Scanner(markerFile);
        sc.useLocale(Locale.US);
        String[] markerNames = sc.nextLine().split(" ");
        
        // fill favourable allele array according to column order in marker matrix
        for(int m=0; m<numMarkers; m++){
            favAlleles[m] = favAlleleMap.get(stripQuotes(markerNames[m]));
        }
        
        // fill marker matrix + individual values according to row order in marker matrix
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
        
        // read G matrix
        sc = new Scanner(GFile);
        sc.useLocale(Locale.US);
        int r = 0, c = 0;
        while(sc.hasNext()){
            G[r][c] = sc.nextDouble();
            c++;
            if(c == numAccessions){
                r++;
                c = 0;
            }
        }
        
        // initialize data
        PopulationData data = new PopulationData(names, values, markers, favAlleles, G);
        
        return data;
        
    }
    
    private String stripQuotes(String str){
        return str.replaceAll("^\"|\"$", "");
    }
    
}
