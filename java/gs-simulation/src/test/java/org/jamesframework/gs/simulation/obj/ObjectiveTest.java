
package org.jamesframework.gs.simulation.obj;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Random;
import java.util.Scanner;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.junit.BeforeClass;

public class ObjectiveTest {

    protected static final Random RG = new Random();
    protected static PopulationData DATA;
    
    @BeforeClass
    public static void beforeClass() throws IOException {
        
        // generate random DATA
        
        int individuals = 200;
        int markers = 800;
        
        String[] names = new String[individuals];
        double[] values = new double[individuals];
        int[][] markerdata = new int[individuals][markers];
        int[] favAlleles = new int[markers];
        
        for(int i=0; i<individuals; i++){
            names[i] = "ind-" + i;
            values[i] = RG.nextDouble();
            for(int m=0; m<markers; m++){
                markerdata[i][m] = RG.nextBoolean() ? 1 : 0;
            }
        }
        for(int m=0; m<markers; m++){
            favAlleles[m] = RG.nextBoolean() ? 1 : 0;
        }
        
        // read arbitrary G matrix
        double[][] G = new double[individuals][individuals];
        Scanner sc = new Scanner(Paths.get(ObjectiveTest.class.getResource("/G.txt").getPath()));
        int i = 0, j = 0;
        while(sc.hasNext()){
            G[i][j] = sc.nextDouble();
            j++;
            if(j == individuals){
                i++;
                j = 0;
            }
        }
        
        DATA = new PopulationData(names, values, markerdata, favAlleles, G);
        
    }
    
}
