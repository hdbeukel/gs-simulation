
package org.jamesframework.gs.simulation.obj;

import java.util.Random;
import org.jamesframework.gs.simulation.data.PopulationData;
import org.junit.BeforeClass;

public class ObjectiveTest {

    protected static final Random RG = new Random();
    protected static PopulationData DATA;
    
    @BeforeClass
    public static void beforeClass() {
        
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
        
        DATA = new PopulationData(names, values, markerdata, favAlleles);
        
    }
    
}
