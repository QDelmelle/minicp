package minicp.examples;

import minicp.engine.core.Solver;

import static minicp.cp.Factory.makeSolver;
import minicp.examples.DARPParser.*;
import minicp.examples.DARPDataModel.*;

public class DARPtest {
    static DARPInstance instance;
    int searchTime;
    //    solPath: Option[String],
//    logPath: Option[String],
    boolean firstSolOnly = true;

    public static void main(String[] args) {
        Solver cp = makeSolver();
        instance = DARPParser.parseInstance("C:\\Users\\Utilisateur\\Documents\\UNIF2020\\TFE\\DARP RK\\DARP\\Cordeau\\a2-16.txt");
        System.out.println(instance.nVehicles);
    }
}
