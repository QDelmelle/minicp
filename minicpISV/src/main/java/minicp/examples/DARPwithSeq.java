
package minicp.examples;

import minicp.cp.Factory;
import minicp.engine.constraints.Circuit;
import minicp.engine.constraints.Element1D;
import minicp.engine.constraints.Element1DVar;
import minicp.engine.constraints.TableCT;
import minicp.engine.core.*;
import minicp.search.DFSearch;
import minicp.search.SearchStatistics;
import minicp.util.io.InputReader;

import java.util.*;
import java.util.stream.IntStream;

import static minicp.cp.BranchingScheme.and;
import static minicp.cp.BranchingScheme.firstFail;
import static minicp.cp.Factory.*;

/**
 * This problem consists in transporting patients to and from a hospital for medical appointments.
 * It is a specific kind of DARP:
 * The objective is to maximize the number of patients transported. A patient can be transported to an appointment,
 * brought back after the end of the appointment or both. In this last case, if one of the travels is done, the other
 * must be also performed but not necessarily by the same vehicle.
 * The patients can only be transported by compatible vehicles. The places from which and where the patient has to be
 * transported or the starting and ending depots of the vehicles are not necessarily the same.
 *
 * This model uses insertion sequence variables
 *
 * @author Quentin Delmelle qdelmelle@gmail.com
 * translated from Charles Thomas (cftmthomas@gmail.com)'s scala version.
 */

public class DARPwithSeq {
    /* data */
    int startTime = 0;
    int endTime = 400;
    int maxWaitTime = 30; //Maximum waiting time for patients before and after appointments.

    int nPatient = 8;
    int[] patientRdv = {60, 60, 120, 150, 240, 300, 300, 300}; //Rdv time
    int[] patientDur = {20, 120, 40, 40, 60, 20, 20, 60}; //Rdv duration
    int[] patientSrv = {2, 3, 2, 2, 5, 2, 2, 2}; //Service time (represents time needed to embark/disembark vehicle)
    int[] patientLoad = {1, 2, 1, 1, 1, 2, 1, 1}; //Number of places taken in the vehicle
    int[] patientCategory = {0, 1, 0, 0, 1, 0, 0, 0};

    int nVehicle = 2;
    int[] vehicleCap = {6, 4}; //Capacity of vehicles
    int[] vehicleStartDepot = {7, 8};
    int[] vehicleEndDepot = {7, 8};
    //Set<Integer>[] vehicleCompatibility = new HashSet<Integer>[2]; //Categories of patients that the vehicle can take
    int[] vehicleStartWindow = {0, 0}; //Start window of vehicle
    int[] vehicleEndWindow = {400, 400}; //End window of vehicle

    //Trip related variables, the forward and backward trips are separated:
    int[] patientForward = {0, 1, 2, 3, 4, 6, 7};
    int[] patientBackward = {0, 1, 2, 3, 4, 5, 6, 7};

    int[] originForward = {9, 10, 11, 13, 14, 16, 17};
    int[] destForward = {0, 1, 2, 3, 4, 6, 6};

    int[] originBackward = {0, 1, 2, 3, 4, 5, 6, 6};
    int[] destBackward = {9, 10, 12, 13, 14, 15, 16, 17};

    //Distance matrix, the places are identified by their id in this array:
    int[][] locationMatrix = {
            {0, 2, 6, 3, 3, 11, 2, 11, 9, 9, 0, 12, 9, 13, 9, 11, 10, 11},
            {2, 0, 7, 1, 0, 11, 2, 10, 8, 11, 2, 13, 11, 15, 6, 13, 9, 13},
            {6, 7, 0, 8, 8, 9, 7, 11, 10, 4, 6, 10, 5, 8, 13, 6, 13, 6},
            {3, 1, 8, 0, 1, 11, 3, 9, 7, 12, 3, 13, 12, 15, 5, 13, 8, 14},
            {3, 0, 8, 1, 0, 11, 2, 10, 7, 11, 2, 12, 11, 15, 6, 13, 8, 13},
            {11, 11, 9, 11, 11, 0, 12, 5, 6, 12, 11, 2, 13, 7, 14, 7, 9, 9},
            {2, 2, 7, 3, 2, 12, 0, 12, 10, 10, 2, 14, 10, 15, 8, 13, 11, 13},
            {11, 10, 11, 9, 10, 5, 12, 0, 2, 16, 11, 6, 17, 13, 10, 12, 4, 14},
            {9, 8, 10, 7, 7, 6, 10, 2, 0, 14, 8, 7, 15, 13, 8, 12, 3, 13},
            {9, 11, 4, 12, 11, 12, 10, 16, 14, 0, 9, 14, 1, 9, 17, 7, 17, 6},
            {0, 2, 6, 3, 2, 11, 2, 11, 8, 9, 0, 12, 9, 13, 8, 11, 10, 11},
            {12, 13, 10, 13, 12, 2, 14, 6, 7, 14, 12, 0, 15, 8, 15, 7, 10, 10},
            {9, 11, 5, 12, 11, 13, 10, 17, 15, 1, 9, 15, 0, 10, 17, 9, 18, 7},
            {13, 15, 8, 15, 15, 7, 15, 13, 13, 9, 13, 8, 10, 0, 20, 2, 16, 4},
            {9, 6, 13, 5, 6, 14, 8, 10, 8, 17, 8, 15, 17, 20, 0, 18, 7, 19},
            {11, 13, 6, 13, 13, 7, 13, 12, 12, 7, 11, 7, 9, 2, 18, 0, 14, 3},
            {10, 9, 13, 8, 8, 9, 11, 4, 3, 17, 10, 10, 18, 16, 7, 14, 0, 16},
            {11, 13, 6, 14, 13, 9, 13, 14, 13, 6, 11, 10, 7, 4, 19, 3, 16, 0}
    };

    public int distance(int i, int j) {
        if (i == -1 || j == -1) return 0;
        else return locationMatrix[i][j];
    }

    public static int[] arrayconcat(int[] a1, int[] a2) {
        int l = a1.length + a2.length;
        int[] ret = new int[l];
        for (int i = 0; i < a1.length; i++)
            ret[i] = a1[i];
        for (int i = a1.length; i < l; i++)
            ret[i] = a2[i];
        return ret;
    }

    protected class Stop {
        int patient;
        int place;
        int winStart;
        int winEnd;
        int operation;

        public Stop(int patient, int place, int winStart, int winEnd, int operation) {
            this.patient = patient;
            this.place = place;
            this.winStart = winStart;
            this.winEnd = winEnd;
            this.operation = operation;
        }

        boolean isDepot() {
            return patient == -1;
        }

        int category() {
            if (isDepot()) return -1;
            else return patientCategory[patient];
        }

        int service() {
            if (isDepot()) return 0;
            else return patientSrv[patient];
        }

        boolean forward() {
            return operation < 2;
        }

        boolean pickup() {
            return operation == 0 || operation == 2;
        }

        boolean isStartDepot() {
            return isDepot() && pickup();
        }

        boolean isEndDepot() {
            return isDepot() && !pickup();
        }

        int load() {
            if (isDepot()) return 0;
            else if (pickup()) return patientLoad[patient];
            else return -patientLoad[patient];
        }
    }

    public void main(String[] args) {
        int[] patientRdvEnd = new int[nPatient];
        for (int i = 0; i < nPatient; i++)
            patientRdvEnd[i] = patientDur[i] + patientRdv[i];

        int[] allPatient = arrayconcat(patientForward, patientBackward); //Patient for each trip
        int[] allOrigin = arrayconcat(originForward, originBackward); //Origin for each trip
        int[] allDest = arrayconcat(destForward, destBackward); //Destination for each trip

        ArrayList<Stop> stopBuffer = new ArrayList<Stop>();
        ArrayList<int[]> travelBuffer = new ArrayList<int[]>(); //stop1 (id), stop2 (id), forward

        //Generating stops for each patient:
        for (int i = 0; i < patientForward.length; i++) {
            int p = patientForward[i];
            int lstFor = patientRdv[p] - patientSrv[p] * 2 - distance(originForward[i], destForward[i]);

            if (patientRdv[p] < startTime || patientRdv[p] > endTime || lstFor < startTime) {
                System.out.println("Warning! Patient " + p + " is infeasible!");
            } else {
                int est = Math.max(startTime, Math.min(patientRdv[p] - maxWaitTime, lstFor));
                int ect = est + patientSrv[p] + distance(originForward[i], destForward[i]);
                int lct = patientRdv[p] - patientSrv[p];
                stopBuffer.add(new Stop(i, originForward[i], est, lstFor, 0));
                stopBuffer.add(new Stop(i, destForward[i], ect, lct, 1));
                travelBuffer.add(new int[]{stopBuffer.size() - 2, stopBuffer.size() - 1, 1});
            }
        }

        for (int i = 0; i < patientForward.length; i++) {
            int p = patientBackward[i];
            int ectBack = patientRdvEnd[p] + patientSrv[p] + distance(originBackward[i], destBackward[i]);

            if (patientRdv[p] < startTime || patientRdv[p] > endTime || ectBack > endTime) {
                System.out.println("Warning! Patient " + p + " is infeasible!");
            } else {
                int est = patientRdvEnd[p];
                int lct = Math.min(endTime, Math.max(est + maxWaitTime, ectBack));
                int lst = lct - distance(originBackward[i], destBackward[i]) - patientSrv[p];
                stopBuffer.add(new Stop(i, originBackward[i], est, lst, 2));
                stopBuffer.add(new Stop(i, destBackward[i], ectBack, lct, 3));
                travelBuffer.add(new int[]{stopBuffer.size() - 2, stopBuffer.size() - 1, 0});
            }
        }

        int nStop = stopBuffer.size();

        //Generating stop for start depot:
        for (int v = 0; v < nVehicle; v++)
            stopBuffer.add(new Stop(-1, vehicleStartDepot[v],
                    Math.max(startTime, vehicleStartWindow[v]),
                    Math.min(endTime, vehicleEndWindow[v]), 0));

        for (int v = 0; v < nVehicle; v++)
            stopBuffer.add(new Stop(-1, vehicleStartDepot[v],
                    Math.max(startTime, vehicleStartWindow[v]),
                    Math.min(endTime, vehicleEndWindow[v]), 3));

        Stop[] sites = (Stop[]) stopBuffer.toArray();
        int[][] travels = (int[][]) travelBuffer.toArray();
        int[] travelPatient = new int[travels.length];
        for (int i = 0; i < travels.length; i++) {
            travelPatient[i] = sites[travels[i][1]].patient;
        }

        int nSite = sites.length;

        int[][] transMatrix = new int[nSite][nSite];
        for (int i = 0; i < nSite; i++) {
            for (int j = 0; j < nSite; j++) {
                transMatrix[i][j] = distance(sites[i].place, sites[j].place);
            }
        }

        /* variables */
        Solver cp = makeSolver();

        BoolVar[] visited = new BoolVarImpl[nPatient];
        for (int i = 0; i < nPatient; i++) {
            visited[i] = makeBoolVar(cp);
        }

        //Stop and travel variables:

        IntVar[] arrival = makeIntVarArray(cp, nSite, startTime, endTime);
        for (int i = 0; i < nSite; i++) {
            arrival[i].removeBelow(sites[i].winStart);
            arrival[i].removeAbove(sites[i].winEnd);
        }
        // dummy var?
        IntVar[] duration = makeIntVarArray(cp, nSite, startTime, endTime);
        for (int i = 0; i < nSite; i++) {
            duration[i] = makeIntVar(cp, sites[i].service(), sites[i].service());
        }

        IntVar[] departure = makeIntVarArray(nSite, i -> {
            return new IntVarViewOffset(arrival[i], duration[i].min());
        });
    }
}
        /*
        IntVar[] siteVehicle = makeIntVarArray(cp, nSite, nVehicle+1);
        for (int i=0; i<nSite; i++) {
            if (sites[i].isDepot()) siteVehicle[i].assign((i-nStop)%nVehicle);
            else {
                for (int v=0; v<nVehicle; v++) {
                    if (vehicleCompatibility[v].)
                }
            }
        }
        val siteVehicle: Array[CPIntVar] = Array.tabulate(nSite)(i => {
        if (sites(i).isDepot) CPIntVar((i - nStop) % 2)
        else CPIntVar((0 to nVehicle).filter(v => v == nVehicle || vehicleCompatibility(v).contains(sites(i).category)))
        })

        val siteLoad: Array[CPIntVar] = sites.map(s => CPIntVar(s.load))

        val travelLoad: Array[CPIntVar] = travelPatient.map(p => CPIntVar(patientLoad(p)))

        //Optional time windows:
        val optionalArrival = Array.tabulate(nVehicle + 1, nSite)((_, s) => CPIntVar(arrival(s).min, arrival(s).max))
        val optionalDuration = Array.tabulate(nVehicle + 1, nSite)((_, s) => CPIntVar(duration(s).min, duration(s).max))
        val optionalDeparture = Array.tabulate(nVehicle + 1, nSite)((v, s) => optionalArrival(v)(s) + optionalDuration(v)(s))

        //Sequence variables:
        val sequences: Array[CPInsertSeqVar] = Array.fill(nVehicle)(CPInsertSeqVar(sites.length))


    }
}







        /* constraints */
        /*
        //Setting depots:
        for (v <- sequences.indices) {
        add(First(sequences(v), nStop + v))
        add(Last(sequences(v), nStop + nVehicle + v))
        }

        //Sequence Allocation
        add(SequenceAllocation(sequences.asInstanceOf[Array[CPSeqVar]], sites.indices, siteVehicle))

        //Linking optional and real time windows:
        for (s <- sites.indices) {
        add(AlternativeActivities(
        arrival(s),
        duration(s),
        departure(s),
        optionalArrival.map(_ (s)),
        optionalDuration.map(_ (s)),
        optionalDeparture.map(_ (s)),
        siteVehicle(s)
        ))
        }

        //Precedence constraints:
        travels.foreach{case (pickup, drop, _) =>
        for(v <- sequences.indices) add(Precedence(sequences(v), pickup, drop, dependent = true))
        }

        //Adding dial and ride constraints:
        val (startSites, endSites) = travels.map{case (s, e, _) => (s, e)}.unzip
        for (v <- sequences.indices) {
        add(TransitionTimes(sequences(v),optionalArrival(v),optionalDuration(v),optionalDeparture(v),transMatrix)) //Transition times constraints

        //Adding cumulative:
        val maxCapVar = CPIntVar(vehicleCap(v))
        val minCapVar = CPIntVar(0)
        add(Cumulative(sequences(v), startSites, endSites, travelLoad, maxCapVar, minCapVar))
        }

        //Linking stop visit:
        for (s <- 0 until nStop) add(visited(sites(s).patient) === (siteVehicle(s) ?!== nVehicle))


        /* Objective function */
        /*
        val nServed: CPIntVar = sum(visited)
        maximize(nServed)


        /* Solution */
        /*
        onSolution {
        println("Patients serviced: " + nServed.value)
        println("Sequences:\n" + sequences.mkString("\n"))

        println("------------")
        }


        /* Search */
        /*
        def isDecided(stop: Int): Boolean = {
        if(siteVehicle(stop).isBound) siteVehicle(stop).value == nVehicle || sequences(siteVehicle(stop).value).isMember(stop)
        else false
        }

        def nInsertionsPossible(i: Int): Int = siteVehicle(i).map(s => {
        if(s == nVehicle) 0
        else sequences(s).nCurrentInsertionsFor(i)
        }).sum

        //Return possible insertions for stop: (seq, elem, pred)
        def computeInsertions(i: Int): Seq[(Int, Int, Int)] = siteVehicle(i).filterNot(_ == nVehicle).flatMap(seq => {
        sequences(seq).allCurrentInsertionsFor(i).map(pred => (seq, i, pred))
        }).toSeq

        search {
        val undecidedElems = sites.indices.filterNot(isDecided) //Filtering undecided stops
        if(undecidedElems.isEmpty) noAlternative //If all stops decided => solution
        else{
        val selectStop = undecidedElems.minBy(nInsertionsPossible) //Selecting stop with less insertions possible
        val inserts = Random.shuffle(if(siteVehicle(selectStop).hasValue(nVehicle)) computeInsertions(selectStop) :+ (nVehicle, selectStop, -1) else computeInsertions(selectStop))
        if(inserts.isEmpty) branchOne(throw Inconsistency) //If no insertion possible for stop => Inconsitency
        else branchAll(inserts){
        case (seq, elem, pred) => if(seq == nVehicle) add(siteVehicle(elem) === nVehicle) else sequences(seq).insertAfter(elem, pred)
        }
        }
        }

        println(start())

}
*/